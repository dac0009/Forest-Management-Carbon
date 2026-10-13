"""
fernow_stats.py — shared helpers for the Fernow carbon notebooks.

Pure numpy / scipy / pandas so it runs in ArcGIS Pro's default `arcgispro-py3`
environment (no statsmodels needed). Import from a notebook with:

    import sys; sys.path.insert(0, NOTEBOOK_DIR)
    from fernow_stats import *
"""
import numpy as np
import pandas as pd
from scipy import stats

ORDER = ['WS4', 'WS2', 'WS5']                         # REF, DLC, STS
LABELS = {'WS4': 'REF', 'WS2': 'DLC', 'WS5': 'STS'}
COLORS = {'WS4': '#4C6444', 'WS2': '#B0733A', 'WS5': '#3A6EA5'}

# ---------------------------------------------------------------- unit constants
FT2AC_TO_M2HA = 0.229568          # basal area ft2/ac -> m2/ha
LBAC_TO_MGHA = 0.00112085         # pounds/acre -> Mg/ha
TONAC_TO_MGHA = 2.241702          # US short tons/acre -> Mg/ha
AC_TO_HA = 0.4046856
IN_TO_CM = 2.54
C_FRACTION = 0.48                 # biomass -> carbon (project convention)

# ---------------------------------------------------------------- allometry
# Jenkins et al. (2003) For. Sci. 49:12-35, Table 4: biomass(kg) = exp(b0 + b1 * ln(dbh_cm))
JENKINS = {
    'mixed hardwood':               (-2.4800, 2.4835),
    'hard maple/oak/hickory/beech': (-2.0127, 2.4342),
    'soft maple/birch':             (-1.9123, 2.3651),
    'true fir/hemlock':             (-2.5384, 2.4814),
}

# FEF species codes (from 'FEF Timber Data WS2_WS5.xlsx') -> Jenkins group.
# Codes not listed default to 'mixed hardwood' (~3% of census basal area; reported in notebook).
FEF_TO_JENKINS = {
    11: 'true fir/hemlock',
    21: 'hard maple/oak/hickory/beech', 32: 'hard maple/oak/hickory/beech',
    34: 'hard maple/oak/hickory/beech', 38: 'hard maple/oak/hickory/beech',
    40: 'hard maple/oak/hickory/beech', 41: 'hard maple/oak/hickory/beech',
    47: 'hard maple/oak/hickory/beech', 76: 'hard maple/oak/hickory/beech',
    29: 'soft maple/birch', 30: 'soft maple/birch', 79: 'soft maple/birch',
    28: 'mixed hardwood', 53: 'mixed hardwood', 54: 'mixed hardwood',
    56: 'mixed hardwood', 66: 'mixed hardwood', 84: 'mixed hardwood',
    85: 'mixed hardwood', 88: 'mixed hardwood', 91: 'mixed hardwood',
}


def jenkins_agb_kg(dbh_cm, group):
    b0, b1 = JENKINS[group]
    return np.exp(b0 + b1 * np.log(np.asarray(dbh_cm, float)))


# ---------------------------------------------------------------- group tests
def welch_anova(groups):
    """Welch's one-way ANOVA (R: oneway.test(var.equal = FALSE))."""
    k = len(groups)
    n = np.array([len(g) for g in groups], float)
    m = np.array([g.mean() for g in groups])
    v = np.array([g.var(ddof=1) for g in groups])
    w = n / v
    grand = np.sum(w * m) / np.sum(w)
    num = np.sum(w * (m - grand) ** 2) / (k - 1)
    tmp = np.sum((1 - w / np.sum(w)) ** 2 / (n - 1))
    F = num / (1 + 2 * (k - 2) / (k ** 2 - 1) * tmp)
    df2 = (k ** 2 - 1) / (3 * tmp)
    return F, float(stats.f.sf(F, k - 1, df2))


def dunnett_vs_control(groups):
    """Dunnett vs groups[0]. SciPy >= 1.11 native; else Welch-t + Holm fallback (labelled)."""
    try:
        from scipy.stats import dunnett as _d
        return list(_d(*groups[1:], control=groups[0]).pvalue), 'Dunnett'
    except ImportError:
        p = np.array([stats.ttest_ind(g, groups[0], equal_var=False).pvalue for g in groups[1:]])
        return list(holm(p)), 'Welch-t+Holm (fallback)'


def holm(p):
    p = np.asarray(p, float); m = len(p); order = np.argsort(p)
    adj = np.empty(m); running = 0.0
    for rank, i in enumerate(order):
        running = max(running, (m - rank) * p[i]); adj[i] = min(running, 1.0)
    return adj


def bh_fdr(p):
    """Benjamini-Hochberg adjusted p-values (NaNs preserved)."""
    p = np.asarray(p, float); out = np.full_like(p, np.nan)
    ok = ~np.isnan(p); q = p[ok]; m = len(q)
    if m == 0:
        return out
    order = np.argsort(q); ranked = q[order] * m / (np.arange(m) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    tmp = np.empty(m); tmp[order] = np.minimum(ranked, 1.0); out[ok] = tmp
    return out


def anova_power(groups, alpha=0.05):
    m = [g.mean() for g in groups]
    pooled = np.mean([g.var(ddof=1) for g in groups])
    f_eff = np.sqrt(np.mean([(x - np.mean(m)) ** 2 for x in m]) / pooled)
    N = sum(len(g) for g in groups); k = len(groups)
    fcrit = stats.f.ppf(1 - alpha, k - 1, N - k)
    return f_eff, float(stats.ncf.sf(fcrit, k - 1, N - k, f_eff ** 2 * N))


def battery(groups):
    """Assumption-checked 3-group battery; groups[0] = REF. NaNs dropped."""
    groups = [np.asarray(g, float) for g in groups]
    groups = [g[~np.isnan(g)] for g in groups]
    F, p = stats.f_oneway(*groups)
    lev_p = stats.levene(*groups).pvalue
    sh_p = stats.shapiro(np.concatenate([g - g.mean() for g in groups])).pvalue
    H, pk = stats.kruskal(*groups)
    Fw, pw = welch_anova(groups)
    dp, dm = dunnett_vs_control(groups)
    f_eff, power = anova_power(groups)
    return dict(n=[len(g) for g in groups], means=[g.mean() for g in groups],
                sds=[g.std(ddof=1) for g in groups], F=F, p=p, lev_p=lev_p, sh_p=sh_p,
                Fw=Fw, pw=pw, H=H, pk=pk, dunnett_p=dp, dunnett_method=dm,
                f_eff=f_eff, power=power)


def battery_row(label, groups, **extra):
    r = battery(groups)
    # Pre-specified lead test (never 'the smaller p'): non-normal residuals -> Kruskal-Wallis;
    # heteroscedastic only -> Welch; otherwise classic ANOVA.
    if r['sh_p'] < 0.05:
        lead, lead_p = 'KW', r['pk']
    elif r['lev_p'] < 0.05:
        lead, lead_p = 'Welch', r['pw']
    else:
        lead, lead_p = 'ANOVA', r['p']
    row = dict(test=label, n=str(r['n']),
               REF=r['means'][0], DLC=r['means'][1], STS=r['means'][2],
               ANOVA_p=r['p'], Welch_p=r['pw'], KW_p=r['pk'],
               Levene_p=r['lev_p'], Shapiro_p=r['sh_p'],
               DLCvREF_p=r['dunnett_p'][0], STSvREF_p=r['dunnett_p'][1],
               lead_test=lead, lead_p=lead_p, power=r['power'])
    row.update(extra)
    return row


# ---------------------------------------------------------------- regression
def ols(y, X, names=None, hc3=True):
    """OLS with intercept added. Returns dict with coef table (classical + HC3 SE), R2, AIC, RSS."""
    y = np.asarray(y, float)
    X = np.column_stack([np.ones(len(y))] + [np.asarray(x, float) for x in X])
    n, k = X.shape
    XtX_inv = np.linalg.pinv(X.T @ X)
    b = XtX_inv @ X.T @ y
    e = y - X @ b
    rss = float(e @ e); df = n - k
    se = np.sqrt(np.diag(XtX_inv) * rss / df)
    if hc3:
        h = np.einsum('ij,jk,ik->i', X, XtX_inv, X)
        meat = X.T @ (X * (e ** 2 / (1 - h) ** 2)[:, None])
        se_r = np.sqrt(np.diag(XtX_inv @ meat @ XtX_inv))
    else:
        se_r = se
    t = b / se_r
    p = 2 * stats.t.sf(np.abs(t), df)
    tss = float(((y - y.mean()) ** 2).sum())
    names = ['Intercept'] + (names or [f'x{i}' for i in range(1, k)])
    tab = pd.DataFrame({'coef': b, 'SE': se, 'SE_HC3': se_r, 't_HC3': t, 'p_HC3': p}, index=names)
    aic = n * np.log(rss / n) + 2 * k
    return dict(table=tab, rss=rss, df=df, k=k, n=n, r2=1 - rss / tss,
                adj_r2=1 - (rss / df) / (tss / (n - 1)), aic=aic, resid=e, fitted=X @ b)


def nested_F(full, reduced):
    """Partial F-test comparing nested OLS fits from ols()."""
    q = full['k'] - reduced['k']
    F = ((reduced['rss'] - full['rss']) / q) / (full['rss'] / full['df'])
    return F, float(stats.f.sf(F, q, full['df']))


def freedman_lane(y, cov, factor_dummies, n_perm=9999, seed=1):
    """Permutation p for the factor (watershed) effect after adjusting for covariate(s).
    Freedman & Lane (1983): permute residuals of the reduced (covariate-only) model."""
    rng = np.random.default_rng(seed)
    red = ols(y, cov, hc3=False); full = ols(y, cov + factor_dummies, hc3=False)
    F_obs, _ = nested_F(full, red)
    fitted, resid = red['fitted'], red['resid']
    count = 0
    for _ in range(n_perm):
        y_star = fitted + rng.permutation(resid)
        F_star, _ = nested_F(ols(y_star, cov + factor_dummies, hc3=False),
                             ols(y_star, cov, hc3=False))
        count += F_star >= F_obs
    return F_obs, (count + 1) / (n_perm + 1)


def within_group_corr(df, x, y, group='watershed', method='pearson'):
    """Correlation of group-demeaned x and y; df reduced by number of groups."""
    d = df[[x, y, group]].dropna()
    dx = d[x] - d.groupby(group)[x].transform('mean')
    dy = d[y] - d.groupby(group)[y].transform('mean')
    if method == 'spearman':
        r = stats.spearmanr(dx, dy).correlation
    else:
        r = np.corrcoef(dx, dy)[0, 1]
    dfree = len(d) - d[group].nunique() - 1
    t = r * np.sqrt(dfree / max(1e-12, 1 - r ** 2))
    return r, float(2 * stats.t.sf(abs(t), dfree)), dfree


# ---------------------------------------------------------------- two-sample / bootstrap
def boot_ci(x, stat=np.median, n=5000, seed=1, alpha=0.05):
    rng = np.random.default_rng(seed); x = np.asarray(x, float)
    bs = np.array([stat(rng.choice(x, len(x))) for _ in range(n)])
    return np.quantile(bs, [alpha / 2, 1 - alpha / 2])


def boot_diff_ci(a, b, stat=np.median, n=5000, seed=1, alpha=0.05):
    rng = np.random.default_rng(seed); a = np.asarray(a, float); b = np.asarray(b, float)
    bs = np.array([stat(rng.choice(a, len(a))) - stat(rng.choice(b, len(b))) for _ in range(n)])
    return np.quantile(bs, [alpha / 2, 1 - alpha / 2])


def cliffs_delta(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    gt = sum((x > b).sum() for x in a); lt = sum((x < b).sum() for x in a)
    return (gt - lt) / (len(a) * len(b))


def fmt_p(p):
    return 'NA' if p is None or (isinstance(p, float) and np.isnan(p)) else ('<0.001' if p < 0.001 else f'{p:.3f}')
