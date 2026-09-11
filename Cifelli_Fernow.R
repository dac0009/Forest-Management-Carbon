#R version 4.4.2 (2024-10-31) -- "Pile of Leaves"

setwd("/Users/emelkangi/Desk_Docs/Research/Cifelli_Fernow")
set.seed(123)

#### Packages ####

library(tidyverse)
library(stringr)
library(ggplot2)
library(multcompView)
library(cowplot)

#### Functions and Themes ####

SE.1 <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  (sqrt(var(x)/length(x)))
}

Mean.SE <- list(mean, SE.1)
names(Mean.SE) <- c("Mean", "SE")
rm(SE.1)

dcPalette3 <- c("#FFFFFF","#808080", "#000000")

#### Elements ####

dc <- read.csv("Cifelli_Fernow.csv")

dc$total_PON_gm2 <- dc$LF_N_gm2+dc$hPOM_N_gm2
dc$MAON.tot_PON <- dc$MAOM_N_gm2/dc$total_PON_gm2

str(dc)

dc <- dc %>% 
  mutate(across(c(watershed, plot, subplot, depth, horizon), as.factor))

dc$depth <- ordered(dc$depth, levels = c("OH", "0 to 20", "20 to 40"))

dc$watershed <- ordered(dc$watershed, 
                        levels = c("Reference", "Diameter Limit", 
                                   "Single Tree Selection"))

dc.sum <- dc %>% #for ggplot2
  group_by(watershed, depth) %>%
  reframe(
    across(fine_root_gm2:soil_CN, Mean.SE, .unpack = TRUE, na.rm = TRUE)
  )

dc.frac.sum <- dc %>%
  subset(!(depth %in% c("OH"))) %>%
  group_by(watershed, depth) %>%
  reframe(
    across(LF_C_gcm3:MAON.tot_PON, Mean.SE, .unpack = TRUE, na.rm = TRUE)
  )

dc.OH <- subset(dc, depth %in% c("OH"))
dc.20 <- subset(dc, depth %in% c("0 to 20"))
dc.40 <- subset(dc, depth %in% c("20 to 40"))

dc.frac <- dc %>%
  subset(horizon %in% c("Mineral")) %>%
  subset(str_detect(subplot, "A"))

dc.frac.20 <- subset(dc.frac, depth %in% c("0 to 20"))
dc.frac.40 <- subset(dc.frac, depth %in% c("20 to 40"))

#### Statistics ####

options(na.action = na.omit)

roots.OH.aov <- aov(fine_root_gm2 ~ watershed, dc[(dc$depth %in% c("OH")),])
summary(roots.OH.aov) #F=2.601, p=0.0844

roots.OH.TK <- TukeyHSD(roots.OH.aov)
roots.OH.letters <- multcompLetters4(roots.OH.aov, roots.OH.TK)
roots.OH.letters$'watershed'$Letters #ref = a, sts = a, dl = a

a <- as.data.frame(roots.OH.letters$'watershed'$Letters)
a$depth <- rep(c("OH"), 3)
a <- rownames_to_column(a)
colnames(a) <- c("watershed", "roots.let", "depth")

roots.20.aov <- aov(fine_root_gm2 ~ watershed, dc[(dc$depth %in% c("0 to 20")),])
summary(roots.20.aov) #n.s.

roots.20.TK <- TukeyHSD(roots.20.aov)
roots.20.letters <- multcompLetters4(roots.20.aov, roots.20.TK)
roots.20.letters$'watershed'$Letters #ref = a, sts = a, dl = a

b <- as.data.frame(roots.20.letters$'watershed'$Letters)
b$depth <- rep(c("0 to 20"), 3)
b <- rownames_to_column(b)
colnames(b) <- c("watershed", "roots.let", "depth")

roots.40.aov <- aov(fine_root_gm2 ~ watershed, dc[(dc$depth %in% c("20 to 40")),])
summary(roots.40.aov) #F=3.818, p=0.0286

roots.40.TK <- TukeyHSD(roots.40.aov)
roots.40.letters <- multcompLetters4(roots.40.aov, roots.40.TK)
roots.40.letters$'watershed'$Letters #ref = a, sts = a, dl = a

c <- as.data.frame(roots.40.letters$'watershed'$Letters)
c$depth <- rep(c("20 to 40"), 3)
c <- rownames_to_column(c)
colnames(c) <- c("watershed", "roots.let", "depth")

roots.let <- bind_rows(a, b, c)
dc.sum <- left_join(dc.sum, roots.let, by = c("watershed", "depth"))

rm(roots.OH.aov, roots.20.aov, roots.40.aov,
   roots.OH.TK, roots.20.TK, roots.40.TK,
   roots.OH.letters, roots.20.letters, roots.40.letters,
   a, b, c, roots.let)

#soil C per m2 (area of sampling) (soilC_gm2)

soilC.OH.aov <- aov(soil_C_gm2 ~ watershed, dc[(dc$depth %in% c("OH")),]) #w/o plot as block
summary(soilC.OH.aov) #F=4.537, p=0.0154

soilC.OH.TK <- TukeyHSD(soilC.OH.aov)
soilC.OH.letters <- multcompLetters4(soilC.OH.aov, soilC.OH.TK)
soilC.OH.letters$'watershed'$Letters #ref = a, sts = ab, dl = b

a <- as.data.frame(soilC.OH.letters$'watershed'$Letters)
a$depth <- rep(c("OH"), 3)
a <- rownames_to_column(a)
colnames(a) <- c("watershed", "soilC.let", "depth")

soilC.20.aov <- aov(soil_C_gm2 ~ watershed, dc.20) #w/o plot as block
summary(soilC.20.aov) #F=2.531, p=0.0895

soilC.20.TK <- TukeyHSD(soilC.20.aov)
soilC.20.letters <- multcompLetters4(soilC.20.aov, soilC.20.TK)
soilC.20.letters$'watershed'$Letters #ref = a, sts = a, dl = a

b <- as.data.frame(soilC.20.letters$'watershed'$Letters)
b$depth <- rep(c("0 to 20"), 3)
b <- rownames_to_column(b)
colnames(b) <- c("watershed", "soilC.let", "depth")

soilC.40.aov <- aov(soil_C_gm2 ~ watershed, dc.40) #w/o plot as block
summary(soilC.40.aov) #n.s. (F=1.75, p=0.184)

soilC.40.TK <- TukeyHSD(soilC.40.aov)
soilC.40.letters <- multcompLetters4(soilC.40.aov, soilC.40.TK)
soilC.40.letters$'watershed'$Letters #ref = a, sts = a, dl = a

c <- as.data.frame(soilC.40.letters$'watershed'$Letters)
c$depth <- rep(c("20 to 40"), 3)
c <- rownames_to_column(c)
colnames(c) <- c("watershed", "soilC.let", "depth")

soilC.let <- bind_rows(a, b, c)
dc.sum <- left_join(dc.sum, soilC.let, by = c("watershed", "depth"))

rm(soilC.OH.aov, soilC.20.aov, soilC.40.aov,
   soilC.OH.TK, soilC.20.TK, soilC.40.TK,
   soilC.OH.letters, soilC.20.letters, soilC.40.letters,
   a, b, c, soilC.let)

#soil N per m2 (area of sampling) (soilN_gm2)

soilN.OH.aov <- aov(soil_N_gm2 ~ watershed, dc.OH) #w/o plot as block
summary(soilN.OH.aov) #F=4.821, p=0.0121

soilN.OH.TK <- TukeyHSD(soilN.OH.aov)
soilN.OH.letters <- multcompLetters4(soilN.OH.aov, soilN.OH.TK)
soilN.OH.letters$'watershed'$Letters #ref = a, sts = ab, dl = b

a <- as.data.frame(soilN.OH.letters$'watershed'$Letters)
a$depth <- rep(c("OH"), 3)
a <- rownames_to_column(a)
colnames(a) <- c("watershed", "soilN.let", "depth")

soilN.20.aov <- aov(soil_N_gm2 ~ watershed, dc.20) #w/o plot as block
summary(soilN.20.aov) #F=3.683, p=0.032

soilN.20.TK <- TukeyHSD(soilN.20.aov)
soilN.20.letters <- multcompLetters4(soilN.20.aov, soilN.20.TK)
soilN.20.letters$'watershed'$Letters #ref = a, sts = b, dl = ab

b <- as.data.frame(soilN.20.letters$'watershed'$Letters)
b$depth <- rep(c("0 to 20"), 3)
b <- rownames_to_column(b)
colnames(b) <- c("watershed", "soilN.let", "depth")

soilN.40.aov <- aov(soil_N_gm2 ~ watershed, dc.40) #w/o plot as block
summary(soilN.40.aov) #n.s.

soilN.40.TK <- TukeyHSD(soilN.40.aov)
soilN.40.letters <- multcompLetters4(soilN.40.aov, soilN.40.TK)
soilN.40.letters$'watershed'$Letters #ref = a, sts = a, dl = a

c <- as.data.frame(soilN.40.letters$'watershed'$Letters)
c$depth <- rep(c("20 to 40"), 3)
c <- rownames_to_column(c)
colnames(c) <- c("watershed", "soilN.let", "depth")

soilN.let <- bind_rows(a, b, c)
dc.sum <- left_join(dc.sum, soilN.let, by = c("watershed", "depth"))

rm(soilN.OH.aov, soilN.20.aov, soilN.40.aov,
   soilN.OH.TK, soilN.20.TK, soilN.40.TK,
   soilN.OH.letters, soilN.20.letters, soilN.40.letters,
   a, b, c, soilN.let)

#soil C-to-N ratio (Soil_CN)

soilCN.OH.aov <- aov(soil_CN ~ watershed, dc.OH) #w/o plot as block
summary(soilCN.OH.aov) #n.s.

soilCN.OH.TK <- TukeyHSD(soilCN.OH.aov)
soilCN.OH.letters <- multcompLetters4(soilCN.OH.aov, soilCN.OH.TK)
soilCN.OH.letters$'watershed'$Letters #ref = a, sts = a, dl = a

a <- as.data.frame(soilCN.OH.letters$'watershed'$Letters)
a$depth <- rep(c("OH"), 3)
a <- rownames_to_column(a)
colnames(a) <- c("watershed", "soilCN.let", "depth")

soilCN.20.aov <- aov(soil_CN ~ watershed, dc.20) #w/o plot as block
summary(soilCN.20.aov) #n.s.

soilCN.20.TK <- TukeyHSD(soilCN.20.aov)
soilCN.20.letters <- multcompLetters4(soilCN.20.aov, soilCN.20.TK)
soilCN.20.letters$'watershed'$Letters #ref = a, sts = a, dl = a

b <- as.data.frame(soilCN.20.letters$'watershed'$Letters)
b$depth <- rep(c("0 to 20"), 3)
b <- rownames_to_column(b)
colnames(b) <- c("watershed", "soilCN.let", "depth")

soilCN.40.aov <- aov(soil_CN ~ watershed, dc.40) #w/o plot as block
summary(soilCN.40.aov) #F=4.939, p=0.0109

soilCN.40.TK <- TukeyHSD(soilCN.40.aov)
soilCN.40.letters <- multcompLetters4(soilCN.40.aov, soilCN.40.TK)
soilCN.40.letters$'watershed'$Letters #ref = ab, sts = b, dl = a

c <- as.data.frame(soilCN.40.letters$'watershed'$Letters)
c$depth <- rep(c("20 to 40"), 3)
c <- rownames_to_column(c)
colnames(c) <- c("watershed", "soilCN.let", "depth")

soilCN.let <- bind_rows(a, b, c)
dc.sum <- left_join(dc.sum, soilCN.let, by = c("watershed", "depth"))

rm(soilCN.OH.aov, soilCN.20.aov, soilCN.40.aov,
   soilCN.OH.TK, soilCN.20.TK, soilCN.40.TK,
   soilCN.OH.letters, soilCN.20.letters, soilCN.40.letters,
   a, b, c, soilCN.let)

#Fractions

#Light Fraction: g of C per area soil (g·m-2)

LF.C.aov <- aov(LF_C_gm2 ~ watershed * depth, dc.frac)
summary(LF.C.aov) #depth = 0.0148, ws effect n.s. (may be masked in depth)

LF.C.20.aov <- aov(LF_C_gm2 ~ watershed, dc.frac[(dc.frac$depth %in% c("0 to 20")),])
summary(LF.C.20.aov) #n.s.

LF.C.20.TK <- TukeyHSD(LF.C.20.aov)
LF.C.20.letters <- multcompLetters4(LF.C.20.aov, LF.C.20.TK)
LF.C.20.letters$'watershed'$Letters #ref = a, sts = a, dl = a

a <- as.data.frame(LF.C.20.letters$'watershed'$Letters)
a$depth <- rep(c("0 to 20"), 3)
a <- rownames_to_column(a)
colnames(a) <- c("watershed", "LFC.let", "depth")

LF.C.40.aov <- aov(LF_C_gm2 ~ watershed, dc.frac[(dc.frac$depth %in% c("20 to 40")),])
summary(LF.C.40.aov) #F = 6.415, p = 0.0127

LF.C.40.TK <- TukeyHSD(LF.C.40.aov)
LF.C.40.letters <- multcompLetters4(LF.C.40.aov, LF.C.40.TK)
LF.C.40.letters$'watershed'$Letters #ref = b, sts = ab, dl = a

b <- as.data.frame(LF.C.40.letters$'watershed'$Letters)
b$depth <- rep(c("20 to 40"), 3)
b <- rownames_to_column(b)
colnames(b) <- c("watershed", "LFC.let", "depth")

LFC.let <- bind_rows(a, b)
dc.frac.sum <- left_join(dc.frac.sum, LFC.let, by = c("watershed", "depth"))

rm(LF.C.aov, LF.C.20.aov, LF.C.40.aov,
   LF.C.20.TK, LF.C.40.TK,
   LF.C.20.letters, LF.C.40.letters,
   a, b, LFC.let)

#Light Fraction: g of N per area soil (g·m-2)

LF.N.aov <- aov(LF_N_gm2 ~ watershed * depth, dc.frac)
summary(LF.N.aov) #depth = 0.0118, ws effect n.s.

LF.N.20.aov <- aov(LF_N_gm2 ~ watershed, dc.frac[(dc.frac$depth %in% c("0 to 20")),])
summary(LF.N.20.aov) #n.s.

LF.N.20.TK <- TukeyHSD(LF.N.20.aov)
LF.N.20.letters <- multcompLetters4(LF.N.20.aov, LF.N.20.TK)
LF.N.20.letters$'watershed'$Letters #ref = a, sts = a, dl = a

a <- as.data.frame(LF.N.20.letters$'watershed'$Letters)
a$depth <- rep(c("0 to 20"), 3)
a <- rownames_to_column(a)
colnames(a) <- c("watershed", "LFN.let", "depth")

LF.N.40.aov <- aov(LF_N_gm2 ~ watershed, dc.frac[(dc.frac$depth %in% c("20 to 40")),])
summary(LF.N.40.aov) #F = 8.284, p = 0.0055

LF.N.40.TK <- TukeyHSD(LF.N.40.aov)
LF.N.40.letters <- multcompLetters4(LF.N.40.aov, LF.N.40.TK)
LF.N.40.letters$'watershed'$Letters #ref = b, sts = ab, dl = a

b <- as.data.frame(LF.N.40.letters$'watershed'$Letters)
b$depth <- rep(c("20 to 40"), 3)
b <- rownames_to_column(b)
colnames(b) <- c("watershed", "LFN.let", "depth")

LFN.let <- bind_rows(a, b)
dc.frac.sum <- left_join(dc.frac.sum, LFN.let, by = c("watershed", "depth"))

rm(LF.N.aov, LF.N.20.aov, LF.N.40.aov,
   LF.N.20.TK, LF.N.40.TK,
   LF.N.20.letters, LF.N.40.letters,
   a, b, LFN.let)

#heavy POM: g of C per area soil (g·m-2)

hPOM.C.aov <- aov(hPOM_C_gm2 ~ watershed * depth, dc.frac)
summary(hPOM.C.aov) #depth = 0.0078, ws effect n.s.

hPOM.C.20.aov <- aov(hPOM_C_gm2 ~ watershed, dc.frac[(dc.frac$depth %in% c("0 to 20")),])
summary(hPOM.C.20.aov) #n.s.

hPOM.C.20.TK <- TukeyHSD(hPOM.C.20.aov)
hPOM.C.20.letters <- multcompLetters4(hPOM.C.20.aov, hPOM.C.20.TK)
hPOM.C.20.letters$'watershed'$Letters #ref = a, sts = a, dl = a

a <- as.data.frame(hPOM.C.20.letters$'watershed'$Letters)
a$depth <- rep(c("0 to 20"), 3)
a <- rownames_to_column(a)
colnames(a) <- c("watershed", "hPOC.let", "depth")

hPOM.C.40.aov <- aov(hPOM_C_gm2 ~ watershed, dc.frac[(dc.frac$depth %in% c("20 to 40")),])
summary(hPOM.C.40.aov) #n.s.

hPOM.C.40.TK <- TukeyHSD(hPOM.C.40.aov)
hPOM.C.40.letters <- multcompLetters4(hPOM.C.40.aov, hPOM.C.40.TK)
hPOM.C.40.letters$'watershed'$Letters #ref = a, sts = a, dl = a

b <- as.data.frame(hPOM.C.40.letters$'watershed'$Letters)
b$depth <- rep(c("20 to 40"), 3)
b <- rownames_to_column(b)
colnames(b) <- c("watershed", "hPOC.let", "depth")

hPOC.let <- bind_rows(a, b)
dc.frac.sum <- left_join(dc.frac.sum, hPOC.let, by = c("watershed", "depth"))

rm(hPOM.C.aov, hPOM.C.20.aov, hPOM.C.40.aov,
   hPOM.C.20.TK, hPOM.C.40.TK,
   hPOM.C.20.letters, hPOM.C.40.letters,
   a, b, hPOC.let)

#heavy POM: g of N per area soil (g·m-2)

hPOM.N.aov <- aov(hPOM_N_gm2 ~ watershed * depth, dc.frac)
summary(hPOM.N.aov) #depth = 0.0046, ws effect n.s.

hPOM.N.20.aov <- aov(hPOM_N_gm2 ~ watershed, dc.frac[(dc.frac$depth %in% c("0 to 20")),])
summary(hPOM.N.20.aov) #n.s.

hPOM.N.20.TK <- TukeyHSD(hPOM.N.20.aov)
hPOM.N.20.letters <- multcompLetters4(hPOM.N.20.aov, hPOM.N.20.TK)
hPOM.N.20.letters$'watershed'$Letters #ref = a, sts = a, dl = a

a <- as.data.frame(hPOM.N.20.letters$'watershed'$Letters)
a$depth <- rep(c("0 to 20"), 3)
a <- rownames_to_column(a)
colnames(a) <- c("watershed", "hPON.let", "depth")

hPOM.N.40.aov <- aov(hPOM_N_gm2 ~ watershed, dc.frac[(dc.frac$depth %in% c("20 to 40")),])
summary(hPOM.N.40.aov) #n.s.

hPOM.N.40.TK <- TukeyHSD(hPOM.N.40.aov)
hPOM.N.40.letters <- multcompLetters4(hPOM.N.40.aov, hPOM.N.40.TK)
hPOM.N.40.letters$'watershed'$Letters #ref = a, sts = a, dl = a

b <- as.data.frame(hPOM.N.40.letters$'watershed'$Letters)
b$depth <- rep(c("20 to 40"), 3)
b <- rownames_to_column(b)
colnames(b) <- c("watershed", "hPON.let", "depth")

hPON.let <- bind_rows(a, b)
dc.frac.sum <- left_join(dc.frac.sum, hPON.let, by = c("watershed", "depth"))

rm(hPOM.N.aov, hPOM.N.20.aov, hPOM.N.40.aov,
   hPOM.N.20.TK, hPOM.N.40.TK,
   hPOM.N.20.letters, hPOM.N.40.letters,
   a, b, hPON.let)

#MAOM: g of C per area soil (g·m-2)

MAOM.C.aov <- aov(MAOM_C_gm2 ~ watershed * depth, dc.frac)
summary(MAOM.C.aov) #depth = 0.0003, ws effect n.s.

MAOM.C.20.aov <- aov(MAOM_C_gm2 ~ watershed, dc.frac[(dc.frac$depth %in% c("0 to 20")),])
summary(MAOM.C.20.aov) #n.s.

MAOM.C.20.TK <- TukeyHSD(MAOM.C.20.aov)
MAOM.C.20.letters <- multcompLetters4(MAOM.C.20.aov, MAOM.C.20.TK)
MAOM.C.20.letters$'watershed'$Letters #ref = a, sts = a, dl = a

a <- as.data.frame(MAOM.C.20.letters$'watershed'$Letters)
a$depth <- rep(c("0 to 20"), 3)
a <- rownames_to_column(a)
colnames(a) <- c("watershed", "MAOC.let", "depth")

MAOM.C.40.aov <- aov(MAOM_C_gm2 ~ watershed, dc.frac[(dc.frac$depth %in% c("20 to 40")),])
summary(MAOM.C.40.aov) #n.s.

MAOM.C.40.TK <- TukeyHSD(MAOM.C.40.aov)
MAOM.C.40.letters <- multcompLetters4(MAOM.C.40.aov, MAOM.C.40.TK)
MAOM.C.40.letters$'watershed'$Letters #ref = a, sts = a, dl = a

b <- as.data.frame(MAOM.C.40.letters$'watershed'$Letters)
b$depth <- rep(c("20 to 40"), 3)
b <- rownames_to_column(b)
colnames(b) <- c("watershed", "MAOC.let", "depth")

MAOC.let <- bind_rows(a, b)
dc.frac.sum <- left_join(dc.frac.sum, MAOC.let, by = c("watershed", "depth"))

rm(MAOM.C.aov, MAOM.C.20.aov, MAOM.C.40.aov,
   MAOM.C.20.TK, MAOM.C.40.TK,
   MAOM.C.20.letters, MAOM.C.40.letters,
   a, b, MAOC.let)

#MAOM: g of N per area soil (g·m-2)

MAOM.N.aov <- aov(MAOM_N_gm2 ~ watershed * depth, dc.frac)
summary(MAOM.N.aov) #depth = 0.0032, ws effect n.s.

MAOM.N.20.aov <- aov(MAOM_N_gm2 ~ watershed, dc.frac[(dc.frac$depth %in% c("0 to 20")),])
summary(MAOM.N.20.aov) #n.s.

MAOM.N.20.TK <- TukeyHSD(MAOM.N.20.aov)
MAOM.N.20.letters <- multcompLetters4(MAOM.N.20.aov, MAOM.N.20.TK)
MAOM.N.20.letters$'watershed'$Letters #ref = a, sts = a, dl = a

a <- as.data.frame(MAOM.N.20.letters$'watershed'$Letters)
a$depth <- rep(c("0 to 20"), 3)
a <- rownames_to_column(a)
colnames(a) <- c("watershed", "MAON.let", "depth")

MAOM.N.40.aov <- aov(MAOM_N_gm2 ~ watershed, dc.frac[(dc.frac$depth %in% c("20 to 40")),])
summary(MAOM.N.40.aov) #n.s.

MAOM.N.40.TK <- TukeyHSD(MAOM.N.40.aov)
MAOM.N.40.letters <- multcompLetters4(MAOM.N.40.aov, MAOM.N.40.TK)
MAOM.N.40.letters$'watershed'$Letters #ref = a, sts = a, dl = a

b <- as.data.frame(MAOM.N.40.letters$'watershed'$Letters)
b$depth <- rep(c("20 to 40"), 3)
b <- rownames_to_column(b)
colnames(b) <- c("watershed", "MAON.let", "depth")

MAON.let <- bind_rows(a, b)
dc.frac.sum <- left_join(dc.frac.sum, MAON.let, by = c("watershed", "depth"))

rm(MAOM.N.aov, MAOM.N.20.aov, MAOM.N.40.aov,
   MAOM.N.20.TK, MAOM.N.40.TK,
   MAOM.N.20.letters, MAOM.N.40.letters,
   a, b, MAON.let)

#total POC: g of C in LF and heavy POM per area soil (g·m-2)

tot.POC.aov <- aov(total_POC_gm2 ~ watershed * depth, dc.frac)
summary(tot.POC.aov) #depth = 0.0003, ws effect n.s.

tot.POC.20.aov <- aov(total_POC_gm2 ~ watershed, dc.frac[(dc.frac$depth %in% c("0 to 20")),])
summary(tot.POC.20.aov) #n.s.

tot.POC.20.TK <- TukeyHSD(tot.POC.20.aov)
tot.POC.20.letters <- multcompLetters4(tot.POC.20.aov, tot.POC.20.TK)
tot.POC.20.letters$'watershed'$Letters #ref = a, sts = a, dl = a

a <- as.data.frame(tot.POC.20.letters$'watershed'$Letters)
a$depth <- rep(c("0 to 20"), 3)
a <- rownames_to_column(a)
colnames(a) <- c("watershed", "tPOC.let", "depth")

tot.POC.40.aov <- aov(total_POC_gm2 ~ watershed, dc.frac[(dc.frac$depth %in% c("20 to 40")),])
summary(tot.POC.40.aov) #F = 4.397, p = 0.0314

tot.POC.40.TK <- TukeyHSD(tot.POC.40.aov)
tot.POC.40.letters <- multcompLetters4(tot.POC.40.aov, tot.POC.40.TK)
tot.POC.40.letters$'watershed'$Letters #ref = ab, sts = b, dl = a

b <- as.data.frame(tot.POC.40.letters$'watershed'$Letters)
b$depth <- rep(c("20 to 40"), 3)
b <- rownames_to_column(b)
colnames(b) <- c("watershed", "tPOC.let", "depth")

tPOC.let <- bind_rows(a, b)
dc.frac.sum <- left_join(dc.frac.sum, tPOC.let, by = c("watershed", "depth"))

rm(tot.POC.aov, tot.POC.20.aov, tot.POC.40.aov,
   tot.POC.20.TK, tot.POC.40.TK,
   tot.POC.20.letters, tot.POC.40.letters,
   a, b, tPOC.let)

#Ratio: MAOC-to-heavy POC

MAOC.hPOC.aov <- aov(MAOC.hPOC ~ watershed * depth, dc.frac)
summary(MAOC.hPOC.aov) #n.s.

MAOC.hPOC.20.aov <- aov(MAOC.hPOC ~ watershed, dc.frac[(dc.frac$depth %in% c("0 to 20")),])
summary(MAOC.hPOC.20.aov) #n.s.

MAOC.hPOC.20.TK <- TukeyHSD(MAOC.hPOC.20.aov)
MAOC.hPOC.20.letters <- multcompLetters4(MAOC.hPOC.20.aov, MAOC.hPOC.20.TK)
MAOC.hPOC.20.letters$'watershed'$Letters #ref = a, sts = a, dl = a

a <- as.data.frame(MAOC.hPOC.20.letters$'watershed'$Letters)
a$depth <- rep(c("0 to 20"), 3)
a <- rownames_to_column(a)
colnames(a) <- c("watershed", "MAOC.hPOC.let", "depth")

MAOC.hPOC.40.aov <- aov(MAOC.hPOC ~ watershed, dc.frac[(dc.frac$depth %in% c("20 to 40")),])
summary(MAOC.hPOC.40.aov) #n.s.

MAOC.hPOC.40.TK <- TukeyHSD(MAOC.hPOC.40.aov)
MAOC.hPOC.40.letters <- multcompLetters4(MAOC.hPOC.40.aov, MAOC.hPOC.40.TK)
MAOC.hPOC.40.letters$'watershed'$Letters #ref = a, sts = a, dl = a

b <- as.data.frame(MAOC.hPOC.40.letters$'watershed'$Letters)
b$depth <- rep(c("20 to 40"), 3)
b <- rownames_to_column(b)
colnames(b) <- c("watershed", "MAOC.hPOC.let", "depth")

MAOC.hPOC.let <- bind_rows(a, b)
dc.frac.sum <- left_join(dc.frac.sum, MAOC.hPOC.let, 
                         by = c("watershed", "depth"))

rm(MAOC.hPOC.aov, MAOC.hPOC.20.aov, MAOC.hPOC.40.aov,
   MAOC.hPOC.20.TK, MAOC.hPOC.40.TK,
   MAOC.hPOC.20.letters, MAOC.hPOC.40.letters,
   a, b, MAOC.hPOC.let)

#Ratio: MAOC-to-total POC

MAOC.tPOC.aov <- aov(MAOC.tot_POC ~ watershed * depth, dc.frac)
summary(MAOC.tPOC.aov) #ws*d = 0.0328, ws = 0.0514, depth = 0.0206
#driven by LF

MAOC.tPOC.20.aov <- aov(MAOC.tot_POC ~ watershed, dc.frac[(dc.frac$depth %in% c("0 to 20")),])
summary(MAOC.tPOC.20.aov) #n.s.

MAOC.tPOC.20.TK <- TukeyHSD(MAOC.tPOC.20.aov)
MAOC.tPOC.20.letters <- multcompLetters4(MAOC.tPOC.20.aov, MAOC.tPOC.20.TK)
MAOC.tPOC.20.letters$'watershed'$Letters #ref = a, sts = a, dl = a

a <- as.data.frame(MAOC.tPOC.20.letters$'watershed'$Letters)
a$depth <- rep(c("0 to 20"), 3)
a <- rownames_to_column(a)
colnames(a) <- c("watershed", "MAOC.tPOC.let", "depth")

MAOC.tPOC.40.aov <- aov(MAOC.tot_POC ~ watershed, dc.frac[(dc.frac$depth %in% c("20 to 40")),])
summary(MAOC.tPOC.40.aov) #F = 3.874, p = 0.044

MAOC.tPOC.40.TK <- TukeyHSD(MAOC.tPOC.40.aov)
MAOC.tPOC.40.letters <- multcompLetters4(MAOC.tPOC.40.aov, MAOC.tPOC.40.TK)
MAOC.tPOC.40.letters$'watershed'$Letters #ref = ab, sts = a, dl = b

b <- as.data.frame(MAOC.tPOC.40.letters$'watershed'$Letters)
b$depth <- rep(c("20 to 40"), 3)
b <- rownames_to_column(b)
colnames(b) <- c("watershed", "MAOC.tPOC.let", "depth")

MAOC.tPOC.let <- bind_rows(a, b)
dc.frac.sum <- left_join(dc.frac.sum, MAOC.tPOC.let, 
                         by = c("watershed", "depth"))

rm(MAOC.tPOC.aov, MAOC.tPOC.20.aov, MAOC.tPOC.40.aov,
   MAOC.tPOC.20.TK, MAOC.tPOC.40.TK,
   MAOC.tPOC.20.letters, MAOC.tPOC.40.letters,
   a, b, MAOC.tPOC.let)

#Ratio: MAON-to-total PON

MAON.tPON.aov <- aov(MAON.tot_PON ~ watershed * depth, dc.frac)
summary(MAON.tPON.aov) #ws*d = n.s., ws = n.s., depth = 0.00998

MAON.tPON.20.aov <- aov(MAON.tot_PON ~ watershed, dc.frac[(dc.frac$depth %in% c("0 to 20")),])
summary(MAON.tPON.20.aov) #n.s.

MAON.tPON.40.aov <- aov(MAON.tot_PON ~ watershed, dc.frac[(dc.frac$depth %in% c("20 to 40")),])
summary(MAON.tPON.40.aov) #n.s.

#### Figures ####

theme_dc <- function() {
  theme_classic() +
    theme(legend.position = "none",
          axis.title.y = element_text(size=12),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=12, angle = 45, hjust = 1),
          axis.text.y = element_text(size=12),
          strip.background = element_blank(),
          strip.text.x = element_blank())
}

#Roots & Soil C, N and C:N

rootplot <- ggplot(dc.sum, 
                   aes(x=watershed, y=fine_root_gm2_Mean, fill=watershed))+
  geom_bar(position=position_dodge(), colour="black", stat="identity") + 
  geom_errorbar(aes(ymin=fine_root_gm2_Mean-fine_root_gm2_SE,
                    ymax=fine_root_gm2_Mean+fine_root_gm2_SE), 
                width=0.2, position=position_dodge(0.9)) +
  scale_fill_manual(values = dcPalette3) + 
  labs(y=bquote('Fine Root Biomass'~(g/m^2))) +
  geom_text(aes(label = roots.let, y = fine_root_gm2_Mean+fine_root_gm2_SE), 
            position = position_dodge(0.9), 
            vjust = -1, size = 10, size.unit = "pt") +
  theme_dc() +
  facet_wrap(~depth, nrow = 3, scales = "free_y") +
  coord_cartesian(clip = "off")
rootplot

#ggsave("DC_Root.svg", plot = rootplot, dpi = 300, height = 9, width = 3)

soilCplot <- ggplot(dc.sum, 
                    aes(x=watershed, y=soil_C_gm2_Mean, fill=watershed))+
  geom_bar(position=position_dodge(), colour="black", stat="identity") + 
  geom_errorbar(aes(ymin=soil_C_gm2_Mean-soil_C_gm2_SE,
                    ymax=soil_C_gm2_Mean+soil_C_gm2_SE), 
                width=0.2, position=position_dodge(0.9)) +
  scale_fill_manual(values = dcPalette3) + 
  labs(y=bquote('Soil Carbon'~(g/m^2))) +
  geom_text(aes(label = soilC.let, y = soil_C_gm2_Mean+soil_C_gm2_SE), 
            position = position_dodge(0.9), 
            vjust = -1, size = 10, size.unit = "pt") +
  theme_dc() +
  facet_wrap(~depth, nrow = 3, scales = "free_y") +
  coord_cartesian(clip = "off")
soilCplot

#ggsave("DC_SoilC.svg", plot = soilCplot, dpi = 300, height = 9, width = 3)

soilNplot <- ggplot(dc.sum, 
                    aes(x=watershed, y=soil_N_gm2_Mean, fill=watershed))+
  geom_bar(position=position_dodge(), colour="black", stat="identity") + 
  geom_errorbar(aes(ymin=soil_N_gm2_Mean-soil_N_gm2_SE,
                    ymax=soil_N_gm2_Mean+soil_N_gm2_SE), 
                width=0.2, position=position_dodge(0.9)) +
  scale_fill_manual(values = dcPalette3) + 
  labs(y=bquote('Soil Nitrogen'~(g/m^2))) +
  geom_text(aes(label = soilN.let, y = soil_N_gm2_Mean+soil_N_gm2_SE), 
            position = position_dodge(0.9), 
            vjust = -1, size = 10, size.unit = "pt") +
  theme_dc() +
  facet_wrap(~depth, nrow = 3, scales = "free_y") +
  coord_cartesian(clip = "off")
soilNplot

#ggsave("DC_SoilN.svg", plot = soilNplot, dpi = 300, height = 9, width = 3)

soilCNplot <- ggplot(dc.sum, 
                     aes(x=watershed, y=soil_CN_Mean, fill=watershed))+
  geom_bar(position=position_dodge(), colour="black", stat="identity") + 
  geom_errorbar(aes(ymin=soil_CN_Mean-soil_CN_SE,
                    ymax=soil_CN_Mean+soil_CN_SE), 
                width=0.2, position=position_dodge(0.9)) +
  scale_fill_manual(values = dcPalette3) + 
  labs(y=bquote('Soil C-to-N Ratio')) +
  geom_text(aes(label = soilCN.let, y = soil_CN_Mean+soil_CN_SE), 
            position = position_dodge(0.9), 
            vjust = -1, size = 10, size.unit = "pt") +
  theme_dc() +
  facet_wrap(~depth, nrow = 3, scales = "free_y") +
  coord_cartesian(clip = "off")
soilCNplot

#ggsave("DC_SoilCN.svg", plot = soilCNplot, dpi = 300, height = 9, width = 3)

#Fracs: C, N and MAOC:POC (heavy & total)

LFCplot <- ggplot(dc.frac.sum, 
                  aes(x=watershed, y=LF_C_gm2_Mean, fill=watershed))+
  geom_bar(position=position_dodge(), colour="black", stat="identity") + 
  geom_errorbar(aes(ymin=LF_C_gm2_Mean-LF_C_gm2_SE,
                    ymax=LF_C_gm2_Mean+LF_C_gm2_SE), 
                width=0.2, position=position_dodge(0.9)) +
  scale_fill_manual(values = dcPalette3) + 
  labs(y=bquote('Carbon per Fraction'~(g/m^2))) +
  geom_text(aes(label = LFC.let, y = LF_C_gm2_Mean+LF_C_gm2_SE), 
            position = position_dodge(0.9),
            vjust = -1, size = 10, size.unit = "pt") +
  theme_dc() +
  facet_wrap(~depth, nrow = 2, scales = "free_y") +
  coord_cartesian(clip = "off")
LFCplot

hPOCplot <- ggplot(dc.frac.sum, 
                   aes(x=watershed, y=hPOM_C_gm2_Mean, fill=watershed))+
  geom_bar(position=position_dodge(), colour="black", stat="identity") + 
  geom_errorbar(aes(ymin=hPOM_C_gm2_Mean-hPOM_C_gm2_SE,
                    ymax=hPOM_C_gm2_Mean+hPOM_C_gm2_SE), 
                width=0.2, position=position_dodge(0.9)) +
  scale_fill_manual(values = dcPalette3) + 
  labs(y=bquote('Carbon per Fraction'~(g/m^2))) +
  geom_text(aes(label = hPOC.let, y = hPOM_C_gm2_Mean+hPOM_C_gm2_SE), 
            position = position_dodge(0.9),
            vjust = -1, size = 10, size.unit = "pt") +
  theme_dc() +
  facet_wrap(~depth, nrow = 2, scales = "free_y") +
  coord_cartesian(clip = "off")
hPOCplot

MAOCplot <- ggplot(dc.frac.sum, 
                   aes(x=watershed, y=MAOM_C_gm2_Mean, fill=watershed))+
  geom_bar(position=position_dodge(), colour="black", stat="identity") + 
  geom_errorbar(aes(ymin=MAOM_C_gm2_Mean-MAOM_C_gm2_SE,
                    ymax=MAOM_C_gm2_Mean+MAOM_C_gm2_SE), 
                width=0.2, position=position_dodge(0.9)) +
  scale_fill_manual(values = dcPalette3) + 
  labs(y=bquote('Carbon per Fraction'~(g/m^2))) +
  geom_text(aes(label = MAOC.let, y = MAOM_C_gm2_Mean+MAOM_C_gm2_SE), 
            position = position_dodge(0.9),
            vjust = -1, size = 10, size.unit = "pt") +
  theme_dc() +
  facet_wrap(~depth, nrow = 2, scales = "free_y") +
  coord_cartesian(clip = "off")
MAOCplot

fracCplot <- plot_grid(
  LFCplot + theme(legend.position = "none", 
                  text = element_text(size=14),
                  axis.title.y = element_blank()), 
  hPOCplot + theme(legend.position = "none", 
                   text = element_text(size=14),
                   axis.title.y = element_blank()), 
  MAOCplot + theme(legend.position = "none", 
                   text = element_text(size=14),
                   axis.title.y = element_blank()),
  nrow = 1,
  ncol = 3,
  align="hv"
)
fracCplot

#ggsave("DC_FracC.svg", plot = fracCplot, dpi = 300, height = 5, width = 6)

LFNplot <- ggplot(dc.frac.sum, 
                  aes(x=watershed, y=LF_N_gm2_Mean, fill=watershed))+
  geom_bar(position=position_dodge(), colour="black", stat="identity") + 
  geom_errorbar(aes(ymin=LF_N_gm2_Mean-LF_N_gm2_SE,
                    ymax=LF_N_gm2_Mean+LF_N_gm2_SE), 
                width=0.2, position=position_dodge(0.9)) +
  scale_fill_manual(values = dcPalette3) + 
  labs(y=bquote('Nitrogen per Fraction'~(g/m^2))) +
  geom_text(aes(label = LFN.let, y = LF_N_gm2_Mean+LF_N_gm2_SE), 
            position = position_dodge(0.9),
            vjust = -1, size = 10, size.unit = "pt") +
  theme_dc() +
  facet_wrap(~depth, nrow = 2, scales = "free_y") +
  coord_cartesian(clip = "off")
LFNplot

hPONplot <- ggplot(dc.frac.sum, 
                   aes(x=watershed, y=hPOM_N_gm2_Mean, fill=watershed))+
  geom_bar(position=position_dodge(), colour="black", stat="identity") + 
  geom_errorbar(aes(ymin=hPOM_N_gm2_Mean-hPOM_N_gm2_SE,
                    ymax=hPOM_N_gm2_Mean+hPOM_N_gm2_SE), 
                width=0.2, position=position_dodge(0.9)) +
  scale_fill_manual(values = dcPalette3) + 
  labs(y=bquote('Nitrogen per Fraction'~(g/m^2))) +
  geom_text(aes(label = hPON.let, y = hPOM_N_gm2_Mean+hPOM_N_gm2_SE), 
            position = position_dodge(0.9),
            vjust = -1, size = 10, size.unit = "pt") +
  theme_dc() +
  facet_wrap(~depth, nrow = 2, scales = "free_y") +
  coord_cartesian(clip = "off")
hPONplot

MAONplot <- ggplot(dc.frac.sum, 
                   aes(x=watershed, y=MAOM_N_gm2_Mean, fill=watershed))+
  geom_bar(position=position_dodge(), colour="black", stat="identity") + 
  geom_errorbar(aes(ymin=MAOM_N_gm2_Mean-MAOM_N_gm2_SE,
                    ymax=MAOM_N_gm2_Mean+MAOM_N_gm2_SE), 
                width=0.2, position=position_dodge(0.9)) +
  scale_fill_manual(values = dcPalette3) + 
  labs(y=bquote('Nitrogen per Fraction'~(g/m^2))) +
  geom_text(aes(label = MAON.let, y = MAOM_N_gm2_Mean+MAOM_N_gm2_SE), 
            position = position_dodge(0.9),
            vjust = -1, size = 10, size.unit = "pt") +
  theme_dc() +
  facet_wrap(~depth, nrow = 2, scales = "free_y") +
  coord_cartesian(clip = "off")
MAONplot

fracNplot <- plot_grid(
  LFNplot + theme(legend.position = "none", 
                  text = element_text(size=14),
                  axis.title.y = element_blank()), 
  hPONplot + theme(legend.position = "none", 
                   text = element_text(size=14),
                   axis.title.y = element_blank()), 
  MAONplot + theme(legend.position = "none", 
                   text = element_text(size=14),
                   axis.title.y = element_blank()),
  nrow = 1,
  ncol = 3,
  align="hv"
)
fracNplot

#ggsave("DC_FracN.svg", plot = fracNplot, dpi = 300, height = 5, width = 6)

rm(rootplot, soilCplot, soilNplot, soilCNplot, LFCplot, hPOCplot, MAOCplot,
   fracCplot) #clean plots

MAOC.hPOCplot <- ggplot(dc.frac.sum, 
                  aes(x=watershed, y=MAOC.hPOC_Mean, fill=watershed))+
  geom_bar(position=position_dodge(), colour="black", stat="identity") + 
  geom_errorbar(aes(ymin=MAOC.hPOC_Mean-MAOC.hPOC_SE,
                    ymax=MAOC.hPOC_Mean+MAOC.hPOC_SE), 
                width=0.2, position=position_dodge(0.9)) +
  scale_fill_manual(values = dcPalette3) + 
  labs(y=bquote('MAOC-to-heavy POC')) +
  geom_text(aes(label = MAOC.hPOC.let, y = MAOC.hPOC_Mean+MAOC.hPOC_SE), 
            position = position_dodge(0.9),
            vjust = -1, size = 10, size.unit = "pt") +
  theme_dc() +
  facet_wrap(~depth, nrow = 2, scales = "free_y") +
  coord_cartesian(clip = "off")
MAOC.hPOCplot

MAOC.tot_POCplot <- ggplot(dc.frac.sum, 
                           aes(x=watershed, y=MAOC.tot_POC_Mean, fill=watershed))+
  geom_bar(position=position_dodge(), colour="black", stat="identity") + 
  geom_errorbar(aes(ymin=MAOC.tot_POC_Mean-MAOC.tot_POC_SE,
                    ymax=MAOC.tot_POC_Mean+MAOC.tot_POC_SE), 
                width=0.2, position=position_dodge(0.9)) +
  scale_fill_manual(values = dcPalette3) + 
  labs(y=bquote('MAOC-to-total POC')) +
  geom_text(aes(label = MAOC.tPOC.let, y = MAOC.tot_POC_Mean+MAOC.tot_POC_SE), 
            position = position_dodge(0.9),
            vjust = -1, size = 10, size.unit = "pt") +
  theme_dc() +
  facet_wrap(~depth, nrow = 2, scales = "free_y") +
  coord_cartesian(clip = "off")
MAOC.tot_POCplot

MAOC.POCplot <- plot_grid(
  MAOC.hPOCplot + theme(legend.position = "none", 
                        text = element_text(size=14),
                        axis.title.y = element_blank()), 
  MAOC.tot_POCplot + theme(legend.position = "none", 
                           text = element_text(size=14),
                           axis.title.y = element_blank()),
  nrow = 1,
  ncol = 2,
  align="hv"
)
MAOC.POCplot

ggsave("DC_MAOC.POC.svg", plot = MAOC.POCplot, dpi = 300, height = 5, width = 3.5)


MAON.tot_PONplot <- ggplot(dc.frac.sum, 
                           aes(x=watershed, y=MAON.tot_PON_Mean, fill=watershed))+
  geom_bar(position=position_dodge(), colour="black", stat="identity") + 
  geom_errorbar(aes(ymin=MAON.tot_PON_Mean-MAON.tot_PON_SE,
                    ymax=MAON.tot_PON_Mean+MAON.tot_PON_SE), 
                width=0.2, position=position_dodge(0.9)) +
  scale_fill_manual(values = dcPalette3) + 
  labs(y=bquote('MAON-to-total PON')) +
  theme_dc() +
  facet_wrap(~depth, nrow = 2, scales = "free_y") +
  coord_cartesian(clip = "off")
MAON.tot_PONplot


#### Linear Mixed Effects Modeling ####

#Fracs only available for some samples -> running without frac data first

dc.NA <- dc[,1:11]
dc.NA <- na.omit(dc.NA) #4 NA rows (samples) removed

#Identify best model that explains total soil C

soilClme <- lmer(soil_C_gm2 #response
                 ~ fine_root_gm2 + #fixed
                   soil_CN +
                   soil_N_gm2 +
                   (1|watershed/depth), #random
                 data = dc.NA,
                 na.action = "na.fail",
                 REML = FALSE)

#Examine all possible models for soil C

soilClme_output <- dredge(soilClme) #singularity warning
#may be collinearity (which we will check) or
#when a random effect variance is estimated very near zero and the data is not 
#sufficiently informative to drag the estimate away from the zero starting value

soilClme_best <- get.models(soilClme_output, 1)[[1]]
soilClme_best

#Model summary for soil C

summary(soilClme_best)
check_collinearity(soilClme_best) #VIF < 2 Low Correlation

r.squaredGLMM(soilClme_best) #marginal = 0.96, conditional = 0.97

#Standardizing coefficients:

#Method 1:

stdCoef.merMod <- function(object) {
  sdy <- sd(getME(object,"y"))
  sdx <- apply(getME(object,"X"), 2, sd)
  sc <- fixef(object)*sdx/sdy
  se.fixef <- coef(summary(object))[,"Std. Error"]
  se <- se.fixef*sdx/sdy
  return(data.frame(stdcoef=sc, stdse=se))
} #function

stdCoef.merMod(soilClme_best)

#Method 2:

std.coef(soilClme_best, partial = TRUE) #MuMIn package

#Both methods produce a similar ratio & relationships stay the same
#Going with MuMIn package values (Method 2)

fig.lmer <- as.data.frame(std.coef(soilClme_best, partial = TRUE))[-1,-3]
fig.lmer <- rownames_to_column(fig.lmer)
fig.lmer$depths <- rep(c("All"), times = 3)
fig.lmer$response <- rep(c("Soil C"), times = 3)

r.squaredGLMM(soilClme_best)[1] #marginal R^2 only
summary(soilClme_best)$AICtab[1] #AIC only
#extractAIC(OH.soilClme_best) #alternative AIC

#Raw Outputs:

raw.soilC <- as.data.frame(summary(soilClme_best)[[10]])[-1,-4]
raw.soilC <- rownames_to_column(raw.soilC)
vifs <- as.data.frame(check_collinearity(soilClme_best))[,1:2] #VIFs
colnames(vifs) <- c("rowname", "VIF")
raw.soilC <- left_join(raw.soilC, vifs, by = "rowname")

#Figure 2: Std coefficients & St Errors for soil C, MAOC, and POC

str(fig.lmer)

fig.lmer <- fig.lmer %>% 
  mutate(across(c(rowname, depths, response), as.factor))

colnames(fig.lmer) <- c("Effect", "Std_Coef", "SE", "Depth", "Response")
levels(fig.lmer$Effect) <- c("Fine Root Biomass", "Soil C:N", "Soil N")

fig.lmer.g <- ggplot(fig.lmer, aes(x=Effect, y=Std_Coef))+
  geom_bar(position=position_dodge(),stat="identity") +
  geom_errorbar(aes(ymin=Std_Coef-SE,ymax=Std_Coef+SE), 
                width=0.2,position=position_dodge(0.9)) +
  theme_dc() + 
  scale_fill_discrete(type = dcPalette3[2]) +
  labs(y="Standardized Coefficient", x = element_blank())
fig.lmer.g

figtbl <- as.data.frame(matrix(,2,2))
colnames(figtbl) <- c("Variable", "Total Soil C")
figtbl$Variable <- c("AIC", "Marginal R2")

figtbl[1,2] <- summary(soilClme_best)$AICtab[1] #fill AIC
figtbl[2,2] <- r.squaredGLMM(soilClme_best)[1] #fill marginal R2

figtbl.plot <- tableGrob(figtbl, theme = ttheme_minimal(), rows = NULL)
grid.newpage()
grid.draw(figtbl.plot)

LMEplot <- plot_grid(
  fig.lmer.g, 
  figtbl.plot,
  nrow = 2,
  ncol = 1,
  rel_heights = c(4,1)
)

ggsave("DC_LME.png", plot = LMEplot, dpi = 300, height = 5, width = 3.55)


#TODO: 
# NMDS plots? Potentially


#### April 2026 ####
#For Eddie: Soil C per m2 (area of sampling) (soilC_gm2) (*across depths*)
#For me: Run Dunnett (instead of Tukey)

soilC.aov <- aov(soil_C_gm2 ~ watershed, dc)
summary(soilC.aov) #F=1.389, p=0.252

soilC.dun <- glht(soilC.aov, linfct=mcp(watershed="Dunnett"))
summary(soilC.dun) #DL: 0.953, ST: 0.211

#Total C: ns, dunnett for tot c: ns

soilC.OH.dun <- glht(soilC.OH.aov, linfct=mcp(watershed="Dunnett"))
summary(soilC.OH.dun) #DL: 0.00884, ST: 0.10103

soilC.20.dun <- glht(soilC.20.aov, linfct=mcp(watershed="Dunnett"))
summary(soilC.20.dun) #DL: 0.7020, ST: 0.0588* the only change

soilC.40.dun <- glht(soilC.40.aov, linfct=mcp(watershed="Dunnett"))
summary(soilC.40.dun) #DL: 0.875, ST: 0.301

#Dunnett's for other vars

roots.OH.dun <- glht(roots.OH.aov, linfct=mcp(watershed="Dunnett"))
summary(roots.OH.dun) #DL: 0.124, ST: 0.955

roots.20.dun <- glht(roots.20.aov, linfct=mcp(watershed="Dunnett"))
summary(roots.20.dun) #DL: 0.996, ST: 0.994

roots.40.dun <- glht(roots.40.aov, linfct=mcp(watershed="Dunnett"))
summary(roots.40.dun) #DL: 0.0378, ST: 0.0360 * the only change



soilN.OH.dun <- glht(soilN.OH.aov, linfct=mcp(watershed="Dunnett"))
summary(soilN.OH.dun) #DL: 0.00698, ST: 0.08490

soilN.20.dun <- glht(soilN.20.aov, linfct=mcp(watershed="Dunnett"))
summary(soilN.20.dun) #DL: 0.1231, ST: 0.0205

soilN.40.dun <- glht(soilN.40.aov, linfct=mcp(watershed="Dunnett"))
summary(soilN.40.dun) #DL: 0.999, ST: 0.790



soilCN.OH.dun <- glht(soilCN.OH.aov, linfct=mcp(watershed="Dunnett"))
summary(soilCN.OH.dun) #DL: 0.142, ST: 0.968

soilCN.20.dun <- glht(soilCN.20.aov, linfct=mcp(watershed="Dunnett"))
summary(soilCN.20.dun) #DL: 0.281, ST: 0.905

soilCN.40.dun <- glht(soilCN.40.aov, linfct=mcp(watershed="Dunnett"))
summary(soilCN.40.dun) #DL: 0.5514, ST: 0.0687



MAOC.tPOC.20.dun <- glht(MAOC.tPOC.20.aov, linfct=mcp(watershed="Dunnett"))
summary(MAOC.tPOC.20.dun) #DL: 0.673, ST: 0.650

MAOC.tPOC.40.dun <- glht(MAOC.tPOC.40.aov, linfct=mcp(watershed="Dunnett"))
summary(MAOC.tPOC.40.dun) #DL: 0.8620, ST: 0.0859



MAON.tPON.20.dun <- glht(MAON.tPON.20.aov, linfct=mcp(watershed="Dunnett"))
summary(MAON.tPON.20.dun) #DL: 0.387, ST: 0.570

MAON.tPON.40.dun <- glht(MAON.tPON.40.aov, linfct=mcp(watershed="Dunnett"))
summary(MAON.tPON.40.dun) #DL: 0.663, ST: 0.589




LF.C.20.dun <- glht(LF.C.20.aov, linfct=mcp(watershed="Dunnett"))
summary(LF.C.20.dun) #n.s.

LF.C.40.dun <- glht(LF.C.40.aov, linfct=mcp(watershed="Dunnett"))
summary(LF.C.40.dun) ##DL: 0.00735, ST: n.s.
