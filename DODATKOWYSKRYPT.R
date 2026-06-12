library(tidyverse)
library(psych)
library(effsize)
library(EnvStats)
library(lme4)
library(lmerTest)
library(emmeans)
library(datasets)
library(broom)
library(ggpubr)
library(rstatix)
library(gridExtra)
library(multcompView)
library(lavaan)
library(foreign) 
library(GPArotation)

df3 <- read.csv("df3.csv")

model <- 'KWS_SS =~ KWS_1+KWS_4+KWS_8
          KWS_WOD =~ KWS_2+KWS_5+KWS_13+KWS_9+KWS_16
          KWS_W =~ KWS_6R+KWS_14R+KWS_11+KWS_18+KWS_19
          KWS_JZ =~ KWS_3R+KWS_7R+KWS_12R+KWS_15R
          KWS_CYN =~ KWS_10+KWS_17
          KWS_OG =~ KWS_SS+KWS_WOD+KWS_W+KWS_JZ+KWS_CYN'
fit <- cfa(model, data=df3, estimator="WLSMV")
summary(fit, standardized = TRUE, ci = TRUE, fit.measures=TRUE)

KWS <- df3 %>%
  select(KWS_1,KWS_2,KWS_3R,KWS_4,KWS_5,KWS_6R,KWS_7R,KWS_8,KWS_9,KWS_10,KWS_11,KWS_12R,KWS_13,KWS_14R,KWS_15R,KWS_16,KWS_17,KWS_18,KWS_19)
omega(KWS, nfactors = 5)
KWS_SS <- df3%>%
  select(KWS_1,KWS_4,KWS_8)
KWS_WOD <- df3%>% select(KWS_2,KWS_5,KWS_13,KWS_9,KWS_16)
KWS_W <-df3%>% select(KWS_6R,KWS_14R,KWS_11,KWS_18,KWS_19)
KWS_JZ <- df3%>% select(KWS_3R,KWS_7R,KWS_12R,KWS_15R)
KWS_CYN = df3%>% select(KWS_10,KWS_17)
alpha(KWS)
alpha(KWS_SS)
alpha(KWS_WOD)
alpha(KWS_W)
alpha(KWS_JZ)
alpha(KWS_CYN)

model <- 'KWS_SS =~ KWS_4+KWS_8
          KWS_WOD =~ KWS_2+KWS_4+KWS_5+KWS_8+KWS_13+KWS_9+KWS_16
          KWS_W =~ KWS_11+KWS_18+KWS_19
          KWS_JZ =~ KWS_3R+KWS_7R+KWS_12R+KWS_15R
          KWS_OG =~ KWS_SS+KWS_WOD+KWS_W+KWS_JZ'
fit <- cfa(model, data=df3, estimator="WLSMV")
summary(fit, standardized = TRUE, ci = TRUE, fit.measures=TRUE)

KWS <- df3 %>%
  select(KWS_2,KWS_3R,KWS_4,KWS_5,KWS_7R,KWS_8,KWS_9,KWS_11,KWS_12R,KWS_13,KWS_15R,KWS_16,KWS_18,KWS_19)
omega(KWS, nfactors = 4)
KWS_SS <- df3%>%
  select(KWS_4,KWS_8)
KWS_WOD <- df3%>% select(KWS_2,KWS_5,KWS_13,KWS_16)
KWS_W <-df3%>% select(KWS_11,KWS_18,KWS_19)
KWS_JZ <- df3%>% select(KWS_3R,KWS_7R,KWS_12R,KWS_15R)

alpha(KWS)
alpha(KWS_SS)
alpha(KWS_WOD)
alpha(KWS_W)
alpha(KWS_JZ)

efa(KWS, nfactors = 4)
omega(KWS, nfactors = 4)
efa(KWS)
scree(KWS, FA)
scree(KWS,factors=TRUE,pc=TRUE,main="Scree plot",hline=NULL,add=FALSE,sqrt=FALSE)

FitEFA_6 <- fa(r = KWS, 
               nfactors = 3, 
               rotate = "varimax")
print(FitEFA_6)
library(psy)
scree.plot(KWS)
KMO(KWS)
cortest.bartlett(KWS)
ev <- eigen(cor(KWS)) # get eigenvalues
ev$values
scree(KWS, pc=FALSE)
fa.parallel(KWS, fa="fa")
Nfacs <- 3
fit <- factanal(KWS, Nfacs, rotation="varimax")
print(fit, digits=2, cutoff=0.3, sort=TRUE)
loads <- fit$loadings
fa.diagram(loads)

model1 <- '
          KWS_WOD =~ KWS_2+KWS_4+KWS_5+KWS_8+KWS_9+KWS_13+KWS_16
          KWS_W =~ KWS_11+KWS_18+KWS_19
          KWS_JZ =~ KWS_3R+KWS_7R+KWS_12R+KWS_15R
          KWS_OG =~ KWS_WOD+KWS_W+KWS_JZ'
fit <- cfa(model1, data=df3, estimator="WLSMV")
summary(fit, standardized = TRUE, ci = TRUE, fit.measures=TRUE)

KWSost <- df3%>% select(KWS_2,KWS_4,KWS_5,KWS_8,KWS_9,KWS_13,KWS_16,
                        KWS_11,KWS_18,KWS_19,
                        KWS_3R,KWS_7R,KWS_12R,KWS_15R)
KWS_WOD <- df3%>% select(KWS_2,KWS_4,KWS_5,KWS_8,KWS_13,KWS_16)
KWS_W <- df3%>% select(KWS_11,KWS_18,KWS_19)
KWS_JZ <- df3%>% select(KWS_3R,KWS_7R,KWS_12R,KWS_15R)                        
omega(KWSost, nfactors = 3)
alpha(KWSost)
alpha(KWS_WOD)
alpha(KWS_W)
alpha(KWS_JZ)

df4 <- df3 %>%
  mutate(KWS_EES = KWS_2+KWS_4+KWS_5+KWS_8+KWS_13+KWS_16,
         KWS_A = KWS_11+KWS_18+KWS_19,
         KWS_QoL = KWS_3R+KWS_7R+KWS_12R+KWS_15R,
         KWS_G = KWS_WOD+KWS_W+KWS_JZ)
write.csv(df4, "df4.csv")
