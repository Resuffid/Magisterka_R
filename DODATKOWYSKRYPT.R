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
library(car)

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

cor(df3$Class, df4$KWS_G, method ='pearson')
cor.test(df3$Class, df4$KWS_G, method ='pearson')
cor.test(df3$Class, df4$KWS_EES, method ='pearson')
cor.test(df3$Class, df4$KWS_A, method ='pearson')
cor.test(df3$Class, df4$KWS_QoL, method ='pearson')

ggplot(df4, mapping = aes(x=Class, y=KWS_G))+
  geom_jitter()+
  geom_smooth(method = "lm")
fit.classkwsg <- lm(Class ~ KWS_G, df4)
summary(fit.classkwsg)
plot(fit.classkwsg)

ggplot(df4, mapping = aes(x=Class, y=KWS_EES))+
  geom_jitter()+
  geom_smooth(method = "lm")
fit.classkwsees <- lm(Class ~ KWS_EES, df4)
summary(fit.classkwsees)
plot(fit.classkwsees)
cohen.d(fit.classkwsees)


ggplot(df4, mapping = aes(x=Class, y=KWS_A))+
  geom_jitter()+
  geom_smooth(method = "lm")
fit.classkwsa <- lm(Class ~ KWS_A, df4)
summary(fit.classkwsa)
plot(fit.classkwsa)

ggplot(df4, mapping = aes(x=Class, y=KWS_QoL))+
  geom_jitter()+
  geom_smooth(method = "lm")
fit.classkwsqol <- lm(Class ~ KWS_QoL, df4)
summary(fit.classkwsqol)
plot(fit.classkwsqol)

ggplot(df4, mapping = aes(x=Class, y=KWS_QoL))+
  geom_jitter()+
  geom_smooth(method = "lm")

fit.classkwsqol <- lm(Class ~ KWS_G+KWS_EES+KWS_A+KWS_QoL, df4)
summary(fit.classkwsqol)
plot(fit.classkwsqol)

contrasts(df4$klasa)
df4$klasa<-as.factor(df4$Class)
modelklasakmsg <- lm(KWS_G ~ klasa, data = df4)
summary(modelklasakmsg)
anova(modelklasakmsg)
k_wielomianowe1 <- contrast(emmeans(modelklasakmsg, "klasa"), method = "poly")
summary(k_wielomianowe1)
ggplot(df4, aes(x = Class, y = KWS_G)) +
  geom_boxplot(aes(group = klasa, fill = klasa), alpha = 0.4, show.legend = FALSE) +
  
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), 
              color = "darkred", size = 1.2, se = TRUE, alpha = 0.15) +
  scale_x_continuous(breaks = 1:length(levels(df4$klasa)), labels = levels(df4$klasa)) +
  theme_minimal() +
  labs(
    title = "Box-and-whiskers plot with trendline representing relationship betweeen general score of KWS and year of secondary school",
    x = "Year of secondary school",
    y = "KWS G"
  )

modelklasakwsees <- lm(KWS_EES ~ klasa, data = df4)
summary(modelklasakwsees)
anova(modelklasakwsees)
k_wielomianowe2 <- contrast(emmeans(modelklasakwsees, "klasa"), method = "poly")
summary(k_wielomianowe2)
ggplot(df4, aes(x = Class, y = KWS_EES)) +
  geom_boxplot(aes(group = klasa, fill = klasa), alpha = 0.4, show.legend = FALSE) +
  
  geom_smooth(method = "lm", formula = y ~ poly(x, 4), 
              color = "darkred", size = 1.2, se = TRUE, alpha = 0.15) +
  scale_x_continuous(breaks = 1:length(levels(df4$klasa)), labels = levels(df4$klasa)) +
  theme_minimal() +
  labs(
    title = "Box-and-whiskers plot with trendline representing relationship betweeen EES subscale score of KWS and year of secondary school",
    x = "Year of secondary school",
    y = "KWS EES"
  )


modelklasakwsa <- lm(KWS_A ~ klasa, data = df4)
summary(modelklasakwsa)
anova(modelklasakwsa)
k_wielomianowe3 <- contrast(emmeans(modelklasakwsa, "klasa"), method = "poly")
summary(k_wielomianowe3)



modelklasakwsqol <- lm(KWS_QoL ~ klasa, data = df4)
summary(modelklasakwsqol)
anova(modelklasakwsqol)
k_wielomianowe4 <- contrast(emmeans(modelklasakwsqol, "klasa"), method = "poly")
summary(k_wielomianowe4)
ggplot(df4, aes(x = Class, y = KWS_QoL)) +
  geom_boxplot(aes(group = klasa, fill = klasa), alpha = 0.4, show.legend = FALSE) +
  
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), 
              color = "darkred", size = 1.2, se = TRUE, alpha = 0.15) +
  scale_x_continuous(breaks = 1:length(levels(df4$klasa)), labels = levels(df4$klasa)) +
  theme_minimal() +
  labs(
    title = "Box-and-whiskers plot with trendline representing relationship betweeen QoL subscale score of KWS and year of secondary school",
    x = "Year of secondary school",
    y = "KWS G"
  )



#ANOVA
KWS_GAOV <- aov(KWS_G ~ Living, data = df4)
KWS_GAOV
summary(KWS_GAOV)
TukeyHSD(KWS_GAOV)
KWS1 <- ggplot(data=df4, mapping = aes(x=Living, y=KWS_G)) +
  geom_boxplot()+
  labs(x="Place of Living", y="KWS G")
eta_squared(KWS_GAOV)

KWS_EESAOV <- aov(KWS_EES ~ Living, data = df4)
KWS_EESAOV
summary(KWS_EESAOV)
TukeyHSD(KWS_EESAOV)
KWS2 <- ggplot(data=df4, mapping = aes(x=Living, y=KWS_EES)) +
  geom_boxplot()+
  labs(x="Place of Living", y="KWS EES")
eta_squared(KWS_EESAOV)

KWS_AAOV <- anova(lm(KWS_A ~ Living, data = df4))
KWS_AAOV
eta_squared(KWS_AAOV)

KWS_QoLAOV <- aov(KWS_QoL ~ Living, data = df4)
KWS_QoLAOV
summary(KWS_QoLAOV)
TukeyHSD(KWS_QoLAOV)
KWS3 <- ggplot(data=df4, mapping = aes(x=Living, y=KWS_QoL)) +
  geom_boxplot()+
  labs(x="Place of Living", y="KWS QoL")
eta_squared(KWS_QoLAOV)
combined_KWS <- grid.arrange(KWS1, KWS2, KWS3, ncol=3)

FRAS_FCPSAOV <- aov(FRAS_FCPS ~ Living, data = df4)
FRAS_FCPSAOV
summary(FRAS_FCPSAOV)
TukeyHSD(FRAS_FCPSAOV)
ggplot(data=df4, mapping = aes(x=Living, y=FRAS_FCPS)) +
  geom_boxplot()+
  labs(x="Place of Living", y="FRAS FCPS")
eta_squared(FRAS_FCPSAOV)

FRAS_FCPSAOV <- aov(FRAS_FCPS ~ Living, data = df4)
FRAS_FCPSAOV
summary(FRAS_FCPSAOV)
TukeyHSD(FRAS_FCPSAOV)
ggplot(data=df4, mapping = aes(x=Living, y=FRAS_FCPS)) +
  geom_boxplot()+
  labs(x="Place of Living", y="FRAS FCPS")
eta_squared(FRAS_FCPSAOV)

FRAS_USERAOV <- aov(FRAS_USER ~ Living, data = df4)
FRAS_USERAOV
summary(FRAS_USERAOV)
TukeyHSD(FRAS_USERAOV)
FRAS1 <- ggplot(data=df4, mapping = aes(x=Living, y=FRAS_USER)) +
  geom_boxplot()+
  labs(x="Place of Living", y="FRAS USER")
eta_squared(FRAS_USERAOV)

FRAS_MPOAOV <- aov(FRAS_MPO ~ Living, data = df4)
FRAS_MPOAOV
summary(FRAS_MPOAOV)
TukeyHSD(FRAS_MPOAOV)
ggplot(data=df4, mapping = aes(x=Living, y=FRAS_MPO)) +
  geom_boxplot()+
  labs(x="Place of Living", y="FRAS MPO")
eta_squared(FRAS_MPOAOV)

FRAS_FCAOV <- aov(FRAS_FC ~ Living, data = df4)
FRAS_FCAOV
summary(FRAS_FCAOV)
TukeyHSD(FRAS_FCAOV)
ggplot(data=df4, mapping = aes(x=Living, y=FRAS_FC)) +
  geom_boxplot()+
  labs(x="Place of Living", y="FRAS FC")
eta_squared(FRAS_FCAOV)

FRAS_FSAOV <- aov(FRAS_FS ~ Living, data = df4)
FRAS_FSAOV
summary(FRAS_FSAOV)
TukeyHSD(FRAS_FSAOV)
ggplot(data=df4, mapping = aes(x=Living, y=FRAS_FS)) +
  geom_boxplot()+
  labs(x="Place of Living", y="FRAS FS")
eta_squared(FRAS_FSAOV)

FRAS_AMMAAOV <- aov(FRAS_AMMA ~ Living, data = df4)
FRAS_AMMAAOV
summary(FRAS_AMMAAOV)
TukeyHSD(FRAS_AMMAAOV)
ggplot(data=df4, mapping = aes(x=Living, y=FRAS_AMMA)) +
  geom_boxplot()+
  labs(x="Place of Living", y="FRAS AMMA")
eta_squared(FRAS_AMMAAOV)

FRAS_ALLAOV <- aov(FRAS_ALL ~ Living, data = df4)
FRAS_ALLAOV
summary(FRAS_ALLAOV)
TukeyHSD(FRAS_ALLAOV)
FRAS2 <- ggplot(data=df4, mapping = aes(x=Living, y=FRAS_ALL)) +
  geom_boxplot()+
  labs(x="Place of Living", y="FRAS ALL")
eta_squared(FRAS_ALLAOV)

grid.arrange(FRAS1, FRAS2, ncol=2)

JiMSLAOV <- aov(JiMSz_L ~ Living, data = df4)
JiMSLAOV
summary(JiMSLAOV)
TukeyHSD(JiMSLAOV)
ggplot(data=df4, mapping = aes(x=Living, y=JiMSz_L)) +
  geom_boxplot()+
  labs(x="Place of Living", y="JiMSz A")
eta_squared((JiMSLAOV))

JiMSMAOV <- aov(JiMSz_M ~ Living, data = df4)
JiMSMAOV
summary(JiMSMAOV)
TukeyHSD(JiMSMAOV)
ggplot(data=df4, mapping = aes(x=Living, y=JiMSz_M)) +
  geom_boxplot()+
  labs(x="Place of Living", y="JiMSz M")

JiMSKAOV <- aov(JiMSz_K ~ Living, data = df4)
JiMSKAOV
summary(JiMSKAOV)
TukeyHSD(JiMSKAOV)
ggplot(data=df4, mapping = aes(x=Living, y=JiMSz_K)) +
  geom_boxplot()+
  labs(x="Place of Living", y="JiMSz L")
eta_squared(JiMSKAOV)

#H3
cor_mat(KWS_G, KWS_EES, KWS_A, KWS_QoL, FRAS_ALL, FRAS_FCPS, FRAS_USER, FRAS_MPO, FRAS_FC, FRAS_FS, FRAS_AMMA, data=df4)
cor_pmat(KWS_G, KWS_EES, KWS_A, KWS_QoL, FRAS_ALL, FRAS_FCPS, FRAS_USER, FRAS_MPO, FRAS_FC, FRAS_FS, FRAS_AMMA, data=df4)

model_KWSG <- lm(KWS_G ~ FRAS_FCPS+FRAS_USER+FRAS_MPO+FRAS_FC+FRAS_FS+FRAS_AMMA, data = df4)
summary(model_KWSG)
plot(model_KWSG)
sigma(model_KWSG)
vif(model_KWSG)
avPlots(model_KWSG)

model_KWSEES <- lm(KWS_EES ~ FRAS_FCPS+FRAS_USER+FRAS_MPO+FRAS_FC+FRAS_FS+FRAS_AMMA, data = df4)
summary(model_KWSEES)
plot(model_KWSEES)
sigma(model_KWSEES)
vif(model_KWSEES)

model_KWSA <- lm(KWS_A ~ FRAS_FCPS+FRAS_USER+FRAS_MPO+FRAS_FC+FRAS_FS+FRAS_AMMA, data = df4)
summary(model_KWSA)
plot(model_KWSA)
sigma(model_KWSA)
vif(model_KWSA)

model_KWSQoL <- lm(KWS_QoL ~ FRAS_FCPS+FRAS_USER+FRAS_MPO+FRAS_FC+FRAS_FS+FRAS_AMMA, data = df4)
summary(model_KWSQoL)
plot(model_KWSQoL)
sigma(model_KWSQoL)
vif(model_KWSQoL)

summary(lm(KWS_G~FRAS_ALL, data=df4))
