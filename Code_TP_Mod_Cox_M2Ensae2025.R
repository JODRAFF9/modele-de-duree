
#Modele de Cox ###################

### Packages

library(survival) 
library(CoxR2)
library(timereg)
library(survminer) 
library(stringi)
library(muhaz)    
library(knitr)
library(gtsummary)
library(GGally)
library(riskRegression)
library(conflicted)
conflict_prefer("prop", "timereg")
conflict_prefer("prop","questionr") 
library(dplyr)
#install.packages("forestmodel")  # une seule fois
library(forestmodel)

#Donnees de cancer

data("lung")
attach(lung)

Data <- lung
View(Data)

head(lung)
?lung
str(lung)

#inst: Code Institution 
#time: temps de survie en jours
#status: indicateur de censure: 0=censure, 1=deces (evt)
#age: âge en annees
#sex: Male=1 Femele=2
#ph.ecog: ECOG performance score (0=good 5=dead)
#ph.karno: Karnofsky performance score (bad=0-good=100) rated by physician
#pat.karno: Karnofsky performance score as rated by patient
#meal.cal: Calories consumed at meals
#wt.loss: Weight loss in last six months

dim(lung)
# 228 observations et 10 variables 

table(status)

status=ifelse(status==2,1,0)

table(status)

attach(lung)

args(coxph)  # donnent les arguments de la commande
# Cox proportional hazards regression model

## 1) modele 1 univarie: effet variable sex 

#fonctions de survie

hist(time)
time

survie = Surv(time,status)
survie
?Surv


class(sex)
genre=as.factor(sex)
class(genre)
genre=ifelse(genre=="1","Hom","Fem")
table(genre)

mod_bres = coxph(survie~genre,method="breslow")
mod_bres
summary(mod_bres)

mod_exa = coxph(survie~genre, method="exact")
mod_exa
summary(mod_exa)

mods = coxph(survie~genre, method="efron")
mods
summary(mods)

coxr2(mods)$rsq # % de variation de variance expliquee par la variable
coxr2(mod_bres)$rsq
coxr2(mod_exa)$rsq 

## 1) modele 1 univarie: effet variable age --------

moda = coxph(survie~age, method="efron")
moda 
summary(moda)

summary(age)

ageb = cut(age,breaks=c(38,56,69,82))

table(ageb)
class(ageb)
levels(ageb)

levels(ageb) = c("C1","C2","C3")

moda1 = coxph(survie~ageb, method="efron")
summary(moda1)

agebb = relevel(ageb, ref = "C2")
table(agebb)

moda2 = coxph(survie~agebb, method="efron")
summary(moda2)

agebbb = relevel(ageb, ref = "C3")
table(agebbb)

moda3 = coxph(survie~agebbb, method="efron")
summary(moda3)

coxr2(moda1)$rsq # Pseudo R2
coxr2(moda2)$rsq
coxr2(moda3)$rsq


# modèle multivarie

table(ph.ecog)
sum(is.na(ph.ecog))
lung[is.na(ph.ecog),]

ph.ecog1=as.factor(ph.ecog)

table(ph.ecog1)

mod1=coxph(survie~age+genre+ph.ecog+ph.karno
           +meal.cal+wt.loss, method="efron")
summary(mod1)

mod2=coxph(survie~age+genre+ph.ecog1+ph.karno
           +meal.cal+wt.loss, method="efron")
summary(mod2)

# selection de variables

mod3=coxph(survie~age+genre+ph.ecog1+ph.karno
           +meal.cal+wt.loss, method="efron")
summary(mod3)

mod3.bis=coxph(survie~age+genre+ph.ecog1+ph.karno
           +meal.cal*wt.loss, method="efron")
summary(mod3.bis)

# On enlève meal.cal vu sa p-value

mod4=coxph(survie~age+genre+ph.ecog1+ph.karno
           +wt.loss, method="efron")
summary(mod4)

#  On enlève wt.loss vu sa p-value
mod5=coxph(survie~age+genre+ph.ecog1+ph.karno, method="efron")
summary(mod5)

# On enlève ph.karno vu sa p-value
mod6=coxph(survie~age+genre+ph.ecog1, method="efron")
summary(mod6)

mod6bis=coxph(survie~age*genre+ph.ecog1, method="efron")
summary(mod6bis)

# On enlève age vu sa p-value
mod7=coxph(survie~genre+ph.ecog1, method="efron")
summary(mod7)

# Comparaison des modeles

anova(mod7,mod6, test="Chisq")

AIC(mod6)
AIC(mod7)

BIC(mod6)
BIC(mod7)

coxr2(mod6)$rsq
coxr2(mod7)$rsq

# selection automatique backward

selectCox(survie~age+genre+ph.ecog1+ph.karno
  +meal.cal+wt.loss,data=lung)

# affichage

ggcoef_model(mod7, exponentiate = TRUE)
# faire dev.off() si le code précédent ne marche pas

forest_model(mod7) 

# Courbe de survie à partir du modèle

ggsurvplot(survfit(mod7,data=lung), palette= 'green',
           ggtheme = theme_minimal())

# estimation de la Courbe de survie à partir du modèle pour 2 individus

nd=with(lung,data.frame(genre=c("Fem","Hom"),ph.ecog1 = c("2","0")))
nd

m=survfit(mod7, newdata=nd,data=lung)
m
ggsurvplot(m, conf.int = TRUE, legend.labs=c("Sex=F", "Sex=H"),
           ggtheme = theme_minimal())

#B- Estimation non paramétrique -----
##1) Estimamtion de la fonction de hasard de base -------

frb = basehaz(mod7)
frb
plot(frb,type="l",col=5,lwd=2)


##2) Adequation du modele 

### a. hypothese de risque proportionnel: 
#condition d'indépendance des variables 
# explicatives en fonction du temps se basant 
#sur létude des résidus de Schoenfield
# H0: HRP respectée (indépendante du temps)  vs  
#H1: HRP violée (variation significative au cours du temps)


test.hrp = cox.zph(mod7)
test.hrp

### b. Diagnostic des résidus --------
#Vérification graphique des residus de Shoenfield ----
# Les écarts systématiques par rapport à la ligne 
# horizontale indiquent 
# des risques proportionnels

ggcoxzph(test.hrp) 

#### * Diagnostic des patients ayant une influence 
#forte sur les coefficients 

moddb = residuals(mod7, type="dfbeta")
moddb

par(mfrow=c(1,2))

for (j in 1:2){
  plot(moddb[,j], ylab=names(coef(mod7))[j])
  abline(h=0, lty=4)
}

par(mfrow=c(1,1))
plot(moddb[,2], ylab=names(coef(mod7))[2])
abline(v=9, lty=4)

ggcoxdiagnostics(mod7, type = "dfbeta", 
                 linear.predictions = FALSE, 
                 ggtheme = theme_bw())

survie1=survie[-48]
mod8=coxph(survie1~genre[-48]+ph.ecog1[-48], method="efron")
summary(mod8)

ggcoxdiagnostics(mod8, type = "dfbeta", # Modèle sans les individus 9 et 23
                 linear.predictions = FALSE, 
                 ggtheme = theme_bw())

AIC(mod7)
AIC(mod8)
coxr2(mod8)$rsq
coxr2(mod7)$rsq
BIC(mod7)
BIC(mod8)

#### * Residus de déviance (transformation normalisée 
#des résidu de martingale) ------
# Ces résidus doivent être distribués de façon à peu près 
#symétrique autour de zéro, avec un écart type de 1

ggcoxdiagnostics(mod7, type = "deviance",
                 linear.predictions = FALSE, 
                 ggtheme = theme_bw())

####Residus de martingale cumulés ---------
# erreur entre le modèle et les données cumulées au cours du temps

# Non linearite pour variable continue (log-linéarité)
# Variables dependant du temps

# ajustement du modèle nul
modnul=coxph(survie ~ 1, data = lung)
summary(modnul)

# stockage des résidus de Martinagle
rm0=residuals(modnul, type = "martingale")
rm7=residuals(mod7, type = "martingale")

# Tracé des résidus de Martingale en fonction de genre
plot(rm0~lung$age) 
abline(h=0,col=4)
g=genre[-1]
plot(rm7~g) 
