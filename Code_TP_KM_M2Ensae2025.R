# Packages

library(survival) 
library(CoxR2)
library(timereg)
library(survminer) 
library(stringi)
library(muhaz)    
library(knitr)
library(Hmisc)

### donnees de travail
# Survie dans un essai randomisé comparant 
# deux traitements contre le cancer de l'ovaire

data(ovarian)
str(ovarian)
?ovarian

# futime: temps de survenu de l'événement décès en semaine
# fustat: statut (0 = survie ou censure / 1 = décès)
# age: âge des patientes en années
# resid.ds: présence résiduel de maladies (1=non / 2=oui)
# rx: groupe de traitement (1=placebo / 2=traitement)
# ecog.ps: échelle de performance du patient (1=moyenne / 2 = faible)

class(ovarian$rx)
# Reconversion des donnees

ovarian$rx = as.factor(ovarian$rx)
table(ovarian$rx)
str(ovarian$rx)
levels(ovarian$rx)

prop.table(table(ovarian$rx))

ovarian$resid.ds = as.factor(ovarian$resid.ds)
table(ovarian$resid.ds)
str(ovarian$resid.ds)

ovarian$ecog.ps = as.factor(ovarian$ecog.ps)
table(ovarian$ecog.ps)
str(ovarian$ecog.ps)

str(ovarian)
attach(ovarian)

# statistiques descriptives

summary(ovarian)

kable(summary(ovarian[,c(1,3)]))

describe(age)

describe(futime)

ggplot(ovarian, aes( x=rx,y=futime, fill=rx, colour=rx)) +
  geom_jitter(width=0.25) +
  geom_boxplot(alpha=0.5)+ 
  xlab(label = "Traitement") + 
  ylab(label = "Duree") + 
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))+
  theme(legend.position="none")+
  theme_classic()+ 
  ggtitle("Comparaison des traitements") 

table(fustat)
boxplot(futime)
### Estimation de la fonction de survie

cbind(futime,fustat)
survie = Surv(futime, fustat)
survie

data.frame(fustat,survie)

# sans stratification

fit = survfit(survie~1, data=ovarian)  
fit
summary(fit)

summary(futime) # temps min et temps max

# Fonction de survie S(t)
par(mfrow=c(1,1))

plot(fit,col=4, xlab="Temps (semaine)", 
     ylab="S(t)",lwd=2)
abline(v=638,col=2,lwd=2)
abline(h=0.5,col="green",lwd=2)
ggsurvplot(fit,conf.int = T)

# survie cumulee

ggsurvplot(fit, fun = "event", risk.table = TRUE, 
           surv.scale = "percent", break.time.by = 12)

# avec stratification (variable rx !: placebo, 2: traitement)

fit2 = survfit(survie~rx,data=ovarian)  
fit2
summary(fit2)

plot(fit2,lty=c(2,1),lwd=c(2,2),col=c("black", "blue"), 
     conf.int=F, xlab="Temps (semaine)", ylab="S(t)")

legend("bottomright", c("Placebo", "Traitement"), 
       lty=c(2,1), lwd=c(1,2), 
       col=c("black","blue"), title="Groupes")

ggsurvplot(fit2)   # sans IC
ggsurvplot(fit2,conf.int = T) # avec IC

levels(resid.ds)
levels(rx)

fit2b = survfit(survie ~ resid.ds + rx , data=ovarian)
summary(fit2b)

ggsurvplot(fit2b,conf.int = F)
ggsurvplot(fit2b,conf.int = T)

# decomposer la variable "age" en classe

summary(age)

age2 = cut(age, breaks = 3)
str(age2)
table(age2)

age3=cut(age,breaks=c(38,50,62,75))
table(age3)

## modalites de age2
levels(age2)

levels(age2) = c("c1","c2","c3")
table(age2)

# pour choisir la categorie de reference

age2ref = relevel(age2,ref="c2")
levels(age2ref)

# estimation des fonctions de survie

fit3 = survfit(survie~age2,data=ovarian)  # avec stratification
fit3
summary(fit3)

plot(fit3,lty=c(2,2,2),lwd=c(2,2,2),col=c("black", "blue","red"), conf.int=F, 
     xlab=c("Temps (semaine)"), ylab=c("S(t)"))

legend("bottomright", c("(38.9,50.8]", "(50.8,62.6]","(62.6,74.5]"), lwd=c(1,1,1), 
       col=c("black","blue","red"),title="Groupes")

ggsurvplot(fit3)   # sans IC
ggsurvplot(fit3,conf.int = T)   # avec IC

### Estimateur de Nelson-Aalen du risque cumule

alen=summary(fit)
alen
risk=alen$n.event/alen$n.risk
risk
cbind(alen$time,round(risk,3))
p=cumsum(risk)
p
plot(alen$time,p, type="s", ylab="FRCNA",
     main= "Estimateur de Nelson-Alen du risque cumulé")

### Estimateur de Breslow du risque cumule H(t)=-log(S(t))

plot(fit$time,-log(fit$surv),type="s",col="red",
     xlab='time',
     main="Estimateur de Breslow du risque cumulé")

par(mfrow=c(2,1))
plot(alen$time,p, type="s", ylab="FRCNA",
     main= "Estimateur de Nelson-Alen du risque cumulé")

plot(fit$time,-log(fit$surv),type="s",col="blue",
     xlab='time',
     main="Estimateur de Breslow du risque cumulé")

# deux courbes dans le meme graphe

par(mfrow=c(1,1))

plot(alen$time,p, type="s", ylab="risque cumulé",
     main= "Estimateur de Nelson-Alen du risque cumulé")

lines(fit$time,-log(fit$surv),type="s",col="blue",
      xlab='time', xlim=c(50,600),
      main="Estimateur de Breslow du risque cumulé")

# risque cumule pour deux groupes

model1 = survfit(survie~rx,data=ovarian, subset=rx=="1")
model2 = survfit(survie~rx,data=ovarian, subset=rx=="2") 

plot(model1$time,-log(model1$surv),type="s",col="red",
     main="Estimateur de Breslow du risque cumulé des deux groupes",ylab="Risque cumulé") 
lines(model2$time,-log(model2$surv),type="s",col="blue")

legend("bottomright", legend=c("Placebo", "Traitement"), 
       col=c("red", "blue"), lty=c(1,1), title = "Groupes")


### Comparaison de fonctions de survie

test.survie = survdiff(survie~rx)
test.survie

ggsurvplot(fit2,
           conf.int=TRUE, # ajoutes les IC
           pval=TRUE, 	   # donne la p-value du test de log-rank
           risk.table=TRUE, #  tableau de risques sous le graphique
           legend.labs=c("Placebo", "Traitement"), # labels des groupes
           legend.title="Impact traitement", 
           palette=c("dodgerblue4", "orchid2"), # couleurs
           title="Fonctions de survie de Kaplan-Meier", 
           risk.table.height=.2)

ts = survdiff(survie~age2)
ts

ggsurvplot(fit3,
           conf.int=TRUE, # ajoutes les IC
           pval=TRUE, 	   # donne la p-value du test de log-rank
           risk.table=TRUE, #  tableau de risques sous le graphique
           legend.labs=c("c1","c2","c3"), # labels des groupes
           legend.title="Impact traitement", 
           palette=c("red", "green","orange"), # couleurs
           title="Fonctions de survie de Kaplan-Meier", 
           risk.table.height=.2)


####################################################

test.survie3 = survdiff(survie~resid.ds)
test.survie3

fit4 = survfit(survie~resid.ds,data=ovarian)  # avec stratification
fit4
summary(fit4)

ggsurvplot(fit4,
           conf.int=TRUE, # ajoutes les IC
           pval=TRUE, 	   # donne la p-value du test de log-rank
           risk.table=TRUE, #  tableau de risques sous le graphique
           legend.labs=c("Non malade", "Malade"), # labels des groupes
           legend.title="Impact de la maladie", 
           palette=c("dodgerblue4", "orchid2"), # couleurs
           title="Fonctions de survie de Kaplan-Meier", 
           risk.table.height=.2)

# Estimateur par noyau de convolution de la fonction de hasard

fr = muhaz(futime,fustat)
fr
plot(fr, xlab = "Temps",ylab = "h(t)",col=2,lwd=2)

# fonction de hasard par groupe

futime1=futime[rx==1] # placebo
length(futime1)
fustat1=fustat[rx==1]  # traitement

futime2=futime[rx==2]
length(futime2)
fustat2=fustat[rx==2]

fr1 = muhaz(futime1,fustat1)
fr2 = muhaz(futime2,fustat2)

plot(fr1, xlab = "Temps",ylab = "h(t)",col=1,lwd=2,
     main="fonction de risque des deux groupes")
lines(fr2, xlab = "Temps",ylab = "h(t)",col=2,lwd=2)
legend("topleft", legend=c("Placebo", "Traitement"), 
       col=c(1,2), lty=c(1,1), title = "Groupes")

# estimation du risque pour l'echelle
#performance du patient

futimec1=futime[ecog.ps=="1"]
fustatc1=fustat[ecog.ps=="1"]

futimec2=futime[ecog.ps=="2"]
fustatc2=fustat[ecog.ps=="2"]


frc1 = muhaz(futimec1,fustatc1)
frc2 = muhaz(futimec2,fustatc2)

plot(frc1, xlab = "Temps",ylab = "lambda1(t)",col=1,lwd=2,
     main="fonction de risque des deux groupes")

lines(frc2, xlab = "Temps",ylab = "lambda2(t)",col=2,lwd=2,
     main="fonction de risque des deux groupes")
legend("topleft", legend=c("Moyenne", "Eleve"), 
       col=c(1,2), lty=c(1,1), title = "Groupes")

# variable age
table(age2)

d1=futime[age2=="c1"]
s1=fustat[age2=="c1"]

d2=futime[age2=="c2"]
s2=fustat[age2=="c2"]

d3=futime[age2=="c3"]
s3=fustat[age2=="c3"]

frc1 = muhaz(d1,s1)
frc2 = muhaz(d2,s2)
frc3 = muhaz(d3,s3)


plot(frc1, xlab = "Temps",ylab = "lambda1(t)",col=1,lwd=2,
     main="fonction de risque des deux groupes")

lines(frc2, xlab = "Temps",ylab = "lambda2(t)",col=2,lwd=2,
      main="fonction de risque des deux groupes")
legend("topleft", legend=c("Moyenne", "Eleve"), 
       col=c(1,2), lty=c(1,1), title = "Groupes")
