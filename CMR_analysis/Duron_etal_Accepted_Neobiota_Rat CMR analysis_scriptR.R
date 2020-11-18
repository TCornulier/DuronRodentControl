### Duron et al. Neobiota ###
### CMR data analysis ###
### Nov 2020 ###
### Dataset CMR: "Duron_etal_Accepted_Neobiota_Rat CMR dataset.txt" ###
### Trap position: "Duron_etal_Accepted_Neobiota_Trap_position.txt" ###

################################################################################
# Library and dataset used
################################################################################

library(secr)
trap<- read.traps("Duron_etal_Accepted_Neobiota_Trap_position.txt")
plot(trap)

################################################################################
# Code to fix covariates
################################################################################

y<- read.capthist("Duron_etal_Accepted_Neobiota_Rat CMR dataset.txt","Duron_etal_Accepted_Neobiota_Trap_position.txt", detector="multi",fmt="trapID",covnames=c("removal","session","sessgr","species","sex","age","indgr"), verify="TRUE")
for (i in 1:length(y))
  covariates(y[[i]])$species <- factor(covariates(y[[i]])$species,
                                       levels = c('RR','RE'))
lapply(covariates(y), function(x) levels(x$species))
for (i in 1:length(y))
  covariates(y[[i]])$age <- factor(covariates(y[[i]])$age,
                                   levels = c('ad','juv'))
lapply(covariates(y), function(x) levels(x$age))
for (i in 1:length(y))
  covariates(y[[i]])$removal <- factor(covariates(y[[i]])$removal,
                                       levels = c('before','after'))
lapply(covariates(y), function(x) levels(x$removal))
for (i in 1:length(y))
  covariates(y[[i]])$session <- factor(covariates(y[[i]])$session,
                                       levels = c('winter','spring','summer'))
lapply(covariates(y), function(x) levels(x$session))
for (i in 1:length(y))
  covariates(y[[i]])$sessgr <- factor(covariates(y[[i]])$sessgr,
                                       levels = c('LowD','HighD'))
lapply(covariates(y), function(x) levels(x$sessgr))
for (i in 1:length(y))
  covariates(y[[i]])$indgr <- factor(covariates(y[[i]])$indgr,
                                       levels = c('RE','RRjuv','RRadM','RRadF'))
lapply(covariates(y), function(x) levels(x$indgr))

verify(y)
summary(y, terse=T)
summary(y)

################################################################################
# CMR models evaluated 
# with b: learned response to trapping; 
# species: R. exulans, R. rattus; 
# age: juveniles, adults; 
# sex: male, female; 
# indgr: R. exulans, R. rattus juveniles, R. rattus adult males, R. rattus adult females;
# session: the 6 sessions of CMR trapping; 
# sessgr: group of sessions 1-2-6 and group of sessions 3-4-5; 
# removal: effect of rat removal.
################################################################################

modelgobindgrsession_sigmaindgrsessgrCL<- secr.fit(y, model = list(D~1,g0~b+indgr+session,sigma~indgr+sessgr), CL=TRUE, buffer=60)###AICc=4124.429

modelgobindgrsessgr_sigmaindgrsessgrCL<- secr.fit(y, model = list(D~1,g0~b+indgr+sessgr,sigma~indgr+sessgr), CL=TRUE, buffer=60)###AICc=4143.143

modelgoindgrsession_sigmaindgrsessgrCL<- secr.fit(y, model = list(D~1,g0~indgr+session,sigma~indgr+sessgr), CL=TRUE, buffer=60)###AICc= 4122.769

modelgospeciessexagesession_sigmaindgrsessgrCL<- secr.fit(y, model = list(D~1,g0~species*sex*age+session,sigma~indgr+sessgr), CL=TRUE, buffer=60)###AICc=4127.628

modelgobspeciesagesession_sigmaindgrsessgrCL<- secr.fit(y, model = list(D~1,g0~b+species+age+session,sigma~indgr+sessgr), CL=TRUE, buffer=60)###AICc=4128.119

modelgospeciesagesession_sigmaindgrsessgrCL<- secr.fit(y, model = list(D~1,g0~species+age+session,sigma~indgr+sessgr), CL=TRUE, buffer=60)###AICc=4126.17

modelgospeciessession_sigmaindgrsessgrCL<- secr.fit(y, model = list(D~1,g0~species+session,sigma~indgr+sessgr), CL=TRUE, buffer=60)###AICc=4147.727

modelgoagesession_sigmaindgrsessgrCL<- secr.fit(y, model = list(D~1,g0~age+session,sigma~indgr+sessgr), CL=TRUE, buffer=60)###AICc=4124.014

modelgosession_sigmaindgrsessgrCL<- secr.fit(y, model = list(D~1,g0~session,sigma~indgr+sessgr), CL=TRUE, buffer=60)###AICc= 4145.788

modelgo1_sigmaindgrsessgrCL<- secr.fit(y, model = list(D~1,g0~1,sigma~indgr+sessgr), CL=TRUE, buffer=60)###AICc=4185.559

modelgoagesession_sigmaindgrCL<- secr.fit(y, model = list(D~1,g0~age+session,sigma~indgr), CL=TRUE, buffer=60)###AICc=4121.872

modelgoagesession_sigmasessgrCL<- secr.fit(y, model = list(D~1,g0~age+session,sigma~sessgr), CL=TRUE, buffer=60)###AICc= 4161.962

modelgoagesession_sigma1CL<- secr.fit(y, model = list(D~1,g0~age+session,sigma~1), CL=TRUE, buffer=60)###AICc= 4156.558

modelgoagesession_sigmasexCL<- secr.fit(y, model = list(D~1,g0~age+session,sigma~sex), CL=TRUE, buffer=60)###AICc=4139.841

modelgoindgrsession_sigmaindgrCL<- secr.fit(y, model = list(D~1,g0~indgr+session,sigma~indgr), CL=TRUE, buffer=60)###AICc= 4120.639

modelgoindgrsession_sigmaindgrremovalCL<- secr.fit(y, model = list(D~1,g0~indgr+session,sigma~indgr+removal), CL=TRUE, buffer=60)###AICc= 4121.056

modelgoindgrsession_sigmaremovalCL<- secr.fit(y, model = list(D~1,g0~indgr+session,sigma~removal), CL=TRUE, buffer=60)### AICc=4155.966

modelgoageETsexsession_sigmaindgrCL<- secr.fit(y, model = list(D~1,g0~age+sex+session,sigma~indgr), CL=TRUE, buffer=60)### AICc=4123.268

################################################################################
# Selection of the best model with AICc
################################################################################

TableAIC<- AIC(modelgobindgrsession_sigmaindgrsessgrCL,
  modelgobindgrsessgr_sigmaindgrsessgrCL,
  modelgoindgrsession_sigmaindgrsessgrCL,
  modelgospeciessexagesession_sigmaindgrsessgrCL,
  modelgobspeciesagesession_sigmaindgrsessgrCL,
  modelgospeciesagesession_sigmaindgrsessgrCL,
  modelgospeciessession_sigmaindgrsessgrCL,
  modelgoagesession_sigmaindgrsessgrCL,
  modelgosession_sigmaindgrsessgrCL,
  modelgo1_sigmaindgrsessgrCL,
  modelgoagesession_sigmaindgrCL, 
  modelgoagesession_sigmasessgrCL,
  modelgoagesession_sigma1CL,
  modelgoagesession_sigmasexCL,
  modelgoindgrsession_sigmaindgrCL,
  modelgoindgrsession_sigmaindgrremovalCL,
  modelgoindgrsession_sigmaremovalCL, 
  modelgoageETsexsession_sigmaindgrCL)

################################################################################
# Determination of Go and Sigma values from the best model : 
# modelgoindgrsession_sigmaindgrCL
################################################################################
print(modelgoindgrsession_sigmaindgrCL) 
# to obtain the results of the models including beta parameters

newdata2<- expand.grid(indgr= c("RE", "RRadF","RRadM","RRjuv"), session=c("1","2","3","4","5","6")) 
# function "expandgrid" to have go and sigma for all the categories and capture session

pred<- predict(modelgoindgrsession_sigmaindgrCL, newdata2)

################################################################################# Estimated rat density by individual groups and sessions
################################################################################

densindgrCL<- derived(modelgoindgrsession_sigmaindgrCL,groups="indgr")




