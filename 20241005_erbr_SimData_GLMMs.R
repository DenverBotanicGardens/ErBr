## April Goebl & Dan Doak
## Script modified 2024-10-01
## Collaboration with CU Boulder and Denver Botanic Gardens on Eriogonum brandegeii modeling 
## Run GLMM models on simulated data 


rm(list=ls())
graphics.off()

#install.packages("Matrix") 
#install.packages("TMB", type = "source")
#install.packages("glmmTMB", type = "source")
library(MASS)
library(dplyr)
library(matrixStats)
library(corrplot)
library(car)
library(coda)
library(lme4)
library(resample)
library(gplots)
library(stringr)
library(plotrix)
library(glmmTMB) # installed from tar for not the last version, but 1.1.9
library(TMB)
## ------------------------------------------------------------------------------------------------







## OBTAIN GLMM ESTIMATES ----------------------------------------------------------------

## Loop over datasets 
dateGLM <- as.character("20241005")
## ASSIGN NAME VARIABLE FOR DESIRED DATASETS 
name <- as.character("SimDat20yrMedGr")


#modList.grwth <- NULL        #List variable to store all model summaries
#modList.surv <- NULL        #List variable to store all models summaries
paramsMM.grwth <- NULL
seMM.grwth <- NULL    #Use SE in error bars in GLMM vs JAGS plots 
paramsMM.surv <- NULL
seMM.surv <- NULL
n.datset <- 10


for (dd in 1:n.datset) {
  
  ## Read in data
  noMiss <- readRDS(file=paste(dateGLM, "_erbr_", name, "NoMiss.srvCor.sdlgCor.grwthCor.",dd,".4GLM",".rds", sep=""))
  
  ## Lag the climate variables
  nn=length(noMiss$TempFall)
  noMiss$TempFall.lag	   =c(noMiss$TempFall[2:(nn)], NA)
  noMiss$TempSummer.lag	 =c(noMiss$TempSummer[2:(nn)], NA)
  noMiss$TempWinter.lag	 =c(noMiss$TempWinter[2:(nn)], NA)
  noMiss$PptFall.lag	   =c(noMiss$PptFall[2:(nn)], NA)
  noMiss$PptWinter.lag	 =c(noMiss$PptWinter[2:(nn)], NA)
  noMiss$PptSummer.lag	 =c(noMiss$PptSummer[2:(nn)], NA)
  
  ## Add t+1 climate, sz, & tag into erbr data 
  noMiss <- noMiss %>% mutate(TagNew1=lead(TagNew), RosNew1=lead(RosNew), Surv1=lead(surv))  
  #Remove lines with mis-matched individuals 
  noMiss <- noMiss[which(noMiss$TagNew == noMiss$TagNew1),]  
  
  ## Note: Log size in all models to match JAGS; this makes ending size a linear function of starting size
  
  plot(noMiss$RosNew,noMiss$RosNew1)
  print(t(table(noMiss$Year)))
  
  
  
  #lagged models ----------------------------------
  ## Growth (param order is same as for JAGS output)
  glmm.grwth <- glmmTMB(RosNew1 ~ log(RosNew) + TempFall.lag + TempSummer.lag + TempWinter.lag +     
                        PptFall.lag + PptSummer.lag + PptWinter.lag + (1|TransectNew), family=truncated_nbinom2(link="log"), data=noMiss)
  
  summary(glmm.grwth)
  
  
  ## Survival  (param order is same as for JAGS output)
  glmm.surv <- glmer(Surv1 ~ log(RosNew) + PptWinter.lag + TempFall.lag + TempSummer.lag + 
                       TempWinter.lag + (1|TransectNew), family=binomial(link='logit'), data=noMiss)
  
  summary(glmm.surv)
  ## ---------------------------------------------
  
  
  
  ## Save model summaries in list
  #modList.grwth[[length(modList.grwth) + 1]] <- summary(glmm.grwth)
  #modList.surv[[length(modList.surv) + 1]] <- summary(glmm.surv)
  
  
  ## Extract parameter estimates, including intercept, & SEs from GLMMs
  ## Growth
  paramsMM.grwthTmp <- fixef(glmm.grwth)
  seMM.grwthTmp <- as.data.frame(summary.sdreport(glmm.grwth$sdr))$`Std. Error`[1:8]
  
  paramsMM.grwth <- as.data.frame(rbind(paramsMM.grwth, paramsMM.grwthTmp$cond))
  seMM.grwth <- as.data.frame(rbind(seMM.grwth, seMM.grwthTmp))
  
  
  ## Survival
  paramsMM.survTmp <- as.data.frame(summary(glmm.surv)$coefficients[1:nrow(summary(glmm.surv)$coefficients),1])
  seMM.survTmp <- as.data.frame(summary(glmm.surv)$coefficients[1:nrow(summary(glmm.surv)$coefficients),2])
  
  colnames(paramsMM.survTmp) <- paste("GLMM.",dd,sep="")
  colnames(seMM.survTmp) <- paste("SE.",dd,sep="")
  paramsMM.surv <- as.data.frame(c(paramsMM.surv, paramsMM.survTmp))
  seMM.surv <- as.data.frame(c(seMM.surv, seMM.survTmp))

}
## ------ END GLMM LOOP -------------------------------------------------------





## RE-FORMAT AS NEEDED AND ADD VARIABLE LABELS COLUMN -------------------------
paramsMM.grwth
paramsMM.grwthT <- as.data.frame(t(paramsMM.grwth))
colnames(paramsMM.grwthT) <- c("GLMM.1", "GLMM.2","GLMM.3","GLMM.4","GLMM.5","GLMM.6","GLMM.7","GLMM.8","GLMM.9","GLMM.10")
seMM.grwthT <- as.data.frame(t(seMM.grwth))
colnames(seMM.grwthT) <- c("GLMM.1", "GLMM.2","GLMM.3","GLMM.4","GLMM.5","GLMM.6","GLMM.7","GLMM.8","GLMM.9","GLMM.10")

paramsMM.surv
seMM.surv
names.paramTitles <- c("Grwth Intercept","Grwth Size","Grwth Fall Temp","Grwth Summer Temp","Grwth Winter Temp",
                       "Grwth Fall Precip","Grwth Summer Precip","Grwth Winter Precip","Surv Intercept","Surv Size",
                       "Surv Winter Precip","Surv Fall Temp","Surv Summer Temp","Surv Winter Temp")


paramsMM.grwthT$ParamTitle <- names.paramTitles[1:8]
seMM.grwthT$ParamTitle <- names.paramTitles[1:8]
paramsMM.surv$ParamTitle <- names.paramTitles[9:14]
seMM.surv$ParamTitle <- names.paramTitles[9:14]


## Save GLMM param and se results 
dateRes <- as.character("20241005")
name

## Save GLMM model summaries
#saveRDS(modList.grwth, file=paste(dateRes, "_erbr_GLMMtmbSummGrwth_", name, ".rds", sep=""))
#saveRDS(modList.surv, file=paste(dateRes, "_erbr_GLMMsummSurv_", name, ".rds", sep=""))

## Note: lag climate vars correction made 10-05 **
saveRDS(paramsMM.grwthT, file=paste(dateRes, "_erbr_paramMMtmbGrwth_", name, ".rds", sep=""))
saveRDS(seMM.grwthT, file=paste(dateRes, "_erbr_seMMtmbGrwth_", name, ".rds", sep=""))
saveRDS(paramsMM.surv, file=paste(dateRes, "_erbr_paramMMsurv_", name, ".rds", sep=""))
saveRDS(seMM.surv, file=paste(dateRes, "_erbr_seMMsurv_", name, ".rds", sep=""))
## --------------------------------------
