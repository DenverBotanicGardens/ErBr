## April Goebl
## Started 06-08-20
## R script to run dredge (on command line if needed)
## for model selection for each vital rate




rm(list=ls())

## Load packages --------------------------------------------------------------------
library(MuMIn)
library(lme4)
library(dplyr)
## ----------------------------------------------------------------------------------




## Load data ------------------------------------------------------------------------
erbr <- read.csv("erbr_tagClust_200708.csv", header=TRUE)
## ----------------------------------------------------------------------------------




## Check structure of, and modify variables as needed -------------------------------
str(erbr)
erbr$Year <- as.factor(erbr$Year)


## Add in modified size variables 
erbr <- erbr %>% mutate(RosClust.sq=RosClust^2, RosClust.log=log(RosClust))


## Scale climate variables
pvars <- c("PptFall","PptWinter","PptSummer","TempFall","TempWinter","TempSummer",
           "TmaxFall","TmaxWinter","TmaxSummer","TminFall","TminWinter","TminSummer",
           "PptFall1","PptWinter1","PptSummer1","TempFall1","TempWinter1","TempSummer1",
           "TmaxFall1","TmaxWinter1","TmaxSummer1","TminFall1","TminWinter1","TminSummer1")

erbr.sc <- erbr
erbr.sc[pvars] <- lapply(erbr.sc[pvars], scale)
## ----------------------------------------------------------------------------------






## SURVIVAL ---------------------------------------------------------------------------
## Scale additional, relevent predictors
erbr.surv <- erbr.sc %>% mutate(RosClust.log=scale(RosClust.log))

## Remove NA values from response variable since dredge won't allow
erbr.surv <- erbr.surv[which(is.na(erbr.surv$SurvClust1)==FALSE),]


## Full global model
global.surv <- glmer(SurvClust1 ~ RosClust.log * (PptFall1 + PptWinter1 + PptSummer1 + TempFall1 + 
                    TempWinter1 + TempSummer1 + TmaxFall1 + TmaxWinter1 + TmaxSummer1 + TminFall1 + 
                    TminWinter1 + TminSummer1) + (1|Year) + (1|TransectNew),
                    family=binomial(link='logit'), data=erbr.surv, na.action='na.fail')

## Model selection (subset rule allows only 1 temperature variable (either mean, max, or min) per season)
dredge.surv <- dredge(global.surv, fixed='RosClust.log', beta='none', trace=2,
                      subset= !(TempFall1 && TmaxFall1) && !(TempFall1 && TminFall1) && !(TmaxFall1 && TminFall1) 
                      && !(TempWinter1 && TmaxWinter1) && !(TempWinter1 && TminWinter1) && !(TmaxWinter1 && TminWinter1) 
                      && !(TempSummer1 && TmaxSummer1) && !(TempSummer1 && TminSummer1) && !(TmaxSummer1 && TminSummer1))


## Look at top models
top.surv <- get.models(dredge.surv, subset=c(1:20))
AICc(top.surv[[1]])
AICc(top.surv[[20]])


## Save output
out <- capture.output(top.surv)
cat(out, file="200708_topSurvOutput.txt", sep="\t")

saveRDS(top.surv, "200708_topSurv.rds")
saveRDS(dredge.surv, "200708_dredgeSurv.rds")
## ----------------------------------------------------------------------------------






## GROWTH ------------------------------------------------------------------------------
## Scale additional, relevent predictors
erbr.grwth <- erbr.sc %>% mutate(RosClust=scale(RosClust), RosClust.sq=scale(RosClust.sq))

## Remove NA values from response variable since dredge won't allow
erbr.grwth <- erbr.grwth[which(is.na(erbr.grwth$RosClust1.new)==FALSE),]


## Full global model ** check if size and size squared are coded correctly for intetactions with climate **
global.grwth <- glmer.nb(RosClust1.new ~ RosClust.sq + (RosClust * (PptFall1 + PptWinter1 + PptSummer1 + TempFall1 + 
                         TempWinter1 + TempSummer1 + TmaxFall1 + TmaxWinter1 + TmaxSummer1 + TminFall1 + TminWinter1 + 
                         TminSummer1)) + (1|Year) + (1|TransectNew), data=erbr.grwth, na.action='na.fail')

## Model selection (subset rule allows only 1 temperature variable (either mean, max, or min) per season)
dredge.grwth <- dredge(global.grwth, fixed=c('RosClust','RosClust.sq'), beta='none', trace=2,
                      subset= !(TempFall1 && TmaxFall1) && !(TempFall1 && TminFall1) && !(TmaxFall1 && TminFall1) 
                      && !(TempWinter1 && TmaxWinter1) && !(TempWinter1 && TminWinter1) && !(TmaxWinter1 && TminWinter1) 
                      && !(TempSummer1 && TmaxSummer1) && !(TempSummer1 && TminSummer1) && !(TmaxSummer1 && TminSummer1))


## Look at top models
top.grwth <- get.models(dredge.grwth, subset=c(1:20))
AICc(top.grwth[[1]])
AICc(top.grwth[[20]])


## Save output
out <- capture.output(top.grwth)
cat(out, file="200708_topGrwthOutput.txt", sep="\t")

saveRDS(top.grwth, "200708_topGrwth.rds")
saveRDS(dredge.grwth, "200708_dredgeGrwth.rds")
## ----------------------------------------------------------------------------------






## PROBABILITY OF REPRODUCTION ------------------------------------------------------
## Scale additional, relevent predictors
erbr.reproYesNo <- erbr.sc %>% mutate(RosClust.log=scale(RosClust.log))

## Remove NA values from response variable since dredge won't allow
erbr.reproYesNo <- erbr.reproYesNo[which(is.na(erbr.reproYesNo$InflClustYesNo)==FALSE),]


## Full global model
global.reproYesNo <- glmer(InflClustYesNo ~ RosClust.log * (PptFall + PptWinter + PptSummer + TempFall + TempWinter + TempSummer
                           + TmaxFall + TmaxWinter + TmaxSummer + TminFall + TminWinter + TminSummer) + (1|Year) + (1|TransectNew), 
                           family=binomial(link=logit), data=erbr.reproYesNo, na.action='na.fail')

## Model selection (subset rule allows only 1 temperature variable (either mean, max, or min) per season)
dredge.reproYesNo <- dredge(global.reproYesNo, fixed='RosClust.log', beta='none', trace=2,
                      subset= !(TempFall && TmaxFall) && !(TempFall && TminFall) && !(TmaxFall && TminFall) 
                      && !(TempWinter && TmaxWinter) && !(TempWinter && TminWinter) && !(TmaxWinter && TminWinter) 
                      && !(TempSummer && TmaxSummer) && !(TempSummer && TminSummer) && !(TmaxSummer && TminSummer))


## Look at top models
top.reproYesNo <- get.models(dredge.reproYesNo, subset=c(1:20))
AICc(top.reproYesNo[[1]])
AICc(top.reproYesNo[[20]])


## Save output
out <- capture.output(top.reproYesNo)
cat(out, file="200708_topReproYesNoOutput.txt", sep="\t")

saveRDS(top.reproYesNo, "200708_topReproYesNo.rds")
saveRDS(dredge.reproYesNo, "200708_dredgeReproYesNo.rds")
## ----------------------------------------------------------------------------------






## REPRODUCTION ---------------------------------------------------------------------
## Scale additional, relevent predictors
erbr.repro <- erbr.sc %>% mutate(RosClust=scale(RosClust))

## Use only reproductive plts
erbr.repro <- subset(erbr.repro, InflClust>0)

## Remove NA values from response variable since dredge won't allow
erbr.repro <- erbr.repro[which(is.na(erbr.repro$InflClust)==FALSE),]


## Full global model (subset rule allows only 1 temperature variable (either mean, max, or min) per season)
global.repro <- lmer(InflClust ~ RosClust * (PptFall + PptWinter + PptSummer + TempFall + TempWinter + TempSummer
                     + TmaxFall + TmaxWinter + TmaxSummer + TminFall + TminWinter + TminSummer) + (1|Year) + 
                     (1|TransectNew), REML=FALSE, data=erbr.repro, na.action='na.fail')

## Model selection
dredge.repro <- dredge(global.repro, fixed='RosClust', beta='none', trace=2,
                       subset= !(TempFall && TmaxFall) && !(TempFall && TminFall) && !(TmaxFall && TminFall) 
                       && !(TempWinter && TmaxWinter) && !(TempWinter && TminWinter) && !(TmaxWinter && TminWinter) 
                       && !(TempSummer && TmaxSummer) && !(TempSummer && TminSummer) && !(TmaxSummer && TminSummer))


## Look at top models
top.repro <- get.models(dredge.repro, subset=c(1:20))
AICc(top.repro[[1]])
AICc(top.repro[[20]])


## Save output
out <- capture.output(top.repro)
cat(out, file="200708_topReproOutput.txt", sep="\t")

saveRDS(top.repro, "200708_topRepro.rds")
saveRDS(dredge.repro, "200708_dredgeRepro.rds")
## ----------------------------------------------------------------------------------






## VARIANCE IN GROWTH ---------------------------------------------------------------

## ** Add model selection for variance in growth ... *** 

