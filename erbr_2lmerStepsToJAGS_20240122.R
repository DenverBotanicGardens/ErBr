## Dan Doak & April Goebl 
## Script modified 20-05-04
## Collaboration with Denver Botanic Gardens
## Modify data and assign variables needed for JAGS model with data lags (missing years)
## Associated JAGS script models growth, survival, repro, and recruitment 
## Use run.jags to run associated JAGS script
## Note: Details on mixed models given here: https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
## Other relevant links:
# https://bayesball.github.io/BOOK/bayesian-multiple-regression-and-logistic-models.html#bayesian-logistic-regression
# http://www.columbia.edu/~cjd11/charles_dimaggio/DIRE/styled-4/styled-11/code-8/#logistic-regression



rm(list=ls())
graphics.off()


## LOAD PACKAGES AND FUNCTIONS --------------------------------------------------------------------
library(lme4)
library(ggplot2)
library(rjags)
library(runjags)
library(dplyr)
library(coda)
library(corrplot)
library(robustbase)
library(resample)
library(gplots)
library(matrixStats)
## ------------------------------------------------------------------------------------------------




## LOAD DATA --------------------------------------------------------------------------------------
dats <- read.csv("erbr_TagClust2022_20230408.csv", header = TRUE)
#dats <- read.csv("erbr_TagClust4to8_210617.csv", header = TRUE)
## ------------------------------------------------------------------------------------------------



## SET WD (WHERE JAGS SCRIPT IS LOCATED) ----------------------------------------------------------
setwd("C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship")
#setwd("~/ErBr")
## ------------------------------------------------------------------------------------------------



## CHANGE NAMES OF CLIMATE COLUMNS ------------------------------------------------------------------
dats <- rename(dats, PptFall=Tot_fall_ppt, PptWinter=Tot_winter_ppt, PptSummer=Tot_summer_ppt,
               TempFall=Mean_fall_temp, TempWinter=Mean_winter_temp, TempSummer=Mean_summer_temp)
## ------------------------------------------------------------------------------------------------




## SETTING UP NECESSARY VARIABLES -----------------------------------------------------------------
## Setting up the jags model with lagged values
Nallrows <- length(dats$Site)

numyears <- length(unique(dats$Year))
numtrans <- length(unique(dats$TransectNew))

## Identify rows that are good dependent values (ending sizes) for surv or growth  
lagforsurv <- dats$lagforsurv       #Full list of lags or of -1 for first observation rows
goodrows <- which(dats$lagforsurv > 0)
goodgrowrows <- which(dats$lagsrtsz > 0)

## Identify the lag for these rows: how far back is the last good size measurement?
lagvals <- dats$lagforsurv[which(dats$lagforsurv > 0)]
Ncases <- length(goodrows)
Ngrowcases <- length(goodgrowrows)
Survs <- dats$surv
RosNew <- dats$RosNew
InflNew <-  dats$InflNew 
InflYesNo <- dats$InflYesNo


## Setting up variables for use in repro fitting: 
dats$InflYesNo <- dats$InflNew
dats$InflYesNo[dats$InflNew>1] <- 1
rows.w.sz <- which(is.na(dats$RosNew)==FALSE)
rows.wo.sz <- which(is.na(dats$RosNew)==TRUE)
Ndirectszcases <- length(rows.w.sz)           #Direct measures of sz upon which to base repro
Nindirectszcases <- length(rows.wo.sz)        #No direct measures of sz; repro to be inferred from estimated sz 
rows.w.inflors <- which(dats$InflNew>0)       #Non-zero estimates of infs so repro amt can be estimated, if reproductive
Nrows.w.inflors <- length(rows.w.inflors)


## Add vector to indicate if alive or dead after missing yr(s)
dats$RowNum <- 1:nrow(dats)                           #Add a column to indicate row number
rows.wo.sz.alive <- as.data.frame(matrix(NA, nrow=length(rows.wo.sz), ncol=2))
colnames(rows.wo.sz.alive) <- c("Rows", "Alive")
rows.wo.sz.alive$Rows <- rows.wo.sz
 
for (ww in rows.wo.sz) {                                                  #Loop over all tags with 1 or more missing yrs                      
  tag.val <- dats$TagNew[ww]
  tag.each <- subset(dats, dats$TagNew==tag.val)                          #Process each tag 
  tag.surv <- tag.each$surv[!is.na(tag.each$surv) & tag.each$RowNum>ww]   #Store surv for 1st non-missing yr post each missed yr
  rows.wo.sz.alive$Alive[rows.wo.sz.alive$Rows==ww] <- tag.surv[1] 
}

rows.wo.sz.alive$Alive[is.na(rows.wo.sz.alive$Alive)] <- 0  #Change NAs to 0, these are lines where missed yr was last & recorded as dead
rows.wo.sz.alive <- rows.wo.sz.alive$Alive                  #Change to vector



## Make transect & year values numerical to use in jags as random effects 
dats$TransectNew.num <- as.factor(dats$TransectNew)
dats$TransectNew.num <- as.numeric(dats$TransectNew.num)
TransectNew.num <- dats$TransectNew.num

dats$Year.num <- as.factor(dats$Year)
dats$Year.num <- as.numeric(dats$Year.num)   
Year.num <- dats$Year.num

## Make a linear index of transect-year combos
yrtranscombo=100*dats$TransectNew.num+dats$Year.num

## Set up a logical variable that is whether there is surv or grwth data in a yr (0,1): this is needed for the summing up of infl nums to predict new plts
datayesno <- rep(1,length(dats$surv))
datayesno[which(is.na(dats$surv)==TRUE)] <- 0 
## ------------------------------------------------------------------------------------------------



## Make dataframe w new plts (that are likely recent seedlings) for each transect & yr ------------
## Make df that will hold data containing new plants
years <- unique(dats$Year.num)
years <- years[order(years)]
dats.newPlts <- as.data.frame(rep(unique(dats$TransectNew.num), each=length(years)))
colnames(dats.newPlts) <- "TransectNew.num"
dats.newPlts$Year.num <- rep(years)

## Identify new plants
newPlts <- dats %>% group_by(TagNew) %>% slice(which.min(Year))   #Identify rows with 1st appearance for each plt
newPlts <- newPlts[newPlts$Year!=2004,]                           #Remove 2004 (first year of data collection)
sz.cutoff <- 5                                                    #Sz cutoff, above which plt was likely not recently a seedling 
newPlts <- newPlts[newPlts$RosNew < sz.cutoff,]                   #Remove if >X rosettes (these were likely missed and are not new)
num.newPlts <- newPlts %>% group_by(TransectNew.num, Year.num) %>% summarise(num.newPlts=n())  #Count num new plts per yr & transect

## Add number of new plants to df of each transect & year
dats.newPlts <- left_join(dats.newPlts, num.newPlts, by=c("TransectNew.num", "Year.num"))
dats.newPlts$num.newPlts[is.na(dats.newPlts$num.newPlts)] <- 0   #Change NAs (no new plants) to zeros
dats.newPlts$num.newPlts[dats.newPlts$Year.num==1] <- NA         #Change new plts in 2004 (yr 1) to NA

## Add column so new plts in t+1 match year t
dats.newPlts <- dats.newPlts %>% mutate(num.newPlts1=lead(num.newPlts))  

## Add climate variables to new plants data 
dats.clim <- dats %>% dplyr::select(c(Year.num, PptFall, PptWinter, PptSummer, TempFall, TempWinter, TempSummer)) 
clim <- unique(dats.clim)
dats.newPlts <- dplyr::left_join(dats.newPlts, clim, by="Year.num")
dats.newPlts <- dats.newPlts %>% dplyr::mutate(PptFall1=lead(PptFall), PptWinter1=lead(PptWinter), PptSummer1=lead(PptSummer),
                                        TempFall1=lead(TempFall), TempWinter1=lead(TempWinter), TempSummer1=lead(TempSummer))  
## Note: variables with '1' at the end (e.g. num.newPlts1, PptFall1) represent t+1 data   


## Remove lines that correspond to transect-year combos that are not in the main data file
## Note: GP East Transect 6 & 7 do not have data from 2004, 2005, or 2006. 
dats.newPlts$yrtranscombo=100*dats.newPlts$TransectNew.num+dats.newPlts$Year.num
yrtrans.unq <- unique(yrtranscombo)
dats.newPlts <- dats.newPlts[dats.newPlts$yrtranscombo %in% yrtrans.unq,]

dats.newPlts <- dats.newPlts[is.na(dats.newPlts$num.newPlts1) == FALSE,]

newplts <- dats.newPlts$num.newPlts1
newplt.trans <- dats.newPlts$TransectNew.num
newplt.yr <- dats.newPlts$Year.num

newPltlines <- length(dats.newPlts$TransectNew.num)

## Make a linear index of transect-year combos for new plts
newplt.yrtranscombo=100*newplt.trans+newplt.yr 
## ------------------------------------------------------------------------------------------------


## 20230420 In response to error, make sure any rows with data for RosNew don't have NA for InflNew and InflYesNo
## For future, make change in erbr_1ReformatData_forJAGS script
dats[1217,]$InflNew <- 0
dats[1217,]$InflYesNo <- 0




## RUN ASSOCIATED JAGS MODEL ----------------------------------------------------------------------
#jags.mod <- run.jags('erbr_JAGSmod_210507.R', n.chains=3, data=dats, burnin=10000, thin=10, sample=30000, adapt=500, method='parallel')
#jags.mod <- run.jags('erbr_JAGSmodComplx_noYRE_210827.R', n.chains=3, data=dats, burnin=10000, thin=10, sample=30000, adapt=500, method='parallel')
#jags.mod <- run.jags('Scripts/erbr_JAGSmodBest_noYRE_20230418.R', n.chains=3, data=dats, burnin=10000, thin=10, sample=30000, adapt=500, method='parallel')

#save(jags.mod, file='erbr_JAGSmod_c3t10s20b5_210406.rdata')
#saveRDS(jags.mod, "erbr_JAGSmodBest_c3t10s30b10_noYRE_20230420.rds")
## ------------------------------------------------------------------------------------------------


## LOOK AT MODEL OUTPUT ---------------------------------------------------------------------------
#jags.modNYE <- jags.mod
jags.mod <- readRDS("Results_data_ms_archive/erbr_JAGSmod_c3t10s30b10_noYRE_4to13even_210617.rds")
summary(jags.mod)
#plot(jags.mod)
summ.mod <- summary(jags.mod)
tail(summ.mod[,1:3], n=31)
gelman.diag(jags.mod, confidence = 0.95, transform=FALSE)


## Compare median param estimates b/w models (e.g. with and without RE)
#summ.modNYE <- summary(jags.modNYE)

#summNYE <- tail(summ.modNYE, n=21)
#summ <- tail(summ.mod, n=21)
#prcnt_diff <- ((summ[,2]-summNYE[,2])/summ[,2])*100


## Calculate median and SD of param estimates
#chains <- jags.mod$mcmc
#chains <- bind_rows(lapply(chains, as.data.frame))
#colMeds <- apply(chains,2,median)
#colSDs <- apply(chains,2,sd)


## Look at correlation b/w params 
#chains.1 <- chains %>% dplyr::select(!contains(c("randomeffect", "precision")))
#chains.1 <- chains.1 %>% dplyr::select(!c(deviance, resid.sum.sq))

#cor.chains <- cor(chains.1)
#corrplot(cor.chains, method="circle", type="lower")
## -----------------------------------------------------------------------------------------------






## Make plot comparing median param ests b/w diff datasets ---------------------------------
## Show median param ests b/w diff datasets as points and 90% limits; include glmm estimates

## Load data from previous JAGS runs
#setwd("C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship")
chains <- readRDS("chains.c3t10s30b10_noYRE_20230420.rds")
chains.4to13 <- readRDS("chains.c3t10s30b10_noYRE_4to13_210616.rds")
chains.4to13evn <- readRDS("chains.c3t10s30b10_noYRE_4to13even_210617.rds")
chains.4to13odd <- readRDS("chains.c3t10s30b10_noYRE_4to13odd_210621.rds")
chains.4to8 <- readRDS("chains.c3t10s30b10_noYRE_4to8_210701.rds")
chains.9to13 <- readRDS("chains.c3t10s30b10_noYRE_9to13_210623.rds")



## CALCULATE MEDIAN PARAMETER VALUES 
medParams <- as.data.frame(colMedians(as.matrix(chains)))
medParams.4to13 <- as.data.frame(colMedians(as.matrix(chains.4to13)))
medParams.4to13evn <- as.data.frame(colMedians(as.matrix(chains.4to13evn)))
medParams.4to13odd <- as.data.frame(colMedians(as.matrix(chains.4to13odd)))
medParams.4to8 <- as.data.frame(colMedians(as.matrix(chains.4to8)))
medParams.9to13 <- as.data.frame(colMedians(as.matrix(chains.9to13)))

medComb <- as.data.frame(cbind(medParams, medParams.4to13, medParams.4to13evn, medParams.4to13odd, 
                               medParams.4to8, medParams.9to13))
colnames(medComb) <- c("Full", "4to13", "4to13evn", "4to13odd", "4to8", "9to13")
medComb$Names <- colnames(chains)




## CALCULATE 90th PERCENTILE ON JAGS ESTIMATES ------------------
quantParams <- chains %>% summarise_all(funs(list(quantile(., probs=0.9)))) #%>% transpose
quantParams <- as.data.frame(t(quantParams))
quantParams.4to13 <- chains.4to13 %>% summarise_all(funs(list(quantile(., probs=0.9))))
quantParams.4to13 <- as.data.frame(t(quantParams.4to13))
quantParams.4to13evn <- chains.4to13evn %>% summarise_all(funs(list(quantile(., probs=0.9))))
quantParams.4to13evn <- as.data.frame(t(quantParams.4to13evn))
quantParams.4to13odd <- chains.4to13odd %>% summarise_all(funs(list(quantile(., probs=0.9))))
quantParams.4to13odd <- as.data.frame(t(quantParams.4to13odd))
quantParams.4to8 <- chains.4to8 %>% summarise_all(funs(list(quantile(., probs=0.9))))
quantParams.4to8 <- as.data.frame(t(quantParams.4to8))
quantParams.9to13 <- chains.9to13 %>% summarise_all(funs(list(quantile(., probs=0.9))))
quantParams.9to13 <- as.data.frame(t(quantParams.9to13))

quantComb <- as.data.frame(cbind(quantParams, quantParams.4to13, quantParams.4to13evn, quantParams.4to13odd, 
                                 quantParams.4to8, quantParams.9to13))
colnames(quantComb) <- c("Full", "4to13", "4to13evn", "4to13odd", "4to8", "9to13")
quantComb$Names <- colnames(chains)


## CALCULATE 10th PERCENTILE ON JAGS ESTIMATES
quant10Params <- chains %>% summarise_all(funs(list(quantile(., probs=0.1)))) 
quant10Params <- as.data.frame(t(quant10Params))
quant10Params.4to13 <- chains.4to13 %>% summarise_all(funs(list(quantile(., probs=0.1))))
quant10Params.4to13 <- as.data.frame(t(quant10Params.4to13))
quant10Params.4to13evn <- chains.4to13evn %>% summarise_all(funs(list(quantile(., probs=0.1))))
quant10Params.4to13evn <- as.data.frame(t(quant10Params.4to13evn))
quant10Params.4to13odd <- chains.4to13odd %>% summarise_all(funs(list(quantile(., probs=0.1))))
quant10Params.4to13odd <- as.data.frame(t(quant10Params.4to13odd))
quant10Params.4to8 <- chains.4to8 %>% summarise_all(funs(list(quantile(., probs=0.1))))
quant10Params.4to8 <- as.data.frame(t(quant10Params.4to8))
quant10Params.9to13 <- chains.9to13 %>% summarise_all(funs(list(quantile(., probs=0.1))))
quant10Params.9to13 <- as.data.frame(t(quant10Params.9to13))

quant10Comb <- as.data.frame(cbind(quant10Params, quant10Params.4to13, quant10Params.4to13evn, quant10Params.4to13odd, 
                                   quant10Params.4to8, quant10Params.9to13))
colnames(quant10Comb) <- c("Full", "4to13", "4to13evn", "4to13odd", "4to8", "9to13")
quant10Comb$Names <- colnames(chains)
## ----------------------------------------------------------------------------



## SUBSET OUTPUT TO JUST KEEP PARAMS OF INTEREST
names.param <- colnames(chains)[55:85]
names.param <- names.param[c(2:8,12:16,18:23,25:30)] #Remove intercept plots 
names.paramTitles <- c("Grwth Size","Grwth Fall Temp","Grwth Summer Temp","Grwth Winter Temp",
                       "Grwth Fall Precip","Grwth Summer Precip","Grwth Winter Precip","Surv Size",
                       "Surv Winter Precip","Surv Fall Temp","Surv Summer Temp","Surv Winter Temp",
                       "p(Repro) Size","p(Repro) Fall Precip","p(Repro) Summer Precip",
                       "p(Repro) Fall Temp","p(Repro) Summer Temp","p(Repro) Winter Temp",
                       "Repro Size","Repro Fall Precip","Repro Summer Precip",
                       "Repro Winter Temp","Repro Fall Temp","Repro Summer Temp")

## Re-order parameter names for plotting 
index<-c(1,9,5,13,20,17,23,2,24,10,6,14,3,21,18,11,7,15,4,22,19,16,12,8)
names.paramOrd <- names.param[order(index)]
names.paramTitlesOrd <- names.paramTitles[order(index)]
## ------------------------------------------------------------------------------------------------




## USE A GLMM FOR EACH VR AND COMPARE PARAM ESTIMATES TO JAGS MODELS ---------------
## Add t+1 climate, sz, & tag into erbr data 
dats <- dats %>% mutate(TagNew1=lead(TagNew), RosNew1=lead(RosNew), Surv1=lead(surv))  
dats <- dats %>% mutate(PptFall1=lead(PptFall), PptWinter1=lead(PptWinter), PptSummer1=lead(PptSummer),
                        TempFall1=lead(TempFall), TempWinter1=lead(TempWinter), TempSummer1=lead(TempSummer))

dats <- dats[which(dats$TagNew == dats$TagNew1),]  #Remove lines with mis-matched individuals 

## Log size in all models to match JAGs; this makes ending size a linear function of starting size

## Growth
glmm.grwth <- glmer.nb(RosNew1 ~ log(RosNew) + PptFall + PptWinter + PptSummer + 
                           TempFall + TempWinter + TempSummer + (1|TransectNew), data=dats)

## Survival  
glmm.surv <- glmer(Surv1 ~ log(RosNew) + PptWinter + TempFall + (TempWinter) + 
                       (TempSummer) + (1|TransectNew), family=binomial(link='logit'), data=dats)

## Probability of reproduction
glmm.reproYesNo <- glmer(InflYesNo ~ log(RosNew) + (PptFall) + (PptSummer) + (TempFall) + 
                             (TempWinter) + (TempSummer) + (1|TransectNew), family=binomial(link=logit), data=dats)

## Reproduction 
## ** Modify inf data to only include reproductive plts (i.e. infs>0) in this model *
dats$InflNewMod <- dats$InflNew
dats$InflNewMod[dats$InflNew == 0] <- NA

glmm.reproOld <- glmer.nb(InflNew ~ log(RosNew) + (PptFall) + (PptSummer) + (TempFall) + 
                           (TempWinter) + (TempSummer) + (1|TransectNew), data=dats)
glmm.repro <- glmer.nb(InflNewMod ~ log(RosNew) + (PptFall) + (PptSummer) + (TempFall) + 
                         (TempWinter) + (TempSummer) + (1|TransectNew), data=dats)


## Extract parameter estimates and SEs from GLMMs
paramsMM.grwth <- as.data.frame(cbind(as.data.frame(summary(glmm.grwth)$coefficients[2:nrow(summary(glmm.grwth)$coefficients),1:2]),
                        c("Grwth Size","Grwth Fall Precip","Grwth Winter Precip","Grwth Summer Precip",
                          "Grwth Fall Temp","Grwth Winter Temp","Grwth Summer Temp")))
colnames(paramsMM.grwth) <- c("GLMM","SE","ParamTitle")
paramsMM.surv <- as.data.frame(cbind(as.data.frame(summary(glmm.surv)$coefficients[2:nrow(summary(glmm.surv)$coefficients),1:2]),
                                      c("Surv Size","Surv Winter Precip",
                                        "Surv Fall Temp","Surv Winter Temp","Surv Summer Temp")))
colnames(paramsMM.surv) <- c("GLMM","SE","ParamTitle")
paramsMM.reproYesNo <- as.data.frame(cbind(as.data.frame(summary(glmm.reproYesNo)$coefficients[2:nrow(summary(glmm.reproYesNo)$coefficients),1:2]),
                                     c("p(Repro) Size","p(Repro) Fall Precip","p(Repro) Summer Precip",
                                       "p(Repro) Fall Temp","p(Repro) Winter Temp","p(Repro) Summer Temp")))
colnames(paramsMM.reproYesNo) <- c("GLMM","SE","ParamTitle")
paramsMM.repro <- as.data.frame(cbind(as.data.frame(summary(glmm.repro)$coefficients[2:nrow(summary(glmm.repro)$coefficients),1:2]),
                                           c("Repro Size","Repro Fall Precip","Repro Summer Precip",
                                             "Repro Fall Temp","Repro Winter Temp","Repro Summer Temp")))
colnames(paramsMM.repro) <- c("GLMM","SE","ParamTitle")
paramsMM <- rbind(paramsMM.grwth, paramsMM.surv, paramsMM.reproYesNo, paramsMM.repro)
paramsMM$SE_upr <- paramsMM$GLMM + paramsMM$SE
paramsMM$SE_lwr <- paramsMM$GLMM - paramsMM$SE
## --------------------------------------------------------------




## Make dataframe of selected relevant values for plotting 
medComb.sel <- NULL
for (nn in 1:length(names.paramOrd)) {
  medComb.sel <- rbind(medComb.sel, as.data.frame(medComb[which(medComb$Names == names.paramOrd[nn]),1:6]))
}

## Combine JAGS and GLMM estimates
medComb.sel$ParamTitle <- names.paramTitlesOrd
medComb.sel <- dplyr::left_join(medComb.sel, paramsMM, by="ParamTitle")
medComb.sel <- medComb.sel %>% relocate(ParamTitle, .after=last_col())

quantComb.sel <- NULL
for (uu in 1:length(names.paramOrd)) {
  quantComb.sel <- rbind(quantComb.sel, as.data.frame(quantComb[which(quantComb$Names == names.paramOrd[uu]),1:6]))
}
quantComb.sel$GLMM_SEupr <- medComb.sel$SE_upr

quant10Comb.sel <- NULL
for (uu in 1:length(names.paramOrd)) {
  quant10Comb.sel <- rbind(quant10Comb.sel, as.data.frame(quant10Comb[which(quant10Comb$Names == names.paramOrd[uu]),1:6]))
}
quant10Comb.sel$GLMM_SElwr <- medComb.sel$SE_lwr



## Calculate the min and max values of error bars & set appropriate y-axis values for plotting 
# Change structure to allow for row min and max calcs 
quantComb.selMod <- cbind(as.numeric(quantComb.sel[,1]),as.numeric(quantComb.sel[,2]),
                           as.numeric(quantComb.sel[,3]),as.numeric(quantComb.sel[,4]),
                           as.numeric(quantComb.sel[,5]),as.numeric(quantComb.sel[,6]),
                           as.numeric(quantComb.sel[,7]))
quant.max <- rowMaxs(as.matrix(quantComb.selMod[,c(1:7)]))
yMax <- quant.max + (quant.max*0.05)

quant10Comb.selMod <- cbind(as.numeric(quant10Comb.sel[,1]),as.numeric(quant10Comb.sel[,2]),
                            as.numeric(quant10Comb.sel[,3]),as.numeric(quant10Comb.sel[,4]),
                            as.numeric(quant10Comb.sel[,5]),as.numeric(quant10Comb.sel[,6]),
                            as.numeric(quant10Comb.sel[,7]))
quant10.min <- rowMins(as.matrix(quant10Comb.selMod[,c(1:7)]))
#Since some min values are positive (all sz ests) and some are negative (all the rest), do separately
yMin <- c(quant10.min[1:4] - (quant10.min[1:4]*0.05), quant10.min[5:24] - (quant10.min[5:24]*-0.05))





## Assign colors
colz <- c("#d73027","#fc8d59","#91bfdb","#4575b4","#fdcb44","#fee090","grey60")


#saveRDS(medComb.sel, "20240117_medCombSel.rds")
#saveRDS(quantComb.sel, "20240117_quantCombSel.rds")
#saveRDS(quant10Comb.sel, "20240117_quant10CombSel.rds")
#medComb.sel <- readRDS("20240117_medCombSel.rds")
#quantComb.sel <- readRDS("20240117_quantCombSel.rds")
#quant10Comb.sel <- readRDS("20240117_quant10CombSel.rds")



## PLOT
#tiff('20230122_ErBr_fig2.tiff', res=400, pointsize=6, compression="lzw")
pdf('20240305_ErBr_fig2.pdf', width=6.7, height=9)
par(mfrow=c(7,4), mar=c(1.25,2,1.9,2))  #Plot so 4 VR models are in cols and upto 7 predictor vars are in rows
#bottom, left, top, and right
for (nn in 1:17) {
  plot(c(1:7), medComb.sel[nn,1:7], col=colz, ylab=NA,
       xaxt = "n", main=names.paramTitlesOrd[nn], cex.axis=0.9,cex.main=0.9, pch=19, 
       ylim=c(yMin[nn],yMax[nn]), cex=1.2)
  abline(h=0, col="grey80")
  arrows(1:6,as.numeric(medComb.sel[nn,1:6]),
         1:6, as.numeric(as.matrix(quantComb.sel[nn,1:6])),
         lwd = 1.25, angle = 90, code = 3, length=0, col=colz)
  arrows(1:6, as.numeric(medComb.sel[nn,1:6]),
         1:6, as.numeric(as.matrix(quant10Comb.sel[nn,1:6])), 
         lwd = 1.25, angle = 90, code = 3, length=0, col=colz)
}
plot.new()
for (nn in 18:20) {
  plot(c(1:7), medComb.sel[nn,1:7], col=colz, ylab=NA, 
       xaxt = "n", main=names.paramTitlesOrd[nn], cex.axis=0.9,cex.main=0.9, pch=19, 
       ylim=c(yMin[nn],yMax[nn]), cex=1.2)
  abline(h=0, col="grey80")
  arrows(1:7,as.numeric(medComb.sel[nn,1:7]),
         1:7, as.numeric(as.matrix(quantComb.sel[nn,1:7])),
         lwd = 1.25, angle = 90, code = 3, length=0, col=colz)
  arrows(1:7, as.numeric(medComb.sel[nn,1:7]),
         1:7, as.numeric(as.matrix(quant10Comb.sel[nn,1:7])), 
         lwd = 1.25, angle = 90, code = 3, length=0, col=colz)
}
plot.new()
for (nn in 21:24) {
  plot(c(1:7), medComb.sel[nn,1:7], col=colz, ylab=NA,
       xaxt = "n", main=names.paramTitlesOrd[nn], cex.axis=0.9,cex.main=0.9, pch=19, 
       ylim=c(yMin[nn],yMax[nn]), cex=1.2)
  abline(h=0, col="grey80")
  arrows(1:7,as.numeric(medComb.sel[nn,1:7]),
         1:7, as.numeric(as.matrix(quantComb.sel[nn,1:7])),
         lwd = 1.25, angle = 90, code = 3, length=0, col=colz)
  arrows(1:7, as.numeric(medComb.sel[nn,1:7]),
         1:7, as.numeric(as.matrix(quant10Comb.sel[nn,1:7])), 
         lwd = 1.25, angle = 90, code = 3, length=0, col=colz)
}
plot.new()

#legend("center", colnames(medComb.sel)[1:7], col=colz, pch=19, cex=1,
#       horiz=FALSE, bty="y",seg.len=1, xpd="NA")
legend("center", c("Full dataset","2004-2013","2004-2013 even yrs only","2004-2013 odd yrs only",
                   "2004-2008","2009-2013","GLMM"), col=colz, pch=19, cex=1,
       horiz=FALSE, bty="y",seg.len=1, xpd="NA")

dev.off()
## -----------------------------------------------------------------





## MAKE FIGURE SHOWING ESTIMATES FOR TRANSECT RANDOM EFFECTS FOR EACH MODEL AND DATASET
## Use code to calculate medians and quntiles and GLMMs above 

## SUBSET OUTPUT TO JUST KEEP PARAMS OF INTEREST (I.E. TRANSECT RANDOM EFFECT ESTIMATES)
namesRE <- colnames(chains)[3:50]


## Make dataframe of selected relevant values for plotting 
medComb.selRE <- NULL
for (nn in 1:length(namesRE)) {
  medComb.selRE <- rbind(medComb.selRE, as.data.frame(medComb[which(medComb$Names==namesRE[nn]),1:6]))
}

## Combine JAGS and GLMM estimates
#medComb.sel$ParamTitle <- names.paramTitlesOrd
#medComb.sel <- dplyr::left_join(medComb.sel, paramsMM, by="ParamTitle")
#medComb.sel <- medComb.sel %>% relocate(ParamTitle, .after=last_col())

quantComb.selRE <- NULL
for (uu in 1:length(namesRE)) {
  quantComb.selRE <- rbind(quantComb.selRE, as.data.frame(quantComb[which(quantComb$Names==namesRE[uu]),1:6]))
}
#quantComb.selRe$GLMM_SEupr <- medComb.selRe$SE_upr

quant10Comb.selRE <- NULL
for (uu in 1:length(namesRE)) {
  quant10Comb.selRE <- rbind(quant10Comb.selRE, as.data.frame(quant10Comb[which(quant10Comb$Names==namesRE[uu]),1:6]))
}
#quant10Comb.selRe$GLMM_SElwr <- medComb.selRe$SE_lwr
## ---------------------------------------------



## Calculate the min and max values of error bars & set appropriate y-axis values for plotting 
# Change structure to allow for row min and max calcs 
quantComb.selModRE <- cbind(as.numeric(quantComb.selRE[,1]),as.numeric(quantComb.selRE[,2]),
                          as.numeric(quantComb.selRE[,3]),as.numeric(quantComb.selRE[,4]),
                          as.numeric(quantComb.selRE[,5]),as.numeric(quantComb.selRE[,6]))#,
                          #as.numeric(quantComb.selRE[,7]))
quant.maxRE <- rowMaxs(as.matrix(quantComb.selModRE[,c(1:6)]))#7)]))
yMaxRE <- quant.maxRE + (quant.maxRE*0.05)

quant10Comb.selModRE <- cbind(as.numeric(quant10Comb.selRE[,1]),as.numeric(quant10Comb.selRE[,2]),
                            as.numeric(quant10Comb.selRE[,3]),as.numeric(quant10Comb.selRE[,4]),
                            as.numeric(quant10Comb.selRE[,5]),as.numeric(quant10Comb.selRE[,6]))#,
                            #as.numeric(quant10Comb.selRE[,7]))
quant10.minRE <- rowMins(as.matrix(quant10Comb.selModRE[,c(1:6)]))#7)]))
#Since some min values are positive (all sz ests) and some are negative (all the rest), do separately
yMinRE <- c(quant10.minRE[1:29] - (quant10.minRE[1:29]*-0.05), quant10.minRE[30] - (quant10.minRE[30]*0.05),
          quant10.minRE[31:48] - (quant10.minRE[31:48]*-0.05)) 




## Assign colors
colz <- c("#d73027","#fc8d59","#91bfdb","#4575b4","#fdcb44","#fee090","grey60")



## PLOT
par(mfrow=c(4,3), mar=c(1.25,2,1.9,2))  
#bottom, left, top, and right

# Growth model
for (nn in 1:12) {
  plot(c(1:6), medComb.selRE[nn,1:6], col=colz, ylab=NA,
       xaxt = "n", main=namesRE[nn], cex.axis=0.9,cex.main=0.9, pch=19, 
       ylim=c(yMinRE[nn],yMaxRE[nn]), cex=1.2)
  abline(h=0, col="grey80")
  arrows(1:6,as.numeric(medComb.selRE[nn,1:6]),
         1:6, as.numeric(as.matrix(quantComb.selRE[nn,1:6])),
         lwd = 1.25, angle = 90, code = 3, length=0, col=colz)
  arrows(1:6, as.numeric(medComb.selRE[nn,1:6]),
         1:6, as.numeric(as.matrix(quant10Comb.selRE[nn,1:6])), 
         lwd = 1.25, angle = 90, code = 3, length=0, col=colz)
}

# Survival model
for (nn in 13:24) {
  plot(c(1:6), medComb.selRE[nn,1:6], col=colz, ylab=NA, 
       xaxt = "n", main=namesRE[nn], cex.axis=0.9,cex.main=0.9, pch=19, 
       ylim=c(yMinRE[nn],yMaxRE[nn]), cex=1.2)
  abline(h=0, col="grey80")
  arrows(1:6,as.numeric(medComb.selRE[nn,1:6]),
         1:6, as.numeric(as.matrix(quantComb.selRE[nn,1:6])),
         lwd = 1.25, angle = 90, code = 3, length=0, col=colz)
  arrows(1:6, as.numeric(medComb.selRE[nn,1:6]),
         1:6, as.numeric(as.matrix(quant10Comb.selRE[nn,1:6])), 
         lwd = 1.25, angle = 90, code = 3, length=0, col=colz)
}

# p(Repro) model
for (nn in 25:36) {
  plot(c(1:6), medComb.selRE[nn,1:6], col=colz, ylab=NA,
       xaxt = "n", main=namesRE[nn], cex.axis=0.9,cex.main=0.9, pch=19, 
       ylim=c(yMinRE[nn],yMaxRE[nn]), cex=1.2)
  abline(h=0, col="grey80")
  arrows(1:6,as.numeric(medComb.selRE[nn,1:6]),
         1:6, as.numeric(as.matrix(quantComb.selRE[nn,1:6])),
         lwd = 1.25, angle = 90, code = 3, length=0, col=colz)
  arrows(1:6, as.numeric(medComb.selRE[nn,1:6]),
         1:6, as.numeric(as.matrix(quant10Comb.selRE[nn,1:6])), 
         lwd = 1.25, angle = 90, code = 3, length=0, col=colz)
}

# Repro model
par(mfrow=c(5,3))
for (nn in 37:48) {
  plot(c(1:6), medComb.selRE[nn,1:6], col=colz, ylab=NA,
       xaxt = "n", main=namesRE[nn], cex.axis=0.9,cex.main=0.9, pch=19, 
       ylim=c(yMinRE[nn],yMaxRE[nn]), cex=1.2)
  abline(h=0, col="grey80")
  arrows(1:6,as.numeric(medComb.selRE[nn,1:6]),
         1:6, as.numeric(as.matrix(quantComb.selRE[nn,1:6])),
         lwd = 1.25, angle = 90, code = 3, length=0, col=colz)
  arrows(1:6, as.numeric(medComb.selRE[nn,1:6]),
         1:6, as.numeric(as.matrix(quant10Comb.selRE[nn,1:6])), 
         lwd = 1.25, angle = 90, code = 3, length=0, col=colz)
}
plot.new()
legend("center", c("Full dataset","2004-2013","2004-2013 even yrs only","2004-2013 odd yrs only",
                   "2004-2008","2009-2013"), col=colz[1:6], pch=19, cex=1.7,
       horiz=FALSE, bty="y",seg.len=1, xpd="NA")
## -----------------------------------------------------------------




