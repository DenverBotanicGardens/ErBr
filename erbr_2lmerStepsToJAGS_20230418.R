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
#dats <- read.csv("erbr_TagClust_210510.csv", header = TRUE)
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
#jags.mod <- readRDS("erbr_JAGSmod_c3t10s30b10_210602.rds")
#summary(jags.mod)
#plot(jags.mod)
#summ.mod <- summary(jags.mod)
#tail(summ.mod[,1:3], n=34)
#gelman.diag(jags.mod, confidence = 0.95, transform=FALSE)


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





## Make bar graphs comparing median param ests b/w diff datasets
## Load data from previous JAGS runs
#setwd("C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship")
chains <- readRDS("chains.c3t10s30b10_noYRE_20230420.rds")
chains.4to13 <- readRDS("chains.c3t10s30b10_noYRE_4to13_210616.rds")
chains.4to13evn <- readRDS("chains.c3t10s30b10_noYRE_4to13even_210617.rds")
chains.4to13odd <- readRDS("chains.c3t10s30b10_noYRE_4to13odd_210621.rds")
chains.4to8 <- readRDS("chains.c3t10s30b10_noYRE_4to8_210701.rds")
chains.9to13 <- readRDS("chains.c3t10s30b10_noYRE_9to13_210623.rds")

#jags.mod.4to13 <- readRDS("erbr_JAGSmod_c3t10s30b10_noYRE_4to13_210616.rds")
#chains.4to13 <- jags.mod.4to13$mcmc
#chains.4to13 <- bind_rows(lapply(chains.4to13, as.data.frame))

#jags.mod.4to13evn <- readRDS("erbr_JAGSmod_c3t10s30b10_noYRE_4to13even_210617.rds")
#chains.4to13evn <- jags.mod.4to13evn$mcmc
#chains.4to13evn <- bind_rows(lapply(chains.4to13evn, as.data.frame))

#jags.mod.4to13odd <- readRDS("erbr_JAGSmod_c3t10s30b10_noYRE_4to13odd_210621.rds")
#chains.4to13odd <- jags.mod.4to13odd$mcmc
#chains.4to13odd <- bind_rows(lapply(chains.4to13odd, as.data.frame))

#jags.mod.4to8 <- readRDS("erbr_JAGSmod_c3t10s30b10_noYRE_4to8_210701.rds")
#chains.4to8 <- jags.mod.4to8$mcmc
#chains.4to8 <- bind_rows(lapply(chains.4to8, as.data.frame))

#jags.mod.9to13 <- readRDS("erbr_JAGSmod_c3t10s30b10_noYRE_9to13_210623.rds")
#chains.9to13 <- jags.mod.9to13$mcmc
#chains.9to13 <- bind_rows(lapply(chains.9to13, as.data.frame))


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


## CALCULATE VARIANCE OF PARAMETER VALUES
## ** Could change this to SD with colStdevs ** 
#varParams <- as.data.frame(colVars(as.matrix(chains)))
#varParams.4to13 <- as.data.frame(colVars(as.matrix(chains.4to13)))
#varParams.4to13evn <- as.data.frame(colVars(as.matrix(chains.4to13evn)))
#varParams.4to13odd <- as.data.frame(colVars(as.matrix(chains.4to13odd)))
#varParams.4to8 <- as.data.frame(colVars(as.matrix(chains.4to8)))
#varParams.9to13 <- as.data.frame(colVars(as.matrix(chains.9to13)))

#varComb <- as.data.frame(cbind(varParams, varParams.4to13, varParams.4to13evn, varParams.4to13odd, 
#                               varParams.4to8, varParams.9to13))
#colnames(varComb) <- c("Full", "4to13", "4to13evn", "4to13odd", "4to8", "9to13")
#varComb$Names <- colnames(chains)


## Plot
#names.param <- colnames(chains)[55:85]
#names.param <- names.param[c(2:8,12:16,18:23,25:30)] #Remove intercept plots 
#names.paramTitles <- c("Grwth Sz","Grwth Fall Temp","Grwth Summer Temp","Grwth Winter Temp",
#                       "Grwth Fall Precip","Grwth Summer Precip","Grwth Winter Precip","Surv Size",
#                       "Surv Winter Precip","Surv Fall Temp","Surv Summer Temp","Surv Winter Temp",
#                       "p(Repro) Size","p(Repro) Fall Precip","p(Repro) Summer Precip",
#                       "p(Repro) Fall Temp","p(Repro) Summer Temp","p(Repro) Winter Temp",
#                       "Repro Size","Repro Fall Precip","Repro Summer Precip",
#                       "Repro Winter Temp","Repro Fall Temp","Repro Summer Temp")

#colfunc <- colorRampPalette(c("black", "grey90"))
#ylim.vals <- c()
#par(mfrow=c(4,8))  
#par(mar=c(2,3,3,3))

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
## (grwth, surv, p(repro), repro; sz, sumTemp, fallTemp, winTemp, sumPpt, fallPpt, winPpt)
#index<-c(1,8,13,19, #Note the correct way to use index and order, use this to correct below
#         3,11,17,24,
#         2,10,16,23,
#         4,12,18,22,
#         6,15,21,
#         5,14,20,
#         7,9)
index<-c(1,9,5,13,20,17,23,2,24,10,6,14,3,21,18,11,7,15,4,22,19,16,12,8)

names.paramOrd <- names.param[order(index)]
names.paramTitlesOrd <- names.paramTitles[order(index)]


#colfunc <- colorRampPalette(c("black", "grey90"))
#ylim.vals <- c()
#par(mfrow=c(4,8))  
#par(mar=c(2,3,3,3))


#for (nn in 1:length(names.param)) {
#  plotCI(barplot(as.matrix(medComb[which(medComb$Names == names.param[nn]),1:6]), col=colfunc(6),
#                 main=names.param[nn], cex.axis=1.25, beside=T, border=TRUE, space=c(0,0,0,0,0,0)),
#         as.matrix(medComb[which(medComb$Names == names.param[nn]),1:6]), pch=NA, sfrac=0,
#         uiw=as.matrix(varComb[which(varComb$Names == names.param[nn]),1:6]),
#         liw=as.matrix(varComb[which(varComb$Names == names.param[nn]),1:6]), add=TRUE)
#}

#for (nn in 1:length(names.param)) {
#  plotCI(barplot(as.matrix(medComb[which(medComb$Names == names.param[nn]),1:6]), col=colfunc(6),
#                 main=names.paramTitles[nn], cex.axis=1.25, beside=T, border=TRUE, space=c(0,0,0,0,0,0),
#                 cex.main=0.85, xaxt='n'),
#         as.matrix(medComb[which(medComb$Names == names.param[nn]),1:6]), pch=NA, sfrac=0,
#         uiw=(as.matrix(medComb[which(medComb$Names == names.param[nn]),1:6])+as.matrix(varComb[which(varComb$Names == names.param[nn]),1:6])),
#         liw=(as.matrix(medComb[which(medComb$Names == names.param[nn]),1:6])-as.matrix(varComb[which(varComb$Names == names.param[nn]),1:6])), 
#         add=TRUE)
#}

#plotCI(barplot(as.matrix(medComb[which(medComb$Names == names.param[1]),1:6]), col=colfunc(6),
#               main=NA, cex.axis=1.25, beside=T, border=TRUE, space=c(0,0,0,0,0,0),
#               cex.main=0.85, xaxt='n', ylim=c(0,1.5)),
#       as.matrix(medComb[which(medComb$Names == names.param[1]),1:6]), sfrac=0,
#       uiw=1,
#       liw=0.5,add=TRUE)

## Without error bars until can fix above 
#par(mfrow=c(6,6), mar=c(2,3,2,3))  
#bottom, left, top, and right
#for (nn in 1:length(names.param)) {
#  barplot(as.matrix(medComb[which(medComb$Names == names.param[nn]),1:6]), col=colfunc(6), xaxt = "n",
#          main=names.param[nn], cex.axis=1.2, beside=T, border=TRUE, space=c(0,0,0,0,0,0),
#          cex.main=0.95)
#}
#legend("bottomright", colnames(varComb)[1:6], col=colfunc(6),pch=15,cex=1.2,
#       horiz=FALSE, bty="y",seg.len=1, xpd="NA", inset=c(-1.75,0))

#for (nn in 1:length(names.param)) {
#  barplot(as.matrix(medComb[which(medComb$Names == names.param[nn]),1:6]), col=colfunc(6),
#          main=names.paramTitles[nn], cex.axis=1.25, beside=T, border=TRUE, space=c(0,0,0,0,0,0),
#          cex.main=0.85, xaxt='n')
#}

#legend("center", c("Full dataset","2004 to 2013","2004 to 2013 even yrs","2004 to 2014 odd yrs","2004 to 2008",
#                   "2009 to 2013"), col=colfunc(6),pch=15,cex=1.2,
#       horiz=FALSE, bty="y",seg.len=1, xpd="NA", inset=c(-1.75,0))
## ------------------------------------------------------------------------------------------------




## Make modified graph showing median param ests b/w diff datasets as points and 80-90% limits; include glmm estimates

## USE A GLMM FOR EACH VR AND COMPARE PARAM ESTIMATES TO JAGS MODELS ---------------
## Add t+1 climate, sz, & tag into erbr data 
dats <- dats %>% mutate(TagNew1=lead(TagNew), RosNew1=lead(RosNew), Surv1=lead(surv))  
dats <- dats %>% mutate(PptFall1=lead(PptFall), PptWinter1=lead(PptWinter), PptSummer1=lead(PptSummer),
                        TempFall1=lead(TempFall), TempWinter1=lead(TempWinter), TempSummer1=lead(TempSummer))

dats <- dats[which(dats$TagNew == dats$TagNew1),]  #Remove lines with mis-matched individuals 

## Growth
global.grwth <- glmer.nb(RosNew1 ~ scale(RosNew) + scale(PptFall) + scale(PptWinter) + scale(PptSummer) + 
                         scale(TempFall) + scale(TempWinter) + scale(TempSummer) + (1|TransectNew), data=dats)
global.grwthNS <- glmer.nb(RosNew1 ~ RosNew + PptFall + PptWinter + PptSummer + 
                           TempFall + TempWinter + TempSummer + (1|TransectNew), data=dats)

## Survival #** Should size be logged in surv and prob repro models? **
global.surv <- glmer(Surv1 ~ scale(RosNew) + scale(PptWinter) + scale(TempFall) + scale(TempWinter) + 
                     scale(TempSummer) + (1|TransectNew), family=binomial(link='logit'), data=dats)
global.survNS <- glmer(Surv1 ~ RosNew + PptWinter + TempFall + (TempWinter) + 
                       (TempSummer) + (1|TransectNew), family=binomial(link='logit'), data=dats)

## Probability of reproduction
global.reproYesNo <- glmer(InflYesNo ~ scale(RosNew) + scale(PptFall) + scale(PptSummer) + scale(TempFall) + 
                           scale(TempWinter) + scale(TempSummer) + (1|TransectNew), family=binomial(link=logit), data=dats)
global.reproYesNoNS <- glmer(InflYesNo ~ (RosNew) + (PptFall) + (PptSummer) + (TempFall) + 
                             (TempWinter) + (TempSummer) + (1|TransectNew), family=binomial(link=logit), data=dats)

## Reproduction #** Should data be subset to only include reproductive plts (i.e. infs>0) in this model? **
global.repro <- glmer.nb(InflNew ~ scale(RosNew) + scale(PptFall) + scale(PptSummer) + scale(TempFall) + 
                           scale(TempWinter) + scale(TempSummer) + (1|TransectNew), data=dats)
global.reproNS <- glmer.nb(InflNew ~ (RosNew) + (PptFall) + (PptSummer) + (TempFall) + 
                           (TempWinter) + (TempSummer) + (1|TransectNew), data=dats)


## Extract parameter estimates and SEs from GLMMs
#How to get in comparable units/ scale to JAGS output? 
paramsMM.grwth <- as.data.frame(cbind(as.data.frame(summary(global.grwth)$coefficients[2:nrow(summary(global.grwth)$coefficients),1:2]),
                        c("Grwth Size","Grwth Fall Precip","Grwth Winter Precip","Grwth Summer Precip",
                          "Grwth Fall Temp","Grwth Winter Temp","Grwth Summer Temp")))
colnames(paramsMM.grwth) <- c("GLMM","SE","ParamTitle")
paramsMM.surv <- as.data.frame(cbind(as.data.frame(summary(global.surv)$coefficients[2:nrow(summary(global.surv)$coefficients),1:2]),
                                      c("Surv Size","Surv Winter Precip",
                                        "Surv Fall Temp","Surv Winter Temp","Surv Summer Temp")))
colnames(paramsMM.surv) <- c("GLMM","SE","ParamTitle")
paramsMM.reproYesNo <- as.data.frame(cbind(as.data.frame(summary(global.reproYesNo)$coefficients[2:nrow(summary(global.reproYesNo)$coefficients),1:2]),
                                     c("p(Repro) Size","p(Repro) Fall Precip","p(Repro) Summer Precip",
                                       "p(Repro) Fall Temp","p(Repro) Winter Temp","p(Repro) Summer Temp")))
colnames(paramsMM.reproYesNo) <- c("GLMM","SE","ParamTitle")
paramsMM.repro <- as.data.frame(cbind(as.data.frame(summary(global.repro)$coefficients[2:nrow(summary(global.repro)$coefficients),1:2]),
                                           c("Repro Size","Repro Fall Precip","Repro Summer Precip",
                                             "Repro Fall Temp","Repro Winter Temp","Repro Summer Temp")))
colnames(paramsMM.repro) <- c("GLMM","SE","ParamTitle")
paramsMM <- rbind(paramsMM.grwth, paramsMM.surv, paramsMM.reproYesNo, paramsMM.repro)
paramsMM$SE_upr <- paramsMM$GLMM + paramsMM$SE
paramsMM$SE_lwr <- paramsMM$GLMM - paramsMM$SE


## Adjust param estimates to account for scaling? ***
sd(dats$RosNew, na.rm=TRUE)
paramsMM[1,1]/sd(dats$RosNew, na.rm=TRUE)
sd(dats$PptSummer, na.rm=TRUE)
paramsMM[15,1]/sd(dats$PptSummer, na.rm=TRUE)
## --------------------------------------------------------------




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



## Make dataframe of selected relevant values for plotting 
medComb.sel <- NULL
for (nn in 1:length(names.paramOrd)) {
  medComb.sel <- rbind(medComb.sel, as.data.frame(medComb[which(medComb$Names == names.paramOrd[nn]),1:6]))
}

## Combine JAGS and GLMM estimates
medComb.sel$ParamTitle <- names.paramTitlesOrd
medComb.sel <- dplyr::left_join(medComb.sel, paramsMM, by="ParamTitle")
medComb.sel <- medComb.sel %>% relocate(ParamTitle, .after=last_col())
#medComb.sel <- rbind(medComb.sel[1:17,],NA,medComb.sel[18:24,])


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

#quant10Comb.trunc <- as.data.frame(quant10Comb[55:85,])
#quant10Comb.trunc <- quant10Comb.trunc[c(2:8,12:16,18:23,25:30),1:6] #Remove intercept plots
#quantComb.trunc <- as.data.frame(quantComb[55:85,])
#quantComb.trunc <- quantComb.trunc[c(2:8,12:16,18:23,25:30),1:6] #Remove intercept plots 




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
#colfunc <- colorRampPalette(c("red", "blue"))
colz <- c("#d73027","#fc8d59","#91bfdb","#4575b4","#fdcb44","#fee090","grey60")


 


## PLOT
#tiff
pdf('ErBr_fig2_20231022.pdf', width=6, height=9)
par(mfrow=c(7,4), mar=c(2,2,2,2))  #Plot so 4 VR models are in cols and upto 7 predictor vars are in rows
#bottom, left, top, and right
for (nn in 1:17) {
  plot(c(1:7), medComb.sel[nn,1:7], col=colz, ylab=NA,
       xaxt = "n", main=names.paramTitlesOrd[nn], cex.axis=1.1,cex.main=0.9, pch=19, 
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
       xaxt = "n", main=names.paramTitlesOrd[nn], cex.axis=1.1,cex.main=0.9, pch=19, 
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
       xaxt = "n", main=names.paramTitlesOrd[nn], cex.axis=1.1,cex.main=0.9, pch=19, 
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

legend("center", colnames(medComb.sel)[1:7], col=colz, pch=19, cex=1,
       horiz=FALSE, bty="y",seg.len=1, xpd="NA")

dev.off()
## -----------------------------------------------------------------
#plot.new()




#legend("bottomright", colnames(varComb)[1:6], col=colfunc(6),pch=15,cex=1.2,
#       horiz=FALSE, bty="y",seg.len=1, xpd="NA", inset=c(-1.75,0))

#for (nn in 1:length(names.param)) {
#  barplot(as.matrix(medComb[which(medComb$Names == names.param[nn]),1:6]), col=colfunc(6),
#          main=names.paramTitles[nn], cex.axis=1.25, beside=T, border=TRUE, space=c(0,0,0,0,0,0),
#          cex.main=0.85, xaxt='n')
#}

#legend("center", c("Full dataset","2004 to 2013","2004 to 2013 even yrs","2004 to 2014 odd yrs","2004 to 2008",
#                   "2009 to 2013"), col=colfunc(6),pch=15,cex=1.2,
#       horiz=FALSE, bty="y",seg.len=1, xpd="NA", inset=c(-1.75,0))

## ------------------------------------------------------------------------------------------------





## Make modified graph showing median param ests b/w diff datasets as points and 80-90% limits
## CALCULATE 90th PERCENTILE 
#quantParams <- chains %>% summarise_all(funs(list(quantile(., probs=0.9)))) #%>% transpose
#quantParams <- as.data.frame(t(quantParams))
#quantParams.4to13 <- chains.4to13 %>% summarise_all(funs(list(quantile(., probs=0.9))))
#quantParams.4to13 <- as.data.frame(t(quantParams.4to13))
#quantParams.4to13evn <- chains.4to13evn %>% summarise_all(funs(list(quantile(., probs=0.9))))
#quantParams.4to13evn <- as.data.frame(t(quantParams.4to13evn))
#quantParams.4to13odd <- chains.4to13odd %>% summarise_all(funs(list(quantile(., probs=0.9))))
#quantParams.4to13odd <- as.data.frame(t(quantParams.4to13odd))
#quantParams.4to8 <- chains.4to8 %>% summarise_all(funs(list(quantile(., probs=0.9))))
#quantParams.4to8 <- as.data.frame(t(quantParams.4to8))
#quantParams.9to13 <- chains.9to13 %>% summarise_all(funs(list(quantile(., probs=0.9))))
#quantParams.9to13 <- as.data.frame(t(quantParams.9to13))

#quantComb <- as.data.frame(cbind(quantParams, quantParams.4to13, quantParams.4to13evn, quantParams.4to13odd, 
#                               quantParams.4to8, quantParams.9to13))
#colnames(quantComb) <- c("Full", "4to13", "4to13evn", "4to13odd", "4to8", "9to13")
#quantComb$Names <- colnames(chains)


## CALCULATE 10th PERCENTILE 
#quant10Params <- chains %>% summarise_all(funs(list(quantile(., probs=0.1)))) 
#quant10Params <- as.data.frame(t(quant10Params))
#quant10Params.4to13 <- chains.4to13 %>% summarise_all(funs(list(quantile(., probs=0.1))))
#quant10Params.4to13 <- as.data.frame(t(quant10Params.4to13))
#quant10Params.4to13evn <- chains.4to13evn %>% summarise_all(funs(list(quantile(., probs=0.1))))
#quant10Params.4to13evn <- as.data.frame(t(quant10Params.4to13evn))
#quant10Params.4to13odd <- chains.4to13odd %>% summarise_all(funs(list(quantile(., probs=0.1))))
#quant10Params.4to13odd <- as.data.frame(t(quant10Params.4to13odd))
#quant10Params.4to8 <- chains.4to8 %>% summarise_all(funs(list(quantile(., probs=0.1))))
#quant10Params.4to8 <- as.data.frame(t(quant10Params.4to8))
#quant10Params.9to13 <- chains.9to13 %>% summarise_all(funs(list(quantile(., probs=0.1))))
#quant10Params.9to13 <- as.data.frame(t(quant10Params.9to13))

#quant10Comb <- as.data.frame(cbind(quant10Params, quant10Params.4to13, quant10Params.4to13evn, quant10Params.4to13odd, 
#                                 quant10Params.4to8, quant10Params.9to13))
#colnames(quant10Comb) <- c("Full", "4to13", "4to13evn", "4to13odd", "4to8", "9to13")
#quant10Comb$Names <- colnames(chains)



## PLOT
#quant10Comb.trunc <- as.data.frame(quant10Comb[55:85,])
#quant10Comb.trunc <- quant10Comb.trunc[c(2:8,12:16,18:23,25:30),1:6] #Remove intercept plots
#quant10Comb.trunc <- cbind(as.numeric(quant10Comb.trunc[,1]),as.numeric(quant10Comb.trunc[,2]),
#                           as.numeric(quant10Comb.trunc[,3]),as.numeric(quant10Comb.trunc[,4]),
#                           as.numeric(quant10Comb.trunc[,5]),as.numeric(quant10Comb.trunc[,6])) #Change structure

#quantComb.trunc <- as.data.frame(quantComb[55:85,])
#quantComb.trunc <- quantComb.trunc[c(2:8,12:16,18:23,25:30),1:6] #Remove intercept plots 
#quantComb.trunc <- cbind(as.numeric(quantComb.trunc[,1]),as.numeric(quantComb.trunc[,2]),
#                           as.numeric(quantComb.trunc[,3]),as.numeric(quantComb.trunc[,4]),
#                           as.numeric(quantComb.trunc[,5]),as.numeric(quantComb.trunc[,6]))


#quant10.min <- rowMins(as.matrix(quant10Comb.trunc[,c(1:6)]))
#quant.max <- rowMaxs(as.matrix(quantComb.trunc[,c(1:6)]))


#yMin <- quant10.min - (quant10.min*0.1)
  #c(0.7,-0.22,-0.5,-2,-0.05,-0.03,-0.03,0.2,-0.1,-2,-4,-3,-0.08,-0.02,-3,-1,-1.5,0.7,-0.02,-0.005,0.1,-0.95,-0.4)
#yMax <- quant.max + (quant.max*0.04)
  #as.numeric(quantParams.trunc) + (as.numeric(quantParams.trunc)*0.1)
  #c(0.85,0.3,0.5,1,-0.001,0.01,0.04,1.4,0.08,1,1.5,2.5,1.7,0.08,0.03,2,1,1.5,1.2,0.009,0.009,0.6,-0.2,0.5)


#colfunc <- colorRampPalette(c("red", "blue"))

#tiff
#pdf('ErBr_fig2_20231022.pdf', width=6, height=9)
#par(mfrow=c(6,4), mar=c(2,3,2,2))  
#bottom, left, top, and right
#for (nn in 1:length(names.param)) {
#  plot(c(1:6), as.matrix(medComb[which(medComb$Names == names.param[nn]),1:6]), col=colfunc(6), 
#       xaxt = "n", main=names.paramTitles[nn], cex.axis=1.1,cex.main=0.9, pch=19, 
#       ylim=c(yMin[nn],yMax[nn]), cex=1.2)
#  arrows(1:6,as.matrix(medComb[which(medComb$Names == names.param[nn]),1:6]),
#         1:6,as.numeric(as.matrix(quantComb[which(quantComb$Names == names.param[nn]),1:6])), 
#         lwd = 1.25, angle = 90, code = 3, length=0, col=colfunc(6))
#  arrows(1:6,as.matrix(medComb[which(medComb$Names == names.param[nn]),1:6]),
#         1:6,as.numeric(as.matrix(quant10Comb[which(quant10Comb$Names == names.param[nn]),1:6])), 
#         lwd = 1.25, angle = 90, code = 3, length=0, col=colfunc(6))
#}

#legend("bottomright", colnames(varComb)[1:6], col=colfunc(6),pch=19,cex=1.2,
#       horiz=FALSE, bty="y",seg.len=1, xpd=NA, inset=c(-1.2,0))

#dev.off()


## -----------------------------------------------------------------

#main=names.paramTitles[nn], cex.axis=1.25, beside=T, border=TRUE, space=c(0,0,0,0,0,0),
#                 cex.main=0.85, xaxt='n'







## Make bar graph comparing median param ests & 80%(?) CIs b/w diff datasets
## Load model from each dataset
#modFull <- readRDS("erbr_JAGSmod_c3t10s30b10_noYRE_210613.rds")
#mod4to13 <- readRDS("erbr_JAGSmod_c3t10s30b10_noYRE_4to13_210616.rds")
#mod4to13evn <- readRDS("erbr_JAGSmod_c3t10s30b10_noYRE_4to13even_210617.rds")
#mod4to13odd <- readRDS("erbr_JAGSmod_c3t10s30b10_noYRE_4to13odd_210621.rds")
#mod4to8 <- readRDS("erbr_JAGSmod_c3t10s30b10_noYRE_4to8_210701.rds")
#mod9to13 <- readRDS("erbr_JAGSmod_c3t10s30b10_noYRE_9to13_210623.rds")

## Obtain summaries
#modFull.summTot <- summary(modFull)
#mod4to13.summTot <- summary(mod4to13)
#mod4to13evn.summTot <- summary(mod4to13evn)
#mod4to13odd.summTot <- summary(mod4to13odd)
#mod4to8.summTot <- summary(mod4to8)
#mod9to13.summTot <- summary(mod9to13)

#modFull.summ <- as.data.frame(tail(modFull.summTot[,1:3], n=31))
#mod4to13.summ <- as.data.frame(tail(mod4to13.summTot[,1:3], n=31))
#mod4to13evn.summ <- as.data.frame(tail(mod4to13evn.summTot[,1:3], n=31))
#mod4to13odd.summ <- as.data.frame(tail(mod4to13odd.summTot[,1:3], n=31))
#mod4to8.summ <- as.data.frame(tail(mod4to8.summTot[,1:3], n=31))
#mod9to13.summ <- as.data.frame(tail(mod9to13.summTot[,1:3], n=31))

## Add column with data set name
#modFull.summ$DataSet <- "Full"  
#mod4to13.summ$DataSet <- "4to13"
#mod4to13evn.summ$DataSet <- "4to13evn"
#mod4to13odd.summ$DataSet <- "4to13odd"
#mod4to8.summ$DataSet <- "4to8"
#mod9to13.summ$DataSet <- "9to13"

#modFull.summ$Param <- rownames(modFull.summ)
#mod4to13.summ$Param <- rownames(mod4to13.summ)
#mod4to13evn.summ$Param <- rownames(mod4to13evn.summ)
#mod4to13odd.summ$Param <- rownames(mod4to13odd.summ)
#mod4to8.summ$Param <- rownames(mod4to8.summ)
#mod9to13.summ$Param <- rownames(mod9to13.summ)

## Combine summaries from different data sets
#mod.comb <- rbind(modFull.summ, mod4to13.summ, mod4to13evn.summ, mod4to13odd.summ, 
#                  mod4to8.summ, mod9to13.summ)
#mod.comb.sort <- arrange(mod.comb, Param)


## Subset by vital rate and plot
#colfunc <- colorRampPalette(c("black", "white"))

#mod.comb.grwth <- dplyr::filter(mod.comb.sort, grepl("grwth_", Param))
#plotCI(barplot(mod.comb.grwth$Median, beside=T, ylab="Parameter estimate", xlab=NA, border=FALSE, 
#       ylim=c(-0.01,2), cex.axis=1, cex.lab=1, col=colfunc(6)), 
#       mod.comb.grwth$Median, uiw=mod.comb.grwth$Upper95, liw=mod.comb.grwth$Lower95, 
#       add=TRUE, pch=NA, sfrac=0)
#space=c(rep(0,2),0.3,rep(0,2),0.3,rep(0,2),0.3,rep(0,2))
#col=rep(c(col.crs,col.P1,col.P2,col.P3,col.P4),3),
#space=c(0,rep(c(0,0.5),7)),

#mod.comb.surv <- dplyr::filter(mod.comb.sort, grepl("surv_", Param))
#plotCI(barplot(mod.comb.surv$Median, beside=T, ylab="Parameter estimate", xlab=NA, border=FALSE, 
#       ylim=c(-3,2), cex.axis=1, cex.lab=1, col=colfunc(6)), 
#       mod.comb.surv$Median, uiw=mod.comb.surv$Upper95, liw=mod.comb.surv$Lower95, 
#       add=TRUE, pch=NA, sfrac=0)

#mod.comb.reproYesNo <- dplyr::filter(mod.comb.sort, grepl("reproyesno_", Param))
#plotCI(barplot(mod.comb.reproYesNo$Median, beside=T, ylab="Parameter estimate", xlab=NA, border=FALSE, 
#       ylim=c(-3,2), cex.axis=1, cex.lab=1, col=colfunc(6)), 
#       mod.comb.reproYesNo$Median, uiw=mod.comb.reproYesNo$Upper95, liw=mod.comb.reproYesNo$Lower95, 
#       add=TRUE, pch=NA, sfrac=0)

#mod.comb.repro <- dplyr::filter(mod.comb.sort, grepl("repro_", Param))
#plotCI(barplot(mod.comb.repro$Median, beside=T, ylab="Parameter estimate", xlab=NA, border=FALSE, 
#       ylim=c(-1,2), cex.axis=1, cex.lab=1, col=colfunc(6)), 
#       mod.comb.repro$Median, uiw=mod.comb.repro$Upper95, liw=mod.comb.repro$Lower95, 
#       add=TRUE, pch=NA, sfrac=0)
## ------------------------------------------------------------------------------------------------





 
