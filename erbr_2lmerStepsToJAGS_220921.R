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
library(plotrix)
library(robustbase)
library(resample)
## ------------------------------------------------------------------------------------------------



## LOAD DATA --------------------------------------------------------------------------------------
dats <- read.csv("C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship/erbr_TagClust_210510.csv", header = TRUE)
#dats <- read.csv("erbr_tagClust4to8_210504.csv", header = TRUE)
## ------------------------------------------------------------------------------------------------



## SET WD (WHERE JAGS SCRIPT IS LOCATED) ----------------------------------------------------------
setwd("C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship/R_scripts")
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
dats.newPlts <- dats.newPlts %>% mutate(PptFall1=lead(PptFall), PptWinter1=lead(PptWinter), PptSummer1=lead(PptSummer),
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




## RUN ASSOCIATED JAGS MODEL ----------------------------------------------------------------------
jags.mod <- run.jags('erbr_3JAGSmodBest_noYRE_210608.R', n.chains=3, data=dats, burnin=10000, thin=10, sample=30000, adapt=500, method='parallel')

#save(jags.mod, file='erbr_JAGSmod_c3t10s20b5_210406.rdata')
#saveRDS(jags.mod, "erbr_JAGSmod_c3t10s30b10_noYRE_210613.rds")
## ------------------------------------------------------------------------------------------------




## LOOK AT MODEL OUTPUT ---------------------------------------------------------------------------
jags.mod <- readRDS("erbr_JAGSmod_c3t10s30b10_noYRE_210613.rds")
summary(jags.mod)  #Should show median, mean, etc. for each parameter from JAGS model
plot(jags.mod)     #Should plot trace plots, CDFs, etc. 
summ.mod <- summary(jags.mod)
tail(summ.mod[,1:5], n=32)
gelman.diag(jags.mod, confidence = 0.95, transform=FALSE)

chains <- jags.mod$mcmc
chains <- bind_rows(lapply(chains, as.data.frame))
colMeds <- apply(chains,2,median)
colSDs <- apply(chains,2,sd)


## Look at correlation b/w params 
chains.1 <- chains %>% dplyr::select(!contains(c("randomeffect", "precision")))
chains.1 <- chains.1 %>% dplyr::select(!c(deviance, resid.sum.sq))

cor.chains <- cor(chains.1)
corrplot(cor.chains, method="circle", type="lower")




## Make bar graphs comparing median param ests b/w diff datasets
## Load data from previous JAGS runs
chains <- readRDS("chains.c3t10s30b10_noYRE_210613.rds")

jags.mod.4to13 <- readRDS("erbr_JAGSmod_c3t10s30b10_noYRE_4to13_210616.rds")
chains.4to13 <- jags.mod.4to13$mcmc
chains.4to13 <- bind_rows(lapply(chains.4to13, as.data.frame))

jags.mod.4to13evn <- readRDS("erbr_JAGSmod_c3t10s30b10_noYRE_4to13even_210617.rds")
chains.4to13evn <- jags.mod.4to13evn$mcmc
chains.4to13evn <- bind_rows(lapply(chains.4to13evn, as.data.frame))

jags.mod.4to13odd <- readRDS("erbr_JAGSmod_c3t10s30b10_noYRE_4to13odd_210621.rds")
chains.4to13odd <- jags.mod.4to13odd$mcmc
chains.4to13odd <- bind_rows(lapply(chains.4to13odd, as.data.frame))

jags.mod.4to8 <- readRDS("erbr_JAGSmod_c3t10s30b10_noYRE_4to8_210701.rds")
chains.4to8 <- jags.mod.4to8$mcmc
chains.4to8 <- bind_rows(lapply(chains.4to8, as.data.frame))

jags.mod.9to13 <- readRDS("erbr_JAGSmod_c3t10s30b10_noYRE_9to13_210623.rds")
chains.9to13 <- jags.mod.9to13$mcmc
chains.9to13 <- bind_rows(lapply(chains.9to13, as.data.frame))


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
varParams <- as.data.frame(colVars(as.matrix(chains)))
varParams.4to13 <- as.data.frame(colVars(as.matrix(chains.4to13)))
varParams.4to13evn <- as.data.frame(colVars(as.matrix(chains.4to13evn)))
varParams.4to13odd <- as.data.frame(colVars(as.matrix(chains.4to13odd)))
varParams.4to8 <- as.data.frame(colVars(as.matrix(chains.4to8)))
varParams.9to13 <- as.data.frame(colVars(as.matrix(chains.9to13)))

varComb <- as.data.frame(cbind(varParams, varParams.4to13, varParams.4to13evn, varParams.4to13odd, 
                               varParams.4to8, varParams.9to13))
colnames(varComb) <- c("Full", "4to13", "4to13evn", "4to13odd", "4to8", "9to13")
varComb$Names <- colnames(chains)


## Plot
names.param <- colnames(chains)[55:85]
colfunc <- colorRampPalette(c("black", "white"))
ylim.vals <- c()
par(mfrow=c(4,8))  

for (nn in 1:length(names.param)) {
  plotCI(barplot(as.matrix(medComb[which(medComb$Names == names.param[nn]),1:6]), col=colfunc(6),
          main=names.param[nn], cex.axis=1.25, beside=T, border=TRUE, space=c(0,0,0,0,0,0)),
         as.matrix(medComb[which(medComb$Names == names.param[nn]),1:6]), pch=NA, sfrac=0,
         uiw=as.matrix(varComb[which(varComb$Names == names.param[nn]),1:6]),
         liw=as.matrix(varComb[which(varComb$Names == names.param[nn]),1:6]), add=TRUE)
}
## ------------------------------------------------------------------------------------------------


