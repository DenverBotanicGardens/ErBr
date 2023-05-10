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
## ------------------------------------------------------------------------------------------------


setwd("~/ErBr")


## LOAD DATA --------------------------------------------------------------------------------------
dats <- read.csv("erbr_TagClust2022_20230408.csv", header = TRUE)
#dats <- read.csv("erbr_TagClust_210510.csv", header = TRUE)
#dats <- read.csv("erbr_TagClust4to8_210617.csv", header = TRUE)
## ------------------------------------------------------------------------------------------------



## SET WD (WHERE JAGS SCRIPT IS LOCATED) ----------------------------------------------------------
#setwd("C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship/R_scripts")
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
jags.mod <- run.jags('Scripts/erbr_JAGSmodBest_noYRE_20230418.R', n.chains=3, data=dats, burnin=10000, thin=10, sample=30000, adapt=500, method='parallel')

#save(jags.mod, file='erbr_JAGSmod_c3t10s20b5_210406.rdata')
saveRDS(jags.mod, "erbr_JAGSmodBest_c3t10s30b10_noYRE_20230420.rds")
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


