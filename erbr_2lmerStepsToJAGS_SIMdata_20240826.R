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
library(stringr)
## ------------------------------------------------------------------------------------------------







## LOOP OVER DATASETS -----------------------------------------------------------------------------
n.datset <- 10
for (dd in 1:n.datset) {

  old <- Sys.time() # get start time

## LOAD DATA --------------------------------------------------------------------------------------
#dats <- read.csv("20240904_erbr_SimDat20yrR2_Format4JAGS.csv", header = TRUE)
name <- as.character("SimDat20yrMiss.")
dats <- read.csv(file=paste("20240906", "_erbr_", name, dd, ".Format4JAGS", ".csv", sep=""), header = TRUE)
## ------------------------------------------------------------------------------------------------



## SET WD (WHERE JAGS SCRIPT IS LOCATED) ----------------------------------------------------------
#setwd("C:/Users/april/Dropbox/CU_Boulder_PhD/DBG_Internship")
#setwd("~/ErBr")
## ------------------------------------------------------------------------------------------------



## CHANGE NAMES OF CLIMATE COLUMNS ------------------------------------------------------------------
dats <- rename(dats, PptFall=Tot_fall_ppt, PptWinter=Tot_winter_ppt, PptSummer=Tot_summer_ppt,
               TempFall=Mean_fall_temp, TempWinter=Mean_winter_temp, TempSummer=Mean_summer_temp)
## ------------------------------------------------------------------------------------------------




## COUNT PLANT OBERVATION YEARS FOR EACH DATASET --------------------------------------------------
sum(dats$surv, na.rm=TRUE) #Full dataset
#consec <- dats[dats$Year<2014,]
#sum(consec$surv, na.rm=TRUE) #Consecutive dataset
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






## RUN ASSOCIATED JAGS MODEL ----------------------------------------------------------------------
#jags.mod <- run.jags('erbr_JAGSmodComplx_noYRE_210827.R', n.chains=3, data=dats, burnin=10000, thin=10, sample=30000, adapt=500, method='parallel')

jags.mod <- run.jags('erbr_3JAGSmodBest_noYRE_20230418short.R', n.chains=3, data=dats, burnin=10000, thin=5, sample=10000, adapt=500, method='parallel')
#jags.mod <- run.jags('erbr_3JAGSmodBest_noYRE_20230418.R', n.chains=3, data=dats, burnin=10000, thin=10, sample=30000, adapt=500, method='parallel')



## Save output
date <- Sys.Date()                                #Enter date to be added to file name
date <- str_replace_all(date, "-", "")
saveRDS(jags.mod, file=paste(date, "_erbr_JAGSmodBest_c3t5s10b10_noYRE_SimDat20yr.", dd, ".rds", sep=""))

summ.mod <- summary(jags.mod)
saveRDS(summ.mod, file=paste(date, "_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat20yr.", dd, ".rds", sep=""))


}   #End dataset loop



time.diff <- Sys.time() - old #calculate difference
print(time.diff) #print in nice format
## ------------------------------------------------------------------------------------------------






## LOOK AT MODEL OUTPUT ---------------------------------------------------------------------------
#jags.mod <- readRDS("erbr_JAGSmodBest_SIM20yr_c3t5s10b5_noYRE_20240830.rds")
#jags.mod <- readRDS("erbr_JAGSmodBest_SIM50yr_c3t5s10b5_noYRE_20240831.rds")
jags.mod.sim10 <- readRDS("erbr_JAGSmodBest_SIM10yr_c3t5s10b5_noYRE_20240901.rds")

summary(jags.mod)
plot(jags.mod)
summ.mod <- summary(jags.mod)
tail(summ.mod[,1:3], n=37)
gelman.diag(jags.mod, confidence = 0.95, transform=FALSE)


## Make chains variable
chains <- jags.mod$mcmc
chains <- bind_rows(lapply(chains, as.data.frame))
chains.sim10 <- jags.mod.sim10$mcmc
chains.sim10 <- bind_rows(lapply(chains.sim10, as.data.frame))
#colMeds <- apply(chains,2,median)
#colSDs <- apply(chains,2,sd)
## -----------------------------------------------------------------------------------------------



## Make plot comparing median param ests of diff simulated datasets & compare to param values used on making sim data
## Show median param ests of diff datasets as points and 90% limits

## Load data from previous JAGS runs
saveRDS(chains, "chains_SIM20yr_c3t5s10b5_noYRE_20240901.rds")
saveRDS(chains, "chains_SIM50yr_c3t5s10b5_noYRE_20240902.rds")
saveRDS(chains, "chains_SIM10yr_c3t5s10b5_noYRE_20240902.rds")

chains.sim20 <- readRDS("chains_SIM20yr_c3t5s10b5_noYRE_20240901.rds")
chains.sim50 <- readRDS("chains_SIM50yr_c3t5s10b5_noYRE_20240902.rds")
chains.sim10 <- readRDS("chains_SIM10yr_c3t5s10b5_noYRE_20240902.rds")


# Load data from previous JAGS runs of real data
medParams.realDat <- readRDS("erbrMedParams_noYRE_20240803")
medParams.realDatTr <- as.data.frame(t(medParams.realDat))
medParams.realDatTr <- as.data.frame(cbind(medParams.realDatTr$`colMedians(as.matrix(chains))`,colnames(medParams.realDat)))
colnames(medParams.realDatTr) <- c("realData", "Name")



## CALCULATE MEDIAN PARAMETER VALUES 
medParams.sim20 <- as.data.frame(colMedians(as.matrix(chains.sim20)))
#medParams.sim20wName <- as.data.frame(cbind(medParams.sim20$`colMedians(as.matrix(chains.sim20))`,colnames(chains.sim20)))
#colnames(medParams.sim20wName) <- c("Sim20", "Name")

medParams.sim50 <- as.data.frame(colMedians(as.matrix(chains.sim50)))
#medParams.sim50wName <- as.data.frame(cbind(medParams.sim50$`colMedians(as.matrix(chains.sim50))`,colnames(chains.sim50)))
#colnames(medParams.sim50wName) <- c("sim50", "Name")

medParams.sim10 <- as.data.frame(colMedians(as.matrix(chains.sim10)))
#medParams.sim10wName <- as.data.frame(cbind(medParams.sim10$`colMedians(as.matrix(chains.sim10))`,colnames(chains.sim10)))
#colnames(medParams.sim10wName) <- c("sim10", "Name")

medComb <- as.data.frame(cbind(medParams.sim10$`colMedians(as.matrix(chains.sim10))`, 
                               medParams.sim20$`colMedians(as.matrix(chains.sim20))`, 
                               medParams.sim50$`colMedians(as.matrix(chains.sim50))`,
                               colnames(chains.sim20)))
colnames(medComb) <- c("Sim10", "Sim20", "Sim50","Name")
#medComb$Names <- colnames(chains.sim20)




## CALCULATE 90th PERCENTILE ON JAGS ESTIMATES ------------------
quantParams.sim20 <- chains.sim20 %>% summarise_all(funs(quantile(., probs=0.9)))
quantParams.sim20 <- as.data.frame(t(quantParams.sim20))
#quantParams.sim20wName <- as.data.frame(cbind(as.numeric(quantParams.sim20$V1),colnames(chains.sim20)))
#colnames(quantParams.sim20wName) <- c("Sim20", "Name")

quantParams.sim50 <- chains.sim50 %>% summarise_all(funs(quantile(., probs=0.9)))
quantParams.sim50 <- as.data.frame(t(quantParams.sim50))

quantParams.sim10 <- chains.sim10 %>% summarise_all(funs(quantile(., probs=0.9)))
quantParams.sim10 <- as.data.frame(t(quantParams.sim10))

quantComb <- as.data.frame(cbind(as.numeric(quantParams.sim10$V1), as.numeric(quantParams.sim20$V1), 
                           as.numeric(quantParams.sim50$V1),colnames(chains.sim20)))
colnames(quantComb) <- c("Sim10", "Sim20", "Sim50","Name")
#quantComb$Names <- colnames(chains)


## CALCULATE 10th PERCENTILE ON JAGS ESTIMATES
quant10Params.sim20 <- chains.sim20 %>% summarise_all(funs(quantile(., probs=0.1)))
quant10Params.sim20 <- as.data.frame(t(quant10Params.sim20))
#quant10Params.sim20wName <- as.data.frame(cbind(as.numeric(quant10Params.sim20$V1),colnames(chains.sim20)))
#colnames(quant10Params.sim20wName) <- c("Sim20", "Name")

quant10Params.sim50 <- chains.sim50 %>% summarise_all(funs(list(quantile(., probs=0.1))))
quant10Params.sim50 <- as.data.frame(t(quant10Params.sim50))

quant10Params.sim10 <- chains.sim10 %>% summarise_all(funs(list(quantile(., probs=0.1))))
quant10Params.sim10 <- as.data.frame(t(quant10Params.sim10))

quant10Comb <- as.data.frame(cbind(as.numeric(quant10Params.sim10$V1), as.numeric(quant10Params.sim20$V1), 
                                   as.numeric(quant10Params.sim50$V1),colnames(chains.sim20)))
colnames(quant10Comb) <- c("Sim10", "Sim20", "Sim50","Name")
#quant10Comb$Names <- colnames(chains)
## ----------------------------------------------------------------------------




## SUBSET OUTPUT TO JUST KEEP PARAMS OF INTEREST
names.param <- colnames(chains.sim20)[50:80]
names.param <- names.param[c(2:8,12:16,18:23,25:30)] #Remove intercept and GrwthVar  
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




## Make dataframe of selected relevant values for plotting 
medComb.sel <- NULL
for (nn in 1:length(names.paramOrd)) {
  medComb.sel <- rbind(medComb.sel, as.data.frame(medComb[which(medComb$Name==names.paramOrd[nn]),1:3]))
}

quantComb.sel <- NULL
for (uu in 1:length(names.paramOrd)) {
  quantComb.sel <- rbind(quantComb.sel, as.data.frame(quantComb[which(quantComb$Name==names.paramOrd[uu]),1:3]))
}

quant10Comb.sel <- NULL
for (uu in 1:length(names.paramOrd)) {
  quant10Comb.sel <- rbind(quant10Comb.sel, as.data.frame(quant10Comb[which(quant10Comb$Name==names.paramOrd[uu]),1:3]))
}

medRealDat.sel <- NULL
for (nn in 1:length(names.paramOrd)) {
  medRealDat.sel <- rbind(medRealDat.sel, as.data.frame(medParams.realDatTr[which(medParams.realDatTr$Name==names.paramOrd[nn]),1]))
}






## Calculate the min and max values of error bars & set appropriate y-axis values for plotting 
# Change structure to allow for row min and max calcs 
#quantComb.sel <- as.numeric(as.character(quantComb.sel[,1]))
quantComb.selMod <- cbind(as.numeric(as.character(quantComb.sel$Sim10)),as.numeric(as.character(quantComb.sel$Sim20)),
                          as.numeric(as.character(quantComb.sel$Sim50)))#,as.numeric(quantComb.sel[,4]),
#                           as.numeric(quantComb.sel[,5]),as.numeric(quantComb.sel[,6]),
#                           as.numeric(quantComb.sel[,7]))
#quant.max <- rowMaxs(as.matrix(quantComb.selMod[,c(1:7)]))
quant.max <- rowMaxs(as.matrix(quantComb.selMod[,c(1:3)]))
#Since some min values are positive and some are negative, do manually
#yMax <- quant.max + (quant.max*0.1)
yMax <- c(quant.max[1:8] + (quant.max[1:8]*0.1), quant.max[9] + (quant.max[9]*-0.1),
          quant.max[10:12] + (quant.max[10:12]*0.1), quant.max[13] + (quant.max[13]*-0.1),
          quant.max[14:24] + (quant.max[14:24]*0.1))


#quant10Comb.sel <- as.numeric(as.character(quant10Comb.sel[,1]))
quant10Comb.selMod <- cbind(as.numeric(as.character(quant10Comb.sel$Sim10)),as.numeric(as.character(quant10Comb.sel$Sim20)),
                            as.numeric(as.character(quant10Comb.sel$Sim50)))#,as.numeric(quant10Comb.sel[,4]),
                          #  as.numeric(quant10Comb.sel[,5]),as.numeric(quant10Comb.sel[,6]),
                          #  as.numeric(quant10Comb.sel[,7]))
quant10.min <- rowMins(as.matrix(quant10Comb.selMod[,c(1:3)]))

#Since some min values are positive and some are negative, do manually
yMin <- c(quant10.min[1:4] - (quant10.min[1:4]*0.1), quant10.min[5:16] - (quant10.min[5:16]*-0.1),
          quant10.min[17] - (quant10.min[17]*0.1), quant10.min[18:23] - (quant10.min[18:23]*-0.1),
          quant10.min[24] - (quant10.min[24]*0.1))
          



## Assign colors
colz <- c("black","grey40","grey80")


#saveRDS(medComb.sel, "20240117_medCombSel.rds")
#saveRDS(quantComb.sel, "20240117_quantCombSel.rds")
#saveRDS(quant10Comb.sel, "20240117_quant10CombSel.rds")
#medComb.sel <- readRDS("20240117_medCombSel.rds")
#quantComb.sel <- readRDS("20240117_quantCombSel.rds")
#quant10Comb.sel <- readRDS("20240117_quant10CombSel.rds")



## PLOT
#tiff('20230122_ErBr_fig2.tiff', res=400, pointsize=6, compression="lzw")
pdf('20240902_ErBr_FigSimDat.pdf', width=6.5, height=8.5)
par(mfrow=c(7,4), mar=c(1.25,2,1.9,2))  #Plot so 4 VR models are in cols and upto 7 predictor vars are in rows
#bottom, left, top, and right
for (nn in 1:17) {
  plot(c(1:ncol(medComb.sel)), c(as.numeric(as.character(medComb.sel[nn,1])),
                                 as.numeric(as.character(medComb.sel[nn,2])),
                                 as.numeric(as.character(medComb.sel[nn,3]))),
       col=colz, ylab=NA, xaxt = "n", main=names.paramTitlesOrd[nn], cex.axis=0.9,cex.main=0.9, 
       pch=19, cex=1.2, ylim=c(yMin[nn],yMax[nn]))
  abline(h=as.numeric(as.character(medRealDat.sel[nn,1])), col="grey30", lty="dotted")
  arrows(c(1:ncol(medComb.sel)),c(as.numeric(as.character(medComb.sel[nn,1])),
                                  as.numeric(as.character(medComb.sel[nn,2])),
                                  as.numeric(as.character(medComb.sel[nn,3]))),
         c(1:ncol(medComb.sel)), c(as.numeric(as.character(quantComb.sel[nn,1])),
                                   as.numeric(as.character(quantComb.sel[nn,2])),
                                   as.numeric(as.character(quantComb.sel[nn,3]))),
        lwd = 1.25, angle = 90, code = 3, length=0, col=colz)
  arrows(c(1:ncol(medComb.sel)),c(as.numeric(as.character(medComb.sel[nn,1])),
                                  as.numeric(as.character(medComb.sel[nn,2])),
                                  as.numeric(as.character(medComb.sel[nn,3]))),
         c(1:ncol(medComb.sel)), c(as.numeric(as.character(quant10Comb.sel[nn,1])),
                                   as.numeric(as.character(quant10Comb.sel[nn,2])),
                                   as.numeric(as.character(quant10Comb.sel[nn,3]))),
         lwd = 1.25, angle = 90, code = 3, length=0, col=colz)
  #arrows(1:6,as.numeric(medComb.sel[nn,1:6]),
   #      1:6, as.numeric(as.matrix(quantComb.sel[nn,1:6])),
    #     lwd = 1.25, angle = 90, code = 3, length=0, col=colz)
  #arrows(1:6, as.numeric(medComb.sel[nn,1:6]),
   #      1:6, as.numeric(as.matrix(quant10Comb.sel[nn,1:6])), 
    #     lwd = 1.25, angle = 90, code = 3, length=0, col=colz)
}
plot.new()
for (nn in 18:20) {
  plot(c(1:ncol(medComb.sel)), c(as.numeric(as.character(medComb.sel[nn,1])),
                                 as.numeric(as.character(medComb.sel[nn,2])),
                                 as.numeric(as.character(medComb.sel[nn,3]))), 
       col=colz, ylab=NA, xaxt = "n", main=names.paramTitlesOrd[nn], cex.axis=0.9,cex.main=0.9, 
       pch=19, cex=1.2, ylim=c(yMin[nn],yMax[nn]))
  #abline(h=0, col="grey80")
  abline(h=as.numeric(as.character(medRealDat.sel[nn,1])), col="grey30", lty="dotted")
  arrows(c(1:ncol(medComb.sel)),c(as.numeric(as.character(medComb.sel[nn,1])),
                                  as.numeric(as.character(medComb.sel[nn,2])),
                                  as.numeric(as.character(medComb.sel[nn,3]))),
         c(1:ncol(medComb.sel)), c(as.numeric(as.character(quantComb.sel[nn,1])),
                                   as.numeric(as.character(quantComb.sel[nn,2])),
                                   as.numeric(as.character(quantComb.sel[nn,3]))),
         lwd = 1.25, angle = 90, code = 3, length=0, col=colz)
  arrows(c(1:ncol(medComb.sel)),c(as.numeric(as.character(medComb.sel[nn,1])),
                                  as.numeric(as.character(medComb.sel[nn,2])),
                                  as.numeric(as.character(medComb.sel[nn,3]))),
         c(1:ncol(medComb.sel)), c(as.numeric(as.character(quant10Comb.sel[nn,1])),
                                   as.numeric(as.character(quant10Comb.sel[nn,2])),
                                   as.numeric(as.character(quant10Comb.sel[nn,3]))),
         lwd = 1.25, angle = 90, code = 3, length=0, col=colz)
}
plot.new()
for (nn in 21:24) {
  plot(c(1:ncol(medComb.sel)), c(as.numeric(as.character(medComb.sel[nn,1])),
                                 as.numeric(as.character(medComb.sel[nn,2])),
                                 as.numeric(as.character(medComb.sel[nn,3]))), 
       col=colz, ylab=NA, xaxt = "n", main=names.paramTitlesOrd[nn], cex.axis=0.9,cex.main=0.9, 
       pch=19, cex=1.2, ylim=c(yMin[nn],yMax[nn]))
  abline(h=as.numeric(as.character(medRealDat.sel[nn,1])), col="grey30", lty="dotted")
  arrows(c(1:ncol(medComb.sel)),c(as.numeric(as.character(medComb.sel[nn,1])),
                                  as.numeric(as.character(medComb.sel[nn,2])),
                                  as.numeric(as.character(medComb.sel[nn,3]))),
         c(1:ncol(medComb.sel)), c(as.numeric(as.character(quantComb.sel[nn,1])),
                                   as.numeric(as.character(quantComb.sel[nn,2])),
                                   as.numeric(as.character(quantComb.sel[nn,3]))),
         lwd = 1.25, angle = 90, code = 3, length=0, col=colz)
  arrows(c(1:ncol(medComb.sel)),c(as.numeric(as.character(medComb.sel[nn,1])),
                                  as.numeric(as.character(medComb.sel[nn,2])),
                                  as.numeric(as.character(medComb.sel[nn,3]))),
         c(1:ncol(medComb.sel)), c(as.numeric(as.character(quant10Comb.sel[nn,1])),
                                   as.numeric(as.character(quant10Comb.sel[nn,2])),
                                   as.numeric(as.character(quant10Comb.sel[nn,3]))),
         lwd = 1.25, angle = 90, code = 3, length=0, col=colz)
}
plot.new()

#legend("center", colnames(medComb.sel)[1:7], col=colz, pch=19, cex=1,
#       horiz=FALSE, bty="y",seg.len=1, xpd="NA")
legend("center", c("SimData 10yrs (no missing)","SimData 20yrs","SimData 50yrs"),
                   col=colz, pch=19, cex=1, horiz=FALSE, bty="y",seg.len=1, xpd="NA")

dev.off()
## -----------------------------------------------------------------







## LOAD MODEL SUMMARY OUTPUT -------------------------------------------------
#jags.mod <- readRDS("erbr_JAGSmodBest_SIM20yr_c3t5s10b5_noYRE_20240830.rds")
#summary(jags.mod)
#plot(jags.mod)
#summ.mod <- summary(jags.mod)
#gelman.diag(jags.mod, confidence = 0.95, transform=FALSE)

summ.mod1 <- readRDS("20240906_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat20yr.1.rds")
summ.mod2 <- readRDS("20240906_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat20yr.2.rds")
summ.mod3 <- readRDS("20240906_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat20yr.3.rds")
summ.mod4 <- readRDS("20240906_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat20yr.4.rds")
summ.mod5 <- readRDS("20240906_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat20yr.5.rds")
summ.mod6 <- readRDS("20240907_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat20yr.6.rds")
summ.mod7 <- readRDS("20240907_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat20yr.7.rds")
summ.mod8 <- readRDS("20240907_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat20yr.8.rds")
summ.mod9 <- readRDS("20240907_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat20yr.9.rds")
summ.mod10 <- readRDS("20240907_erbr_JAGSmodBestSUMM_c3t5s10b10_noYRE_SimDat20yr.10.rds")
## -----------------------------------------------------------------------------------------------



## Make plot comparing median param ests of rep simulated datasets & compare GLMM estimate w/out missing data
## Show median param ests of diff datasets as points and upper and lower 95% limits (from JAGS summary)

## SUBSET OUTPUT TO JUST KEEP PARAMS OF INTEREST
#names.param <- colnames(as.matrix(summ.mod1$mcmc[1]))[26:41]
names.param <- rownames(summ.mod1)[26:41]
names.param <- names.param[c(2:8,12:16)] #Remove intercept and GrwthVar  
names.paramTitles <- c("Grwth Size","Grwth Fall Temp","Grwth Summer Temp","Grwth Winter Temp",
                       "Grwth Fall Precip","Grwth Summer Precip","Grwth Winter Precip","Surv Size",
                       "Surv Winter Precip","Surv Fall Temp","Surv Summer Temp","Surv Winter Temp")

## Re-order parameter names for plotting 
#index<-c(1,5,3,7,11,9,13,2,14,6,4,8)
#names.paramOrd <- names.param[order(index)]
#names.paramTitlesOrd <- names.paramTitles[order(index)]
## -----------------------------------------------------------------------------------------------




## OBTAIN GLMM ESTIMATES -------------------------------------------------------------------------
## Loop over datasets 
paramsMM.grwth <- NULL
seMM.grwth <- NULL
paramsMM.surv <- NULL
seMM.surv <- NULL
n.datset <- 10

for (dd in 1:n.datset) {
  
  noMiss <- readRDS(file=paste("20240907_erbr_SimDat20yrNoMiss.",dd,".4GLM",".rds", sep=""))
  #noMiss2 <- readRDS("20240907_erbr_SimDat20yrNoMiss.2.4GLM.rds")
  
  ## Add t+1 climate, sz, & tag into erbr data 
  noMiss <- noMiss %>% mutate(TagNew1=lead(TagNew), RosNew1=lead(RosNew), Surv1=lead(surv))  
  #noMiss <- noMiss %>% mutate(PptFall1=lead(PptFall), PptWinter1=lead(PptWinter), PptSummer1=lead(PptSummer),
  #                         TempFall1=lead(TempFall), TempWinter1=lead(TempWinter), TempSummer1=lead(TempSummer))
  
  noMiss <- noMiss[which(noMiss$TagNew == noMiss$TagNew1),]  #Remove lines with mis-matched individuals 
  
  ## Log size in all models to match JAGs; this makes ending size a linear function of starting size
  
  ## Growth ** REORDER PARAMS FOR FUTURE RUNS ***
  glmm.grwth <- glmer.nb(RosNew1 ~ log(RosNew) + PptFall + PptWinter + PptSummer + 
                           TempFall + TempWinter + TempSummer + (1|TransectNew), data=noMiss)
  
  ## Survival  
  glmm.surv <- glmer(Surv1 ~ log(RosNew) + PptWinter + TempFall + (TempWinter) + 
                       (TempSummer) + (1|TransectNew), family=binomial(link='logit'), data=noMiss)
  
  
  ## Extract parameter estimates and SEs from GLMMs
  paramsMM.grwthTmp <- as.data.frame(summary(glmm.grwth)$coefficients[2:nrow(summary(glmm.grwth)$coefficients),1])
  seMM.grwthTmp <- as.data.frame(summary(glmm.grwth)$coefficients[2:nrow(summary(glmm.grwth)$coefficients),2])
  
  colnames(paramsMM.grwthTmp) <- paste("GLMM.",dd,sep="")
  colnames(seMM.grwthTmp) <- paste("SE.",dd,sep="")#),"ParamTitle")
  
  paramsMM.grwth <- as.data.frame(c(paramsMM.grwth, paramsMM.grwthTmp))
  seMM.grwth <- as.data.frame(c(seMM.grwth, seMM.grwthTmp))
  
  
  paramsMM.survTmp <- as.data.frame(summary(glmm.surv)$coefficients[2:nrow(summary(glmm.surv)$coefficients),1])#,
  seMM.survTmp <- as.data.frame(summary(glmm.surv)$coefficients[2:nrow(summary(glmm.surv)$coefficients),2])#,
  #c("Surv Size","Surv Winter Precip",
  # "Surv Fall Temp","Surv Winter Temp","Surv Summer Temp"))
  colnames(paramsMM.survTmp) <- paste("GLMM.",dd,sep="")
  colnames(seMM.survTmp) <- paste("SE.",dd,sep="")
  paramsMM.surv <- as.data.frame(c(paramsMM.surv, paramsMM.survTmp))
  seMM.surv <- as.data.frame(c(seMM.surv, seMM.survTmp))
  
}

paramsMM.grwth$ParamTitle <- c("Grwth Size","Grwth Fall Precip","Grwth Winter Precip","Grwth Summer Precip",
                               "Grwth Fall Temp","Grwth Winter Temp","Grwth Summer Temp")
seMM.grwth$ParamTitle <- c("Grwth Size","Grwth Fall Precip","Grwth Winter Precip","Grwth Summer Precip",
                           "Grwth Fall Temp","Grwth Winter Temp","Grwth Summer Temp")
paramsMM.surv$ParamTitle <- c("Surv Size","Surv Winter Precip",
                              "Surv Fall Temp","Surv Winter Temp","Surv Summer Temp")
seMM.surv$ParamTitle <- c("Surv Size","Surv Winter Precip",
                          "Surv Fall Temp","Surv Winter Temp","Surv Summer Temp")


## Load save glmm results
paramsMM.grwth <- readRDS(file=paste("20240908_erbr_paramMMgrwth_SimDat20yr",".rds", sep=""))
seMM.grwth <- readRDS(file=paste("20240908_erbr_seMMgrwth_SimDat20yr",".rds", sep=""))
paramsMM.surv <- readRDS(file=paste("20240908_erbr_paramMMsurv_SimDat20yr",".rds", sep=""))
seMM.surv <- readRDS(file=paste("20240908_erbr_seMMsurv_SimDat20yr",".rds", sep=""))

## Re-order parameter names for plotting 
paramsMM.grwthOrd <- paramsMM.grwth[match(names.paramTitles[1:7], paramsMM.grwth$ParamTitle),]
seMM.grwthOrd <- seMM.grwth[match(names.paramTitles[1:7], seMM.grwth$ParamTitle),]
paramsMM.survOrd <- paramsMM.surv[match(names.paramTitles[8:12], paramsMM.surv$ParamTitle),]
seMM.survOrd <- seMM.surv[match(names.paramTitles[8:12], seMM.surv$ParamTitle),]


#rbind(paramsMM.grwth, paramsMM.surv)
saveRDS(paramsMM.grwth, file=paste("20240908", "_erbr_paramMMgrwth_SimDat20yr", ".rds", sep=""))
saveRDS(seMM.grwth, file=paste("20240908", "_erbr_seMMgrwth_SimDat20yr", ".rds", sep=""))
saveRDS(paramsMM.surv, file=paste("20240908", "_erbr_paramMMsurv_SimDat20yr", ".rds", sep=""))
saveRDS(seMM.surv, file=paste("20240908", "_erbr_seMMsurv_SimDat20yr", ".rds", sep=""))
## ------------------------------------------------------------------------------------------------------------




## COMBINE MEDIAN PARAMETER VALUES ---------------------------------------------------------------------------
#summ.mod1$summaries
medParams.1 <- summ.mod1[c(27:33,37:41),2]
medParams.2 <- summ.mod2[c(27:33,37:41),2]
medParams.3 <- summ.mod3[c(27:33,37:41),2]
medParams.4 <- summ.mod4[c(27:33,37:41),2]
medParams.5 <- summ.mod5[c(27:33,37:41),2]
medParams.6 <- summ.mod6[c(27:33,37:41),2]
medParams.7 <- summ.mod7[c(27:33,37:41),2]
medParams.8 <- summ.mod8[c(27:33,37:41),2]
medParams.9 <- summ.mod9[c(27:33,37:41),2]
medParams.10 <- summ.mod10[c(27:33,37:41),2]

low95.1 <- summ.mod1[c(27:33,37:41),1]
low95.2 <- summ.mod2[c(27:33,37:41),1]
low95.3 <- summ.mod3[c(27:33,37:41),1]
low95.4 <- summ.mod4[c(27:33,37:41),1]
low95.5 <- summ.mod5[c(27:33,37:41),1]
low95.6 <- summ.mod6[c(27:33,37:41),1]
low95.7 <- summ.mod7[c(27:33,37:41),1]
low95.8 <- summ.mod8[c(27:33,37:41),1]
low95.9 <- summ.mod9[c(27:33,37:41),1]
low95.10 <- summ.mod10[c(27:33,37:41),1]

up95.1 <- summ.mod1[c(27:33,37:41),3]
up95.2 <- summ.mod2[c(27:33,37:41),3]
up95.3 <- summ.mod3[c(27:33,37:41),3]
up95.4 <- summ.mod4[c(27:33,37:41),3]
up95.5 <- summ.mod5[c(27:33,37:41),3]
up95.6 <- summ.mod6[c(27:33,37:41),3]
up95.7 <- summ.mod7[c(27:33,37:41),3]
up95.8 <- summ.mod8[c(27:33,37:41),3]
up95.9 <- summ.mod9[c(27:33,37:41),3]
up95.10 <- summ.mod10[c(27:33,37:41),3]

## COMBINE *****
## ------------------------------------------------------------------------------------------------




## PLOT GLMM ESTIMATES VS JAGS ESTIMATES ---------------------------------------------------------
#plot(paramsMM.grwth$GLMM[1], as.numeric(as.character(medComb$rep1[1])))
#abline(a=0, b=1)
#plot(as.numeric(paramsMM.grwth[1,c(1,4)]), c(as.numeric(medParams.1[1]),as.numeric(medParams.2[1])),
#     ylim=c(0.5,1), xlim=c(0.5,1))
#abline(a=0, b=1)
#c(5,2)
par(mfrow=c(1,1), mar=c(3.9,1.7,2.3,1.5))  #bottom, left, top and right 
par(pty="s")
#, pch=19

# Grwth Sz
plot(as.numeric(paramsMM.grwthOrd[1,1:10]), 
     c(as.numeric(medParams.1[1]),as.numeric(medParams.2[1]),
       as.numeric(medParams.3[1]),as.numeric(medParams.4[1]),
       as.numeric(medParams.5[1]),as.numeric(medParams.6[1]),
       as.numeric(medParams.7[1]),as.numeric(medParams.8[1]),
       as.numeric(medParams.9[1]),as.numeric(medParams.10[1])),
     ylim=c(0.5,0.9), xlim=c(0.5,0.9), xlab="GLMM estimate", 
     ylab="JAGS estimate", main=paramsMM.grwthOrd[1,11])
abline(a=0, b=1)
arrows(c(as.numeric(paramsMM.grwthOrd[1,1:10])),
       c(as.numeric(medParams.1[1]),as.numeric(medParams.2[1]),
         as.numeric(medParams.3[1]),as.numeric(medParams.4[1]),
         as.numeric(medParams.5[1]),as.numeric(medParams.6[1]),
         as.numeric(medParams.7[1]),as.numeric(medParams.8[1]),
         as.numeric(medParams.9[1]),as.numeric(medParams.10[1])),
       (c(as.numeric(paramsMM.grwthOrd[1,1:10])) + c(as.numeric(seMM.grwthOrd[1,1:10]))),
       lwd = 1.25, angle = 90, code = 3, length=0)
arrows(c(as.numeric(paramsMM.grwthOrd[1,1:10])),
       c(as.numeric(medParams.1[1]),as.numeric(medParams.2[1]),
         as.numeric(medParams.3[1]),as.numeric(medParams.4[1]),
         as.numeric(medParams.5[1]),as.numeric(medParams.6[1]),
         as.numeric(medParams.7[1]),as.numeric(medParams.8[1]),
         as.numeric(medParams.9[1]),as.numeric(medParams.10[1])),
       (c(as.numeric(paramsMM.grwthOrd[1,1:10])) - c(as.numeric(seMM.grwthOrd[1,1:10]))),
       lwd = 1.25, angle = 90, code = 3, length=0)
arrows(c(as.numeric(paramsMM.grwthOrd[1,1:10])),
       c(as.numeric(medParams.1[1]),as.numeric(medParams.2[1]),
         as.numeric(medParams.3[1]),as.numeric(medParams.4[1]),
         as.numeric(medParams.5[1]),as.numeric(medParams.6[1]),
         as.numeric(medParams.7[1]),as.numeric(medParams.8[1]),
         as.numeric(medParams.9[1]),as.numeric(medParams.10[1])),
       c(as.numeric(low95.1[1]),as.numeric(low95.2[1]),
         as.numeric(low95.3[1]),as.numeric(low95.4[1]),
         as.numeric(low95.5[1]),as.numeric(low95.6[1]),
         as.numeric(low95.7[1]),as.numeric(low95.8[1]),
         as.numeric(low95.9[1]),as.numeric(low95.10[1])),
       lwd=1.25, angle=90, code=3, length=0)

##points(mean_PCs$PC1score[5], mean_PCs$PC2score[5], col=cols[3], pch=symbs[3], cex=2.5) #post_non_d
## Add standard error bars 
##plotCI(x=mean_PCs$PC1score[1], y=mean_PCs$PC2score[1], uiw=PC1_SE[1,2],err="x", col=cols, add=T)
##plotCI(x=mean_PCs$PC1score[1], y=mean_PCs$PC2score[1], uiw=PC2_SE[1,2],err="y", col=cols, add=T)

# Surv Sz
plot(as.numeric(paramsMM.survOrd[1,1:10]), 
     c(as.numeric(medParams.1[8]),as.numeric(medParams.2[8]),
       as.numeric(medParams.3[8]),as.numeric(medParams.4[8]),
       as.numeric(medParams.5[8]),as.numeric(medParams.6[8]),
       as.numeric(medParams.7[8]),as.numeric(medParams.8[8]),
       as.numeric(medParams.9[8]),as.numeric(medParams.10[8])),
     ylim=c(0.2,0.8), xlim=c(0.2,0.8), xlab="GLMM estimate", 
     ylab="JAGS estimate", main=paramsMM.survOrd[1,11], pch=19)
abline(a=0, b=1)

# Grwth Fall Temp
plot(as.numeric(paramsMM.grwthOrd[2,1:10]), 
     c(as.numeric(medParams.1[2]),as.numeric(medParams.2[2]),
       as.numeric(medParams.3[2]),as.numeric(medParams.4[2]),
       as.numeric(medParams.5[2]),as.numeric(medParams.6[2]),
       as.numeric(medParams.7[2]),as.numeric(medParams.8[2]),
       as.numeric(medParams.9[2]),as.numeric(medParams.10[2])),
     ylim=c(-0.2,0.3), xlim=c(-0.2,0.3), xlab="GLMM estimate", 
     ylab="JAGS estimate", main=paramsMM.grwthOrd[2,11], pch=19)
abline(a=0, b=1)

# Surv Fall Temp
plot(as.numeric(paramsMM.survOrd[3,1:10]), 
     c(as.numeric(medParams.1[10]),as.numeric(medParams.2[10]),
       as.numeric(medParams.3[10]),as.numeric(medParams.4[10]),
       as.numeric(medParams.5[10]),as.numeric(medParams.6[10]),
       as.numeric(medParams.7[10]),as.numeric(medParams.8[10]),
       as.numeric(medParams.9[10]),as.numeric(medParams.10[10])),
     ylim=c(-1,1), xlim=c(-1,1), xlab="GLMM estimate", 
     ylab="JAGS estimate", main=paramsMM.survOrd[3,11], pch=19)
abline(a=0, b=1)

# Grwth Summer Temp
plot(as.numeric(paramsMM.grwthOrd[3,1:10]), 
     c(as.numeric(medParams.1[3]),as.numeric(medParams.2[3]),
       as.numeric(medParams.3[3]),as.numeric(medParams.4[3]),
       as.numeric(medParams.5[3]),as.numeric(medParams.6[3]),
       as.numeric(medParams.7[3]),as.numeric(medParams.8[3]),
       as.numeric(medParams.9[3]),as.numeric(medParams.10[3])),
     ylim=c(-0.2,0.3), xlim=c(-0.2,0.3), xlab="GLMM estimate", 
     ylab="JAGS estimate", main=paramsMM.grwthOrd[3,11], pch=19)
abline(a=0, b=1)

# Surv Summer Temp
plot(as.numeric(paramsMM.survOrd[4,1:10]), 
     c(as.numeric(medParams.1[11]),as.numeric(medParams.2[11]),
       as.numeric(medParams.3[11]),as.numeric(medParams.4[11]),
       as.numeric(medParams.5[11]),as.numeric(medParams.6[11]),
       as.numeric(medParams.7[11]),as.numeric(medParams.8[11]),
       as.numeric(medParams.9[11]),as.numeric(medParams.10[11])),
     ylim=c(-0.4,0.9), xlim=c(-0.4,0.9), xlab="GLMM estimate", 
     ylab="JAGS estimate", main=paramsMM.survOrd[4,11], pch=19)
abline(a=0, b=1)

# Grwth Winter Temp
plot(as.numeric(paramsMM.grwthOrd[4,1:10]), 
     c(as.numeric(medParams.1[4]),as.numeric(medParams.2[4]),
       as.numeric(medParams.3[4]),as.numeric(medParams.4[4]),
       as.numeric(medParams.5[4]),as.numeric(medParams.6[4]),
       as.numeric(medParams.7[4]),as.numeric(medParams.8[4]),
       as.numeric(medParams.9[4]),as.numeric(medParams.10[4])),
     ylim=c(-0.2,0.3), xlim=c(-0.2,0.3), xlab="GLMM estimate", 
     ylab="JAGS estimate", main=paramsMM.grwthOrd[4,11], pch=19)
abline(a=0, b=1)

# Surv Winter Temp
plot(as.numeric(paramsMM.survOrd[5,1:10]), 
     c(as.numeric(medParams.1[12]),as.numeric(medParams.2[12]),
       as.numeric(medParams.3[12]),as.numeric(medParams.4[12]),
       as.numeric(medParams.5[12]),as.numeric(medParams.6[12]),
       as.numeric(medParams.7[12]),as.numeric(medParams.8[12]),
       as.numeric(medParams.9[12]),as.numeric(medParams.10[12])),
     ylim=c(-0.4,0.7), xlim=c(-0.4,0.7), xlab="GLMM estimate", 
     ylab="JAGS estimate", main=paramsMM.survOrd[5,11], pch=19)
abline(a=0, b=1)

# Grwth Winter Precip
plot(as.numeric(paramsMM.grwthOrd[7,1:10]), 
     c(as.numeric(medParams.1[7]),as.numeric(medParams.2[7]),
       as.numeric(medParams.3[7]),as.numeric(medParams.4[7]),
       as.numeric(medParams.5[7]),as.numeric(medParams.6[7]),
       as.numeric(medParams.7[7]),as.numeric(medParams.8[7]),
       as.numeric(medParams.9[7]),as.numeric(medParams.10[7])),
     ylim=c(-0.01,0.02), xlim=c(-0.01,0.02), xlab="GLMM estimate", 
     ylab="JAGS estimate", main=paramsMM.grwthOrd[7,11], pch=19)
abline(a=0, b=1)

# Surv Winter Precip
plot(as.numeric(paramsMM.survOrd[2,1:10]), 
     c(as.numeric(medParams.1[9]),as.numeric(medParams.2[9]),
       as.numeric(medParams.3[9]),as.numeric(medParams.4[9]),
       as.numeric(medParams.5[9]),as.numeric(medParams.6[9]),
       as.numeric(medParams.7[9]),as.numeric(medParams.8[9]),
       as.numeric(medParams.9[9]),as.numeric(medParams.10[9])),
     ylim=c(-0.015,0.05), xlim=c(-0.015,0.05), xlab="GLMM estimate", 
     ylab="JAGS estimate", main=paramsMM.survOrd[2,11], pch=19)
abline(a=0, b=1)
## ----------------------------------------------------------------

