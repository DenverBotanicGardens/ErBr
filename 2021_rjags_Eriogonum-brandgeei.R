library(R2jags)
library(dplyr)
library(bazar)

setwd("C:/Users/DePrengm/Denver Botanic Gardens/Conservation - General/AllProjectsBySpecies/Eriogonum brandegeei/2020_Eriogonum-brandegeei_AprilGoebl_PVA/2021_Eriogonum-brandegeei_integratedpopulationmodels")
dats <- read.csv("erbr_TagClust_2021-08-11.csv")

# Scale predictors for better model convergence
dats <- dats %>% mutate(RosNew = scale(RosNew),
                         PptFall = scale(PptFall),
                         PptWinter = scale(PptWinter),
                         PptSummer = scale(PptSummer),
                         TempFall = scale(TempFall),
                         TempWinter = scale(TempWinter),
                         TempSummer = scale(TempSummer))

lagforsurv <- dats$lagforsurv       #Full list of lags or of -1 for first observation rows; Fixed so first are -1 
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
numyears <- length(unique(dats$Year))
numtrans <- length(unique(dats$TransectNew))
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


## Make dataframe w new plts (that are likely recent seedlings) for each transect & yr ------------
## Make df that will hold data containing new plants
years <- unique(dats$Year.num)
years <- years[order(years)]
dats.newPlts <- as.data.frame(rep(unique(dats$TransectNew.num), each=length(years)))
colnames(dats.newPlts) <- "TransectNew.num"
dats.newPlts$Year.num <- rep(years)

## Identify new plants
newPlts <- dats %>% group_by(TagNew) %>% slice(which.min(Year))   #Identify rows with 1st appearance for each plt
newPlts <- data.frame(newPlts)
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
newplt.trans <- dats.newPlts$TransectNew.num          # to make newplt.yrtranscombo
newplt.yr <- dats.newPlts$Year.num                    # to make newplt.yrtranscombo
newPltlines <- length(dats.newPlts$TransectNew.num)

## Make a linear index of transect-year combos for new plts
newplt.yrtranscombo=100*newplt.trans+newplt.yr 

# needs data frame dats and the lines and lengths for loops
jags_erbrdata <- with(dats, list(Site = Site,Year=Year,RosNew=RosNew,InflNew=InflNew,TagNew=TagNew,surv=surv,lagsrtsz=lagsrtsz,
                                 lagforsurv=lagforsurv,PptFall=PptFall,PptWinter=PptWinter,PptSummer=PptSummer,
                                 TempFall=TempFall,TempWinter=TempWinter,TempSummer=TempSummer,TransectNew=TransectNew,
                                 InflYesNo=InflYesNo,RowNum=RowNum,TransectNew.num=TransectNew.num,Year.num=Year.num, 
                                 Survs = Survs, RosNew = RosNew, InflNew = InflNew, InflYesNo = InflYesNo,
                                 goodrows = goodrows, goodgrowrows = goodgrowrows,lagvals = lagvals, 
                                 # Reproduction 
                                 rows.w.sz = rows.w.sz, rows.wo.sz = rows.wo.sz, rows.w.inflors= rows.w.inflors, datayesno = datayesno,
                                 # Survival
                                 rows.wo.sz.alive =rows.wo.sz.alive, 
                                 # Recruitment
                                 newPltlines = newPltlines, newplts = newplts,
                                 # Random effects
                                 datayesno= datayesno, Year.num =Year.num, TransectNew.num = TransectNew.num, yrtranscombo = yrtranscombo, 
                                 newplt.yrtranscombo = newplt.yrtranscombo, 
                                 # Loop lengths
                                 Nallrows = length(dats$Site), Nyears = length(unique(dats$Year)), Ntrans = length(unique(dats$TransectNew)),
                                 Ncases = Ncases, Ngrowcases = Ngrowcases, 
                                 Ndirectszcases = Ndirectszcases, Nindirectszcases = Nindirectszcases,
                                 Nrows.w.inflors = length(rows.w.inflors)))


# note that precision (1/variance) 
# ~ to denote that on left is a random variables distributed along the distribution on the right
# first check smaller models in jags_Eriogonum-brande

lm1_jags <- function(){
  # Likelihood
  
  ##############################################################################################
  # Growth ________________________________________________________________________________________________________
  # Size for each plant; estimate survival to each year: chance of survival for unobserved years in the past
  for(i in 1:Ncases){
    # lag from size (RosNew)
    regression_mean[goodrows[i] - lagvals[i]+1] <- exp(grwth_intercept + grwth_RosCoef*log(RosNew[goodrows[i]-lagvals[i]]) +
                                                         grwth_TempFallCoef*TempFall[goodrows[i] - lagvals[i]+1] +
                                                         grwth_Transect_randomeffect[TransectNew.num[goodrows[i] - lagvals[i]]] +
                                                         grwth_Year_randomeffect[Year.num[goodrows[i] - lagvals[i]]])
    
    r.growth[goodrows[i] - lagvals[i]+1] <- exp(grwthvar_intercept + grwthvar_RosCoef*log(RosNew[goodrows[i] - lagvals[i]]))
    
    
  # Survival_______________________________________________________________________________________________________
  # or could write as "logit(Surv_mu) ~ "instead of "~ 1/(1+exp(-lm))"
    Surv_mu[goodrows[i]-lagvals[i]+1] <- 1/(1+exp(-(surv_intercept + surv_RosCoef*log(RosNew[goodrows[i]-lagvals[i]]) + 
                                                      surv_PptWinterCoef*PptWinter[goodrows[i]-lagvals[i]+1] + 
                                                      surv_TempWinterCoef*TempWinter[goodrows[i]-lagvals[i]+1]  + 
                                                      surv_TempFallCoef*TempFall[goodrows[i]-lagvals[i]+1] + 
                                                      surv_TempSummerCoef*TempSummer[goodrows[i]-lagvals[i]+1]+ 
                                                      surv_Transect_randomeffect[TransectNew.num[goodrows[i]-lagvals[i]]] + 
                                                      surv_Year_randomeffect[Year.num[goodrows[i]-lagvals[i]]])))
    for (j in (goodrows[i]-lagvals[i]+1):(goodrows[i]-1)) {
      regression_mean[j+1] <- exp(grwth_intercept + grwth_RosCoef*log(regression_mean[j]) + 
                                    grwth_TempFallCoef*TempFall[j+1] + 
                                    grwth_Transect_randomeffect[TransectNew.num[j]] + 
                                    grwth_Year_randomeffect[Year.num[j]])
      
      r.growth[j+1] <- exp(grwthvar_intercept + grwthvar_RosCoef*log(regression_mean[j]))
      Surv_mu[j+1] <- Surv_mu[j]*1/(1+exp(-(surv_intercept + surv_RosCoef*log(regression_mean[j]) + 
                                              surv_PptWinterCoef*PptWinter[j+1] + 
                                              surv_TempWinterCoef*TempWinter[j+1] + 
                                              surv_TempFallCoef*TempFall[j+1] + 
                                              surv_TempSummerCoef*TempSummer[j+1]+ 
                                              surv_Transect_randomeffect[TransectNew.num[j]] + 
                                              surv_Year_randomeffect[Year.num[j]])))
      } #End lags loop
    } #End of going through cases for growth and survival
  ##############################################################################################
    
  # Reproduction_____________________________________________________
  for(i in 1:Ndirectszcases){            # Based on observed sizes
    InflYesNo[rows.w.sz[i]] ~ dbern(repro_prob[rows.w.sz[i]])
    
    logit(repro_prob[rows.w.sz[i]]) <- reproyesno_intercept + reproyesno_RosCoef*log(RosNew[rows.w.sz[i]]) +
      reproyesno_TempWinterCoef*TempWinter[rows.w.sz[i]] + 
      reproyesno_PptFallCoef*PptFall[rows.w.sz[i]] +
      reproyesno_Transect_randomeffect[TransectNew.num[rows.w.sz[i]]] + 
      reproyesno_Year_randomeffect[Year.num[rows.w.sz[i]]] 
    
    repro_amount[rows.w.sz[i]] <- exp(repro_intercept + repro_RosCoef*log(RosNew[rows.w.sz[i]]) + 
                                        repro_TempFallCoef*TempFall[rows.w.sz[i]] + 
                                        repro_PptSummerCoef*PptSummer[rows.w.sz[i]] + 
                                        repro_TempWinterCoef*TempWinter[rows.w.sz[i]] +
                                        repro_Transect_randomeffect[TransectNew.num[rows.w.sz[i]]] + 
                                        repro_Year_randomeffect[Year.num[rows.w.sz[i]]]) 
    
    inflors[rows.w.sz[i]] = InflNew[rows.w.sz[i]]
  }
  
  for(i in 1:Nindirectszcases){         # Based on inferred sizes, from regression_mean of size (Growth) instead of RosNew
    logit(repro_prob[rows.wo.sz[i]]) <- reproyesno_intercept + reproyesno_RosCoef*log(regression_mean[rows.wo.sz[i]]) +
      reproyesno_TempWinterCoef*TempWinter[rows.wo.sz[i]] + 
      reproyesno_PptFallCoef*PptFall[rows.wo.sz[i]] + 
      reproyesno_Transect_randomeffect[TransectNew.num[rows.wo.sz[i]]] + 
      reproyesno_Year_randomeffect[Year.num[rows.wo.sz[i]]]
    
  repro_amount[rows.wo.sz[i]] <- exp(repro_intercept + repro_RosCoef*log(regression_mean[rows.wo.sz[i]]) + 
                                       repro_TempFallCoef*TempFall[rows.wo.sz[i]] +  
                                       repro_PptSummerCoef*PptSummer[rows.wo.sz[i]] + 
                                       repro_TempWinterCoef*TempWinter[rows.wo.sz[i]] +
                                       repro_Transect_randomeffect[TransectNew.num[rows.wo.sz[i]]] + 
                                       repro_Year_randomeffect[Year.num[rows.wo.sz[i]]]) 
    
  ## Max statement enforces that if a plt must have been alive, 
  # its estimated infs should be counted in the tot for estimating new plts - new plant contributions?
  inflors[rows.wo.sz[i]] = max(Surv_mu[rows.wo.sz[i]], rows.wo.sz.alive[i])*repro_prob[rows.wo.sz[i]]*repro_amount[rows.wo.sz[i]]
  }
  
  
  # compare estimated inflorescence number to observed inflorescence number
  for(i in 1:Nrows.w.inflors){
    p.infls[i] <- r.infls/(r.infls + repro_amount[rows.w.inflors[i]])
    InflNew[rows.w.inflors[i]] ~ dnegbin(p.infls[i], r.infls)
  }
  
  for(i in 1:Ngrowcases){
    regression_residual[i] <- RosNew[goodgrowrows[i]] - regression_mean[goodgrowrows[i]]
    p[goodgrowrows[i]] <- r.growth[goodgrowrows[i]]/(r.growth[goodgrowrows[i]]+regression_mean[goodgrowrows[i]])
    RosNew[goodgrowrows[i]] ~ dnegbin( p[goodgrowrows[i]], r.growth[goodgrowrows[i]])
  } #End of going through cases for sizes
  
  for(i in 1:Ncases){
    ## Fit the survival function:
    Survs[goodrows[i]] ~ dbern(Surv_mu[goodrows[i]])
  } #End loop to do survival predictions

  ## total number of inflorescence and new plants the next year
  for(i in 1:newPltlines){
    pred.tot.inflors[i] = inprod(yrtranscombo==newplt.yrtranscombo[i], inflors) # The inprod & the 1st logical are making 0,1 for whether a given line is part of each trans-yr combo
    mn.new.plts[i] = exp(newplt_intercept + log(pred.tot.inflors[i] + 0.1))
    p.newplts[i] <- r.newplts/(r.newplots+mn.new.plts[i])
    newplots[i] ~ dnegbin(p.newplts[i], r.newplts)
  }
    
  # Priors
  
  ## Growth
    grwth_intercept ~ dunif(-3,3)
    grwth_RosCoef ~ dunif(-3,3)
    grwth_TempFallCoef ~ dunif(-2,2)
    grwthvar_intercept ~ dunif(-3,3)
    grwthvar_RosCoef ~ dunif(-3,3)
    
    for(TransectNew.num_iterator in 1:numtrans){
      grwth_Transect_randomeffect[TransectNew.num_iterator] ~ dnorm(0, grwth_Transect_precision)
    }
    grwth_Transect_sigma ~ dunif(0,100)                                              # standard deviation 
    grwth_Transect_precision ~ 1 / (grwth_Transect_sigma * grwth_Transect_sigma)     # sigma^2 for tau
    
    for(Year_iterator in 1:numyears){
      grwth_Year_randomeffect[Year_iterator] ~ dnorm(0, grwth_Year_precision)
    }
    grwth_Year_sigma ~ dunif(0,100)                                                  # standard deviation 
    grwth_Year_precision ~ 1 / (grwth_Year_sigma * grwth_Year_sigma)                 # sigma^2 for tau
    
    
  ## Survival  
    surv_intercept ~ dnorm(0, 0.00001)  # dbeta(1,1)
    surv_RosCoef ~ dnorm(0, 0.00001)
    surv_PptWinterCoef ~ dnorm(0, 0.00001)
    surv_TempFallCoef ~ dnorm(0, 0.00001)
    surv_TempSummerCoef ~ dnorm(0, 0.00001) 
    surv_TempWinterCoef ~ dnorm(0, 0.00001)
    
    for(TransectNew.num_iterator in 1:numtrans){
      surv_Transect_randomeffect[TransectNew.num_iterator] ~ dnorm(0, surv_Transect_precision)
    }
    surv_Transect_sigma ~ dunif(0,100)                                            # standard deviation 
    surv_Transect_precision ~ 1 / (surv_Transect_sigma * surv_Transect_sigma)     # sigma^2 for tau
    
    for(Year_iterator in 1:numyears){
      surv_Year_randomeffect[Year_iterator] ~ dnorm(0, surv_Year_precision)
    }
    surv_Year_sigma ~ dunif(0,100)                                                # standard deviation 
    surv_Year_precision ~ 1 / (surv_Year_sigma * surv_Year_sigma)                 # sigma^2 for tau
    
    
  ## Probability of reproduction
    reproyesno_intercept ~  dnorm(0, 0.00001)
    reproyesno_RosCoef ~  dnorm(0, 0.00001)
    reproyesno_PptFallCoef ~  dnorm(0, 0.00001)
    reproyesno_TempWinterCoef ~  dnorm(0, 0.00001)
    
    for(TransectNew.num_iterator in 1:numtrans){
      reproyesno_Transect_randomeffect[TransectNew.num_iterator] ~ dnorm(0, reproyesno_Transect_precision)
    }
    reproyesno_Transect_sigma ~ dunif(0,100)                                                        # standard deviation 
    reproyesno_Transect_precision ~ 1 / (reproyesno_Transect_sigma * reproyesno_Transect_sigma)     # sigma^2 for tau
    
    for(Year_iterator in 1:numyears){
      reproyesno_Year_randomeffect[Year_iterator] ~ dnorm(0, reproyesno_Year_precision)
    }
    reproyesno_Year_sigma ~ dunif(0,100)                                                            # standard deviation 
    reproyesno_Year_precision ~ 1 / (reproyesno_Year_sigma * reproyesno_Year_sigma)                 # sigma^2 for tau
    
    
  ## Amount of reproduction
    r.infls ~ dunif(0,10) 
    repro_intercept ~ dunif(-5,5) 
    repro_RosCoef ~ dunif(0.5,1.5) 
    repro_PptSummerCoef ~ dnorm(0, 0.00001)
    repro_TempWinterCoef ~ dnorm(0, 0.00001)   
    repro_TempFallCoef ~ dnorm(0, 0.00001)     
    
    
    for(TransectNew.num_iterator in 1:numtrans){
      repro_Transect_randomeffect[TransectNew.num_iterator] ~ dnorm(0, repro_Transect_precision)
    }
    repro_Transect_sigma ~ dunif(0,100)                                              # standard deviation 
    repro_Transect_precision ~ 1 / (repro_Transect_sigma * repro_Transect_sigma)     # sigma^2 for tau
    
    for(Year_iterator in 1:numyears){
      repro_Year_randomeffect[Year_iterator] ~ dnorm(0, repro_Year_precision)
    }
    repro_Year_sigma ~ dunif(0,100)                                                  # standard deviation 
    repro_Year_precision ~ 1 / (repro_Year_sigma * repro_Year_sigma)                 # sigma^2 for tau
    
    # Negative binomial variance for predicted inflorescence number
    ## *** Why need this? never used is it?
    inflors_precision ~ dgamma(0.001, 0.001)
    var.inflors = 1/inflors_precision
    
    
    ## NEW PLANTS
    newplt_intercept ~ dnorm(0, 10^-6)
    r.newplts ~ dgamma(0.01,0.01)
    
    
    resid.sum.sq <- sum(regression_residual^2)
}


init_values <- function(){
  list(grwth_intercept = runif(1),
       grwth_RosCoef = runif(1),
       grwth_TempFallCoef = runif(1),
       grwth_Year_precision = runif(17),
       grwth_Transect_precision = runif(12),
       grwthvar_intercept = runif(1), 
       grwthvar_RosCoef = runif(1),
       
       surv_intercept = runif(1), #rbeta(1,0.5),
       surv_RosCoef = rnorm(1),
       surv_PptWinterCoef = rnorm(1),
       surv_TempWinterCoef = rnorm(1),
       surv_TempFallCoef = rnorm(1),
       surv_TempSummerCoef = rnorm(1),
       surv_Transect_precision = runif(12),
       surv_Year_precision= runif(17),
       
       reproyesno_intercept = runif(1), #rbeta(1,0.5), # conjugate of a binomial or Bernoulli is beta
       reproyesno_RosCoef = rnorm(1),
       reproyesno_TempWinterCoef = rnorm(1),
       reproyesno_PptFallCoef = rnorm(1),
       reproyesno_Year_precision= runif(17),
       reproyesno_Transect_precision = runif(12),
       
       repro_precision = runif(1),
       repro_intercept = rgamma(1,2),       # conjugate of Poisson is gamma
       repro_RosCoef = runif(1),
       repro_PptSummerCoef = rnorm(1),
       repro_TempWinterCoef = rnorm(1),
       repro_TempFallCoef = rnorm(1),
       repro_Year_precision =runif(17),
       repro_Transect_precision = runif(12),
       
       newplt_intercept = rnorm(1),
       r.newplts = rgamma(1,2))
}

params <- c("grwth_Transect_randomeffect","grwth_Year_randomeffect","surv_Transect_randomeffect","surv_Year_randomeffect",
"reproyesno_Transect_randomeffect","reproyesno_Year_randomeffect","repro_Transect_randomeffect","repro_Year_randomeffect",
"repro_Transect_precision","repro_Year_precision","reproyesno_Transect_precision","reproyesno_Year_precision","surv_Transect_precision",
"surv_Year_precision","grwth_Transect_precision","grwth_Year_precision","grwth_intercept","grwth_RosCoef",
"grwth_TempFallCoef","grwthvar_intercept","grwthvar_RosCoef","surv_intercept","surv_RosCoef","surv_PptWinterCoef",
"surv_TempFallCoef","surv_TempSummerCoef","surv_TempWinterCoef","reproyesno_intercept","reproyesno_RosCoef",
"reproyesno_PptFallCoef","reproyesno_TempWinterCoef","repro_intercept","repro_RosCoef","repro_PptSummerCoef",
"repro_TempWinterCoef","repro_TempFallCoef","newplt_intercept")


fit_lm1 <- jags(data = jags_erbrdata, inits = init_values, parameters.to.save = params, model.file = lm1_jags,
                n.chains = 3, n.iter = 12000, n.burnin = 2000, n.thin = 10, DIC = F)



