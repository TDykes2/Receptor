#Code by Alissa Gordon and Tyler Dykes

## Null and Delta Code ###############

## Required Packages ##
library(e1071)
library(survival)
library(survminer)
library(data.table)
#install.packages('BSDA')
library(BSDA)
##install.packages('fpp3')
library(fpp3)
##install.packages('survRM2')
library(survRM2)
##install.packages('dplyr')
library(dplyr)


T1Reps <- data.frame("RCST" = NA,"KM"= NA,"ST"= NA,"COX"= NA, "Recist"= NA, "KM"= NA, "ST"= NA,"Cox"= NA, "Continue?"= NA)


for (p in 1:1000) {
  


## Key inputs ##
n<- 140 #Using N=170 to distribute the amount of patients evening between both treatments.
pmedian<- 3 #Primary Median Lifespan                        ### Change both to 2.5 for NULL
mmedian<- 2 #Metastatic Median Lifespan                     ### Change both to 2.5 for NULL
plambda<-log(2)/pmedian #finds lambda for which median of exp curve would be pmedian
mlambda<-log(2)/mmedian #finds lambda for which median of exp curve would be mmedian
alpha<-0.1 #Number of times to run entire simulation
yearsrecruiting<-3 #years spent recruiting participants
yearsoftrial<-7 #years of trial
d.sorted <- NA





## Subgroup creation ##
N<-2*n #uses sample size that includes both treatments
subgroup<-0 #establish variable subgroup
#g1-g12 used to count, start off at 0
g1<-0 
g2<-0
g3<-0
g4<-0
g5<-0
g6<-0
g7<-0
g8<-0
g9<-0
g10<-0
g11<-0
g12<-0
allen= data.frame("subgroup"=0) ################
i=0

while(g1<N||g2<N||g3<N||g4<N){
  #subgroup each loop determined by random discrete probability from the given probabilities
  
  
  
  subgroup<-rdiscrete(1, probs=c(0.4, 0.1,0.1,0.1,0.06,0.06,0.05,0.03,0.03,0.03,0.02,0.02))
  
  #for groups 1-4, value is only counted in group size if group size is still less than N
  if (subgroup==1 & g1<N){
    g1<-g1+1; i=i+1 ;allen[i,]= subgroup
  }
  if (subgroup==2 & g2<N){
    g2<-g2+1; i=i+1 ;allen[i,]= subgroup
  }
  if (subgroup==3 & g3<N){
    g3<-g3+1; i=i+1 ;allen[i,]= subgroup
  }
  if (subgroup==4 & g4<N){
    g4<-g4+1; i=i+1 ;allen[i,]= subgroup
  }
  
  #for groups 5-12, value is always counted in group size
  if (subgroup==5){
    g5<-g5+1; i=i+1 ;allen[i,]= subgroup
  }
  if (subgroup==6){
    g6<-g6+1; i=i+1 ;allen[i,]= subgroup
  }
  if (subgroup==7){
    g7<-g7+1; i=i+1 ;allen[i,]= subgroup
  }
  if (subgroup==8){
    g8<-g8+1; i=i+1 ;allen[i,]= subgroup
  }
  if (subgroup==9){
    g9<-g9+1; i=i+1 ;allen[i,]= subgroup
  }
  if (subgroup==10){
    g10<-g10+1; i=i+1 ;allen[i,]= subgroup
  }
  if (subgroup==11){
    g11<-g11+1; i=i+1 ;allen[i,]= subgroup
  }
  if (subgroup==12){
    g12<-g12+1; i=i+1 ;allen[i,]= subgroup
  }
  
  
  
}

groups<-data.frame(g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12)
ototal<-sum(groups) #finds grand total of all subgroup sizes


## if odd add an extra one to g12 and one to total
if((ototal%%2)!=0){
  g12<- g12+1
  ototal<- ototal+1
  allen[nrow(allen) + 1, ]<- c(12)
  groups<-data.frame(g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12)
}



## Data frame creation and censoring ##
total<- ceiling(ototal/2)



arm<-rep(c(1,2),each=total) # primary=1 and metastatic=2; n rows created for each treatment (2n total rows)



enrollment<-runif(total*2, min=0, max=yearsrecruiting) #enrollment randomly generated from a uniform distribution from 0 to yearsrecruiting



survival<-ifelse (arm==1, rexp(total, plambda), rexp(total, mlambda)) #randomly generates survival times based off exponential model and lambda calculated for each group

d<- data.frame(survival,arm, enrollment)
d['newvalue']<-NA
d['status1']<-NA #1=death 0=censored
d['status2']<-NA
d['status']<-NA 
d['time']<-NA
d['censorP']<-NA
d['censorM']<-NA


### CENSORING ###

##establishing censor rate of around 10%##
censoring<- total%/%5 #20% of data potentially censored (will be around 10% after process); evenly rounded in order to use in sample size estimation (decimals interfere with dataframe size)
censoringp<-survival[1:censoring] #the first 20% of primary data marked for potential censoring
censoringm<-survival[total+1:censoring] #the first 20% of metastatic data marked for potential censoring

###creates a new column that generates numbers for the marked data from same exponential dataset##
for (val in censoringp) {
  censornumberp<-rexp(censoring, plambda)
} 
for (val in censoringm){
  censornumberm<-rexp(censoring,mlambda)
}

d$censorP[1:censoring] <- censornumberp  ## Quality Control:
d$censorM[total+1:censoring] <- censornumberm  ## Quality Control: 


##random censoring##
for (i in 1:length(censoringp)){
  ifelse(censoringp[i]<censornumberp[i], d$newvalue[i]<-survival[i], d$newvalue[i]<-censornumberp[i])
} #if the potential censor value is less than the original survival time, subject is censored.  This change in time or lack of change is represented by newvalue

for (i in 1:length(censoringm)){
  ifelse(censoringm[i]<censornumberm[i], d$newvalue[total+i]<-survival[total+i], d$newvalue[total+i]<-censornumberm[i])
}#same process but with metastatic patients

for (i in 1:length(survival)){
  if (is.na(d$newvalue[i])){
    d$newvalue[i]<-survival[i]
  }
}
##marks censored data based off of randomness##
for (i in 1:length(survival)){
  ifelse(d$newvalue[i]<survival[i], d$status1[i]<-0, d$status1[i]<-1)
} 

d$time<-d$newvalue+d$enrollment

##censors if time is more than yearsoftrial##
for (i in 1: length(survival)){
  ifelse(d$time[i]>yearsoftrial, d$status2[i]<-0, d$status2[i]<-1)
} #marks as censored in status2 if time>yearsoftrial

for(i in 1:length(survival)){
  ifelse(d$status2[i]==0,d$newvalue[i]<-yearsoftrial-enrollment[i],d$newvalue[i]<-d$newvalue[i])
} #if censored, newvalue (variable being analyzed in km curve) is yearsoftrial minus enrollment, otherwise newvalue remains initial survival time


##marks as censored in new variable status if censored either due to time or randomness##
for(i in 1:length(survival)){
  ifelse(d$status1[i]==0|d$status2[i]==0, d$status[i]<-0, d$status[i]<-1)
}


##d.sorted creation and recist function##
d.sorted<- data.frame(d)
d.sorted<- d.sorted[order(d.sorted$enrollment),] #sorts from earliest enrolled to latest


#5 and 12 subgroups columns created
d.sorted$twelvegroups <- allen$subgroup
d.sorted$fivegroups <- d.sorted$twelvegroups
d.sorted$fivegroups <- ifelse(d.sorted$fivegroups >= 5, 5, d.sorted$fivegroups)



clock <- d.sorted$enrollment[N]
cut <- clock +.25
d.sorted <- d.sorted[which(d.sorted$enrollment < cut),]


d.sortedgrowth <- d.sorted[d.sorted$enrollment < clock,] ; d.sortedgrowth <- d.sortedgrowth[d.sortedgrowth$status1 == 1,]
d.sortedgrowth$g <- rdiscrete(nrow(d.sortedgrowth), probs=c(0.25,0.25,0.25,0.25)) #generates 1's (Disappeared), 2's (decrease), 3's (No change), 4's (Increase)
d.sortedgrowth$g[d.sortedgrowth$newvalue<=clock] = 0 #Changes values in growth column who did not last 3 months to 0 representing "death" or "lost contact"

rchisq <- chisq.test(d.sortedgrowth$g, d.sortedgrowth$arm) #Performs chi sq test explaining growth by group
RST <- ifelse(rchisq$p.value <= 0.05, "Sig", "NotSig")



##New Censoring
d.sorted$censortest <- d.sorted$newvalue + d.sorted$enrollment
d.sorted$newvalue <- ifelse(d.sorted$censortest > cut,
                            cut-d.sorted$enrollment, d.sorted$newvalue)

d.sorted$statusC <- ifelse(d.sorted$censortest > cut, 0, d.sorted$status)



## km curve output and p value ##
kmcurve<-survfit(Surv(newvalue, statusC)~arm, data=d.sorted[which (d.sorted$fivegroups == 1),])  #creates curve model of newvalue, censored by status, and sorted by group
ggsurvplot(kmcurve, pval=TRUE, break.x.by=0.1, surv.median.line = "v", risk.table = "nrisk_cumevents", cumcensor=TRUE, xscale=1)
## survdiff(Surv(d.sorted$newvalue,d.sorted$statusC)~d.sorted$arm) #plots the curve and shows risk, event, and censor table by year
survdifffxn<- survdiff(Surv(d.sorted$newvalue, d.sorted$statusC)~d.sorted$arm)
KMP <- 1 - pchisq(survdifffxn$chisq, length(survdifffxn$n) - 1) #finds p value using chi square
KM <- ifelse(KMP <= 0.05, "Sig", "NotSig")
## Needs to just be on the first group



## Cox Regression ##
cox.mod <- coxph( Surv(d.sorted$newvalue,d.sorted$statusC) ~ d.sorted$arm*d.sorted$fivegroups)
mallon<- summary(cox.mod)
IT<- ifelse(mallon$coefficients[3,5] <= 0.05, "Sig", "NotSig")
MCox<- mallon$coefficients[3,5]




## Sign Test ##
STBase<- aggregate(d.sorted, by = list(d.sorted$arm, d.sorted$twelvegroups),FUN = mean) %>% 
  select(Group.1,Group.2,survival)
ST <- data.frame(c((STBase[1,3]-STBase[2,3]),(STBase[3,3]-STBase[4,3]), (STBase[5,3]-STBase[6,3]),
                   (STBase[7,3]-STBase[8,3]), (STBase[9,3]-STBase[10,3]),(STBase[11,3]-STBase[12,3]),
                   (STBase[13,3]-STBase[14,3]),(STBase[15,3]-STBase[16,3]),(STBase[17,3]-STBase[18,3]),
                   (STBase[19,3]-STBase[20,3]),(STBase[21,3]-STBase[22,3]),(STBase[23,3]-STBase[24,3])))
ST$Sign <- ifelse(ST[,1]>=0,"+","-")
STTest <- ifelse(max(table(ST$Sign)) > 10,"Sig", "NotSig")


## Key Output ##
T1<- data.frame(rchisq$p.value, KMP, max(table(ST$Sign)), MCox, RST, KM, STTest, IT)


  
  ## Logic  NEW based initial runs
  T1$Continue <- ifelse(IT == "Sig", "Continue",
                        ifelse(IT == "NotSig" && STTest == "Sig" , "Stop","Futility?"))

T1Reps[p,] <- T1

} ## End of for loop
## Visualize the pattern in Sig NotSig data
Output = aggregate(KM ~ Recist +KM.1  +ST.1  +Cox  ,data = T1Reps, FUN = length, drop = FALSE)
