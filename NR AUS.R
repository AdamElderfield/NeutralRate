#-----------------------------------------------------------------------------------------
# Neutral rate estimation for Australia
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# load libraries
#-----------------------------------------------------------------------------------------

library(dplyr)
library(tidyverse)
library(lubridate)
library(zoo)
library(xts)
library(nloptr)
library(dlm)
source("~/GitHub/TST themes/Chart themes.R")
source("~/GitHub/Packages/tst.package/R/tst.macrodata.R")
source("~/GitHub/LW/median.unbiased.estimator.stage1.R")
source("~/GitHub/LW/median.unbiased.estimator.stage2.R")
source("~/GitHub/Neutral rate/HLW/hpfilter.R")


#-----------------------------------------------------------------------------------------
# Prepare data
#-----------------------------------------------------------------------------------------
# Load model database
m <- tst.macrodata()

# We need output, unemployment, real interest rates, and inflation, terms of trade and MTPG

NRdata <- tibble(log.output =    100*log(m$GDPE_r),
                 log.output.l1 = lag(log.output),
                 log.output.l2 = lag(log.output,2),
                 log.output.l3 = lag(log.output,3),
                 log.output.l4 = lag(log.output,4),
                 
                 
                 d_ltot = 100*(log(m$PXGS/m$PMGS)-log(lag(m$PXGS)/lag(m$PMGS))),
                 d_ltot.l1 = lag(d_ltot),
                 d_ltot.l2 = lag(d_ltot,2),
                 d_ltot.l3 = lag(d_ltot,3),
                 d_ltot.l4 = lag(d_ltot,4),
                 
                 unr = m$rate_UNE,
                 unr.l1 = lag(unr),
                 d_unr = m$d_rate_UNE,
                 
                 # Change to underlying
                 Inflation.mgs = m$d4l_PMGS,
                 Inflation =    m$d4l_CPI, 
                 Inflation.h = m$d4l_CPI,
                 Inflation.q = m$dl_CPI,
                 Inflation.l1 = lag(Inflation),
                 Inflation.q.l2 = lag(Inflation.q,2),
                 Inflation.q.l3 = lag(Inflation.q,3),
                 Inflation.q.l4 = lag(Inflation.q,4),
                 
                 # Pull NS INFE in
                 Inflation.q.e = rollapplyr(Inflation.q,list(-(4:1)),mean,fill=NA),
                 Inflation.e = rollapplyr(Inflation,list(-(12:1)),mean,fill=NA),
                 
#                 MTP_GAP =  as.numeric(neverhpfilter::yth_filter(xts(log(m$MTP_GDP),order.by = m$Date))$y.cycle),
                 
                 nominal.r = m$rate_90,
                 
                 Date = m$Date) %>% 
  mutate(real.r = nominal.r-Inflation.h,
         real.r.l1 = lag(real.r),
         real.r.l2 = lag(real.r,2),
         real.r.l8 = lag(real.r,8),
         Inflation.e = Inflation.e, #-Inflation.l1,
         Inflation.mgs = Inflation.mgs-Inflation.l1,
  ) %>% 
  filter(Date >= "1982-12-01")



#-----------------------------------------------------------------------------------------
# DLM1 - SEE LATEX NOTE 
#-----------------------------------------------------------------------------------------


NRDLM <- dlm(
  
  FF = matrix(0,5,17),
  
  V = diag(0.00001, 5),
  
  GG =diag(0,17),
  
  #JGG = diag(1,17),
  
  W = diag(0,17),
  
  m0 = rep(0,17),
  
  C0 = diag(10000000,17)
  
  #X = NRdata[,c("Inflation.e")]
  
)

# Matrix to parametrise VCV matrix W
R <- diag(0,17)

# Set all elements of JGG to zero (will change below)
#NRDLM$JGG <- diag(0,17)

SIG <- sqrt(2)

# build DLM
buildNRDLM <- function(p){
  
  FF(NRDLM)[1,1] <- 1
  FF(NRDLM)[1,2] <- 1
  
  FF(NRDLM)[2,11] <- 1
  FF(NRDLM)[2,12] <- 1
  
  FF(NRDLM)[3,9] <- 1

  FF(NRDLM)[4,14] <- 1
  FF(NRDLM)[5,17] <- 1
  
  GG(NRDLM)[1,1]  <- 1 
  GG(NRDLM)[1,5]  <- 1
  GG(NRDLM)[2,2]  <- p[1] 
  GG(NRDLM)[2,3]  <- p[2]
  GG(NRDLM)[2,7]  <- -p[3]/2
  GG(NRDLM)[2,8]  <- -p[3]/2
  GG(NRDLM)[2,9]  <- p[3]/2
  GG(NRDLM)[2,10]  <- p[3]/2
  GG(NRDLM)[3,2]  <- 1
  GG(NRDLM)[4,3]  <- 1
  
  
  GG(NRDLM)[5,5]  <- 1 
  GG(NRDLM)[6,6]  <- 1 
  GG(NRDLM)[7,5]  <- p[10]  #4 
  GG(NRDLM)[7,6]  <- 1 
  GG(NRDLM)[8,7]  <- 1 
  
  GG(NRDLM)[9,9]  <- 1 
  GG(NRDLM)[10,9]  <- 1 
  
  GG(NRDLM)[11,11] <- 1
  
  GG(NRDLM)[12,2]  <- p[4]*(0.4*p[1]+0.3)        
  GG(NRDLM)[12,3]  <- p[4]*(0.4*p[2]+0.2)        
  GG(NRDLM)[12,4]  <- p[4]*0.1       
  GG(NRDLM)[12,7]  <- p[4]*0.4*-p[3]/2       
  GG(NRDLM)[12,8]  <- p[4]*0.4*-p[3]/2       
  GG(NRDLM)[12,9] <- p[4]*0.4*p[3]/2            
  GG(NRDLM)[12,10] <- p[4]*0.4*p[3]/2           
  GG(NRDLM)[13,12] <- 1           
  
  GG(NRDLM)[14,13] <- p[5]
  GG(NRDLM)[14,14] <- p[6]/3
  GG(NRDLM)[14,15] <- p[6]/3
  GG(NRDLM)[14,16] <- p[6]/3
  GG(NRDLM)[14,17] <- (1-p[6])
  GG(NRDLM)[15,14] <- 1
  GG(NRDLM)[16,15] <- 1
  
  GG(NRDLM)[17,17] <- 1
  
  
  #JGG(NRDLM)[14,17] <- 1
  
  
  # Variance covariance - RR'
  
   R[1,1] <- p[7]
   R[2,2] <- p[8]
   R[5,5] <- p[9]
   R[6,6] <-  SIG   #p[10]
   R[7,6] <-  SIG   #p[10]
   R[7,5] <-  SIG*p[9]  #p[10]*p[9]
   R[7,7] <-  SIG*p[9]+p[11]  #p[10]*p[9]
   
   R[9,9] <- p[11]
   R[11,11] <- p[12]
   R[12,12] <- p[13]+p[8]
   R[12,2] <- p[4]*p[8]
   R[14,14] <- p[14]+p[15]
   R[17,14] <- p[6]*p[15]
   R[17,17] <- p[15]
    
   W(NRDLM) <- R%*%t(R)
  #################################################
  #########
  # USING MATRIx R to PICK OUT SHOCK VARIANCES LIKE IN ESTIMA EXAMPLE
  ########
  # R[1,1] <- 1
  # R[2,2] <- 1
  # R[5,5] <- 1
  # R[6,6] <- 1
  # R[7,6] <- 1
  # R[7,5] <- 4
  # R[9,9] <- 1
  # R[11,11] <- 1
  # R[12,12] <- 1
  # R[12,2] <- 1 #p[4]
  # R[14,14] <- 1
  # R[14,17] <- 1 #p[6]#(1-p[6]) 
  # R[17,17] <- 1
  #SHOCKS <-  matrix(c(p[7]^2,p[8]^2,0,0,p[9]^2,p[10]^2,0,0,p[11]^2,0,p[12]^2,p[13]^2,0,p[14]^2,0,0,p[15]^2), 17, 1, byrow = T)
  #W(NRDLM) <- diag(c(R%*%SHOCKS), 17)
  ###################################################
  
   # PRIORS
   m0(NRDLM) <- c(NRdata$log.output[1],0,0,0,mean(diff(NRdata$log.output[1:4])),0,NRdata$real.r[2],NRdata$real.r[1],NRdata$real.r[2],NRdata$real.r[1],NRdata$unr[1],0,0,NRdata$Inflation[3],NRdata$Inflation[2],NRdata$Inflation[1], NRdata$Inflation.e[1])
  
  
    # C0(NRDLM)[1,1] <- 100
    # #C0(NRDLM)[5,5] <- var(NRdata$log.output)
    # C0(NRDLM)[6,6] <- 100 #var(NRdata$real.r)
    # C0(NRDLM)[7,7] <- 100 #var(NRdata$real.r)
    # C0(NRDLM)[9,9] <- var(NRdata$real.r)
    # C0(NRDLM)[11,11] <-100 
    # C0(NRDLM)[12,12] <- var(NRdata$unr)
    # C0(NRDLM)[14,14] <- var(NRdata$Inflation)
    # C0(NRDLM)[14,14] <- var(NRdata$Inflation.e)
    # 
   
   
  return(NRDLM)
  
}

theta <- c(1.53,-0.54, -0.05, -0.62, -0.32, 0.39, sqrt(0.38), sqrt(0.54), sqrt(0.05), 4, 1.113391 , sqrt(0.15), sqrt(0.07), sqrt(0.79),0.5) # estimates from paper
# Estimate model

NRDLM.est <-  dlmMLE(y = cbind(NRdata$log.output,NRdata$unr,NRdata$real.r,NRdata$Inflation, NRdata$Inflation.e), parm = theta, build = buildNRDLM, 
                     control = list(trace = 6, REPORT = 5, maxit = 1000), debug = FALSE, method = "L-BFGS-B")

#lower =c(rep(-Inf,6),rep(exp(-8),3),-Inf,rep(exp(-8),5)), upper= c(rep(Inf,2),-0.00025,rep(Inf,3),rep(exp(12),3),Inf,rep(exp(12),5)),


#--------------------------------------------------------------------------------------------------------------------------
# filtered and smoothed estimates
#--------------------------------------------------------------------------------------------------------------------------

NRDLMbuilt <- buildNRDLM(NRDLM.est$par)

#NRDLMbuilt <- buildNRDLM(theta)


filtered <- dlmFilter(y =cbind(NRdata$log.output,NRdata$unr,NRdata$real.r,NRdata$Inflation, NRdata$Inflation.e), mod = NRDLMbuilt)

smoothed <- dlmSmooth(y = cbind(NRdata$log.output,NRdata$unr,NRdata$real.r,NRdata$Inflation), mod = NRDLMbuilt)


cbind(filtered$y[-c(2:9),1],filtered$m[-c(1:9),1]) %>% matplot(type ="l")

