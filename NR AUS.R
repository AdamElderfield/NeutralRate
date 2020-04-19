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

EST <- list()

for(i in c( 0.005, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 3.0)){
  

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

SIG <-sqrt(i)  #sqrt(0.4)

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
   R[12,12] <- p[13]+p[8]*p[4]
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

theta <- c(1.53,-0.54, -0.15, -0.62, -0.32, 0.39, sqrt(0.38), sqrt(0.54), sqrt(0.05), 4, 1.113391 , sqrt(0.15), sqrt(0.07), sqrt(0.79),0.5) # estimates from RBA  paper

#theta <- c(1.33,-0.54, -0.15, -0.62, -0.32, 0.39, sqrt(0.38), sqrt(0.54), sqrt(0.05), 4, 1.113391 , sqrt(0.15), sqrt(0.07), sqrt(0.79),0.5) # estimates from Sheen  paper

LC <- rep(-Inf,15)
UC <- rep(Inf,15)
UC[3] <- -0.0025
LC[c(7:9,11:15)] <- exp(-8)
UC[c(7:9,11:15)] <- exp(12)


NRDLM.est <-  dlmMLE(y = cbind(NRdata$log.output,NRdata$unr,NRdata$real.r,NRdata$Inflation, NRdata$Inflation.e), parm = theta, build = buildNRDLM, 
                     lower = LC, upper = UC,
                     control = list(trace = 6, REPORT = 5, maxit = 1000), debug = FALSE, method = "L-BFGS-B", hessian = TRUE)


# Create table of estimates, standard erros, and t-stats
EST[[paste0("sigma z = ",i)]] <-  tibble(Var =c("a1", "a2", "a3","g1","b1","b2", "sigma y*", "sigma ygap", "sigma mu","d", "sigma r","sigma ugap","sigma u*","sigma pi","sigma pie" ),
       Par =round(NRDLM.est$par,5),
       SE = sqrt(diag(solve(NRDLM.est$hessian))),
       `t-stat` = Par/SE) %>% 
  mutate(SE = ifelse(grepl("*sigma", .$Var), NA, SE),
         `t-stat`= ifelse(grepl("*sigma", .$Var), NA, Par/SE),
         Par = ifelse(grepl("*sigma", .$Var), Par, Par),
  )

}
#--------------------------------------------------------------------------------------------------------------------------
# filtered and smoothed estimates
#--------------------------------------------------------------------------------------------------------------------------

#NRDLMbuilt <- buildNRDLM(NRDLM.est$par)

NRDLMbuilt <- buildNRDLM(EST$`sigma z = 1.4`$Par)


filtered <- dlmFilter(y =cbind(NRdata$log.output,NRdata$unr,NRdata$real.r,NRdata$Inflation, NRdata$Inflation.e), mod = NRDLMbuilt)

smoothed <- dlmSmooth(y = cbind(NRdata$log.output,NRdata$unr,NRdata$real.r,NRdata$Inflation,NRdata$Inflation.e), mod = NRDLMbuilt)


cbind(filtered$y[-c(2:9),1],filtered$m[-c(1:9),1]) %>% matplot(type ="l")

cbind(filtered$y[-c(2:9),1],smoothed$s[-c(1:9),1]) %>% matplot(type ="l")

cbind(filtered$y[-c(2:9),3],smoothed$s[-c(1:9),7]) %>% matplot(type ="l")


#--------------------------------------------------------------------------------------------------------------------------
# Using the LW three stage approach
#--------------------------------------------------------------------------------------------------------------------------

  
  stage1NRDLM <- dlm(
    
    FF = matrix(0,4,12),
    
    V = diag(0.00001, 4),
    
    GG =diag(0,12),
    
    W = diag(0,12),
    
    m0 = rep(0,12),
    
    C0 = diag(10000000,12)
    
    
  )
  
  # Matrix to parametrise VCV matrix W
  R <- diag(0,12)
  

  # build DLM
  stage1buildNRDLM <- function(p){
    
    FF(stage1NRDLM)[1,1] <- 1
    FF(stage1NRDLM)[1,2] <- 1
    
    FF(stage1NRDLM)[2,6] <- 1
    FF(stage1NRDLM)[2,7] <- 1
    
    FF(stage1NRDLM)[3,9] <- 1
    
    FF(stage1NRDLM)[4,12] <- 1

    GG(stage1NRDLM)[1,1]  <- 1 
    GG(stage1NRDLM)[1,5]  <- 1
    GG(stage1NRDLM)[2,2]  <- p[1] 
    GG(stage1NRDLM)[2,3]  <- p[2]
    GG(stage1NRDLM)[3,2]  <- 1
    GG(stage1NRDLM)[4,3]  <- 1
  
    GG(stage1NRDLM)[5,5] <- 1
    GG(stage1NRDLM)[6,6] <- 1
    
    GG(stage1NRDLM)[7,2]  <- p[3]*(0.4*p[1]+0.3)        
    GG(stage1NRDLM)[7,3]  <- p[3]*(0.4*p[2]+0.2)        
    GG(stage1NRDLM)[7,4]  <- p[3]*0.1       
    GG(stage1NRDLM)[8,7] <- 1           
    
    GG(stage1NRDLM)[9,8] <- p[4]
    GG(stage1NRDLM)[9,9] <- p[5]/3
    GG(stage1NRDLM)[9,10] <- p[5]/3
    GG(stage1NRDLM)[9,11] <- p[5]/3
    GG(stage1NRDLM)[9,12] <- (1-p[5])
    GG(stage1NRDLM)[10,9] <- 1
    GG(stage1NRDLM)[11,10] <- 1
    
    GG(stage1NRDLM)[12,12] <- 1
    
    
    # Variance covariance - RR'
    
    R[1,1] <- p[6]
    R[2,2] <- p[7]
    R[6,6] <- p[8]
    R[7,7] <- p[9]+p[4]*p[7]
    R[7,2] <- p[4]*p[7]
    R[9,9] <- p[10]+(1-p[5])*p[11]
    R[12,9] <- (1-p[5])*p[11]
    R[12,12] <- p[11]
    
    W(stage1NRDLM) <- R%*%t(R)

    # PRIORS
    m0(stage1NRDLM) <- c(NRdata$log.output[1],0,0,0,mean(diff(NRdata$log.output[1:4])),NRdata$unr[1],0,0,NRdata$Inflation[3],NRdata$Inflation[2],NRdata$Inflation[1], NRdata$Inflation.e[1])
    
    

    return(stage1NRDLM)
    
  }
  
  theta <- c(1.53,-0.54, -0.15, -0.62, 0.32, sqrt(0.38), sqrt(0.54), sqrt(0.05), sqrt(0.15), sqrt(0.07), sqrt(0.79),0.5) # estimates from RBA  paper
  

  LC <- rep(-Inf,15)
  UC <- rep(Inf,15)
  # UC[3] <- -0.0025
  # LC[c(7:9,11:15)] <- exp(-8)
  # UC[c(7:9,11:15)] <- exp(12)
  
  
  stage1NRDLM.est <-  dlmMLE(y = cbind(NRdata$log.output,NRdata$unr,NRdata$Inflation, NRdata$Inflation.e), parm = theta, build = stage1buildNRDLM, 
                       lower = LC, upper = UC,
                       control = list(trace = 6, REPORT = 5, maxit = 1000), debug = FALSE, method = "L-BFGS-B", hessian = TRUE)
  
  
#--------------------------------------------------------------------------------------------------------------------------
# filtered and smoothed estimates
#--------------------------------------------------------------------------------------------------------------------------

stage1NRDLMbuilt <- stage1buildNRDLM(stage1NRDLM.est$par)


stage1filtered <- dlmFilter(y =cbind(NRdata$log.output,NRdata$unr,NRdata$Inflation, NRdata$Inflation.e), mod = stage1NRDLMbuilt)

stage1smoothed <- dlmSmooth(y = cbind(NRdata$log.output,NRdata$unr,NRdata$Inflation,NRdata$Inflation.e), mod = stage1NRDLMbuilt)


cbind(filtered$y[-c(2:9),1],filtered$m[-c(1:9),1]) %>% matplot(type ="l")

cbind(filtered$y[-c(2:9),1],smoothed$s[-c(1:9),1]) %>% matplot(type ="l")

cbind(stage1filtered$m[-c(2:9),5],stage1smoothed$s[-c(1:9),5]) %>% matplot(type ="l")


#--------------------------------------------------------------------------------------------------------------------------
# Median unbiased estimate of the variance for mu (ypot growth) 
#--------------------------------------------------------------------------------------------------------------------------

# Using potential output (smoothed)
# function does a rolling regression (chow break point test) on smoothed potential and a constant, searching for a break and returns lambda g (see paper for more details)
lambda.g <- median.unbiased.estimator.stage1(stage1smoothed$s[,1])



#--------------------------------------------------------------------------------------------------------------------------
# Using the LW three stage approach - stage 2
#--------------------------------------------------------------------------------------------------------------------------

# mu becomes time varying, we introduce inerest rates and a constant


stage2NRDLM <- dlm(
  
  FF = matrix(0,5,15),
  
  V = diag(0.00001, 5),
  
  GG =diag(0,15),
  
  W = diag(0,15),
  
  m0 = rep(0,15),
  
  C0 = diag(10000000,15)
  
  
)

# Matrix to parametrise VCV matrix W
R <- diag(0,15)


# build DLM
stage2buildNRDLM <- function(p){
  
  FF(stage2NRDLM)[1,1] <- 1
  FF(stage2NRDLM)[1,2] <- 1
  
  FF(stage2NRDLM)[2,9] <- 1
  FF(stage2NRDLM)[2,10] <- 1
  
  FF(stage2NRDLM)[3,12] <- 1
  
  FF(stage2NRDLM)[4,15] <- 1

  FF(stage2NRDLM)[5,6] <- 1
  
  
  GG(stage2NRDLM)[1,1]  <- 1 
  GG(stage2NRDLM)[1,5]  <- 1
  GG(stage2NRDLM)[2,2]  <- p[1] 
  GG(stage2NRDLM)[2,3]  <- p[2]
  GG(stage2NRDLM)[2,6]  <- p[3]/2
  GG(stage2NRDLM)[2,7]  <- p[3]/2
  GG(stage2NRDLM)[2,8]  <- 1
  GG(stage2NRDLM)[2,5]  <- p[14]
  GG(stage2NRDLM)[3,2]  <- 1
  GG(stage2NRDLM)[4,3]  <- 1
  
  GG(stage2NRDLM)[5,5] <- 1
  GG(stage2NRDLM)[6,6] <- 1
  GG(stage2NRDLM)[7,6] <- 1
  GG(stage2NRDLM)[8,8] <- 1
  
  GG(stage2NRDLM)[9,9] <- 1
  
  GG(stage2NRDLM)[10,2]  <- p[4]*(0.4*p[1]+0.3)        
  GG(stage2NRDLM)[10,3]  <- p[4]*(0.4*p[2]+0.2)        
  GG(stage2NRDLM)[10,4]  <- p[4]*0.1       
  GG(stage2NRDLM)[10,6]  <- p[4]*0.4*p[3]/2       
  GG(stage2NRDLM)[10,7]  <- p[4]*0.4*p[3]/2       
  GG(stage2NRDLM)[10,8] <- p[4]*0.4            
  GG(stage2NRDLM)[10,5] <- p[4]*0.4*p[14]            
  
  GG(stage2NRDLM)[11,10] <- 1           
  
  GG(stage2NRDLM)[12,11] <- p[5]
  GG(stage2NRDLM)[12,12] <- p[6]/3
  GG(stage2NRDLM)[12,13] <- p[6]/3
  GG(stage2NRDLM)[12,14] <- p[6]/3
  GG(stage2NRDLM)[12,15] <- (1-p[6])
  GG(stage2NRDLM)[13,12] <- 1
  GG(stage2NRDLM)[14,13] <- 1

  GG(stage2NRDLM)[15,15] <- 1
  
  
  # Variance covariance - RR'
  
  R[1,1] <- p[7]
  R[2,2] <- p[8]
  R[5,5] <-lambda.g*p[7]
  R[9,9] <- p[9]
  R[10,10] <- p[10]+p[4]*p[8]
  R[10,2] <- p[4]*p[8]
  R[12,12] <- p[11]+(1-p[6])*p[12]
  R[15,12] <- (1-p[6])*p[12]
  R[15,15] <- p[12]
  R[6,6] <- p[13]
  
  
  W(stage2NRDLM) <- R%*%t(R)
  
  # PRIORS
  m0(stage2NRDLM) <- c(NRdata$log.output[1],0,0,0,mean(diff(NRdata$log.output[1:4])),NRdata$real.r[2],NRdata$real.r[1],0,NRdata$unr[1],0,0,NRdata$Inflation[3],NRdata$Inflation[2],NRdata$Inflation[1], NRdata$Inflation.e[1])
  
  
  
  return(stage2NRDLM)
  
}

theta <- c(1.53,-0.54 ,-0.15, -0.62, -0.32, 0.5, sqrt(0.38), sqrt(0.54), sqrt(0.15), sqrt(0.07),sqrt(0.79),0.5,0.5, 0) # estimates from RBA  paper
#theta <- c(1.53,-0.54 ,-0.15, -0.62, -0.32, 0.5, 0.5, 0.5, 0.5, 0.5,0.5,0.5) #

LC <- rep(-Inf,14)
UC <- rep(Inf,14)
UC[3] <- -0.0025
UC[5] <- -0.25   # Parameter on Phillips curve is problematic
LC[c(7:12)] <- exp(-8)
UC[c(7:12)] <- exp(12)


stage2NRDLM.est <-  dlmMLE(y = cbind(NRdata$log.output,NRdata$unr,NRdata$Inflation, NRdata$Inflation.e, NRdata$real.r), parm = theta, build = stage2buildNRDLM, 
                           lower = LC, upper = UC,
                           control = list(trace = 6, REPORT = 5, maxit = 1000), debug = FALSE, method = "L-BFGS-B", hessian = TRUE)


#--------------------------------------------------------------------------------------------------------------------------
# filtered and smoothed estimates
#--------------------------------------------------------------------------------------------------------------------------

stage2NRDLMbuilt <- stage2buildNRDLM(stage2NRDLM.est$par)


stage2filtered <- dlmFilter(y =cbind(NRdata$log.output,NRdata$unr,NRdata$Inflation, NRdata$Inflation.e,NRdata$real.r), mod = stage2NRDLMbuilt)

stage2smoothed <- dlmSmooth(y = cbind(NRdata$log.output,NRdata$unr,NRdata$Inflation,NRdata$Inflation.e,NRdata$real.r), mod = stage2NRDLMbuilt)


cbind(stage2filtered$y[-c(2:9),1],stage2filtered$m[-c(1:9),1]) %>% matplot(type ="l")

cbind(stage2filtered$y[-c(2:9),1],stage2smoothed$s[-c(1:9),1]) %>% matplot(type ="l")

cbind(stage2filtered$m[-c(1:9),5],stage2smoothed$s[-c(1:9),5]) %>% matplot(type ="l")


#--------------------------------------------------------------------------------------------------------------------------
# Median unbiased estimate of the variance for mu (ypot growth) 
#--------------------------------------------------------------------------------------------------------------------------
output.gap.smoothed <- stage2smoothed$s[,2] 

## Inputs for median.unbiased.estimator.stage2.R
y <- output.gap.smoothed[3:length(output.gap.smoothed)]
x <- cbind(output.gap.smoothed[2:(length(output.gap.smoothed)-1)],
           output.gap.smoothed[1:(length(output.gap.smoothed)-2)],
           (NRdata$real.r[2:length(NRdata$real.r)]+NRdata$real.r.l1[2:length(NRdata$real.r)])/2,
           stage2smoothed$s[2:(length(output.gap.smoothed)-1),5],
           rep(1,length(NRdata$log.output)-1))

lambda.z <- median.unbiased.estimator.stage2(y,x)

#--------------------------------------------------------------------------------------------------------------------------
# Stage 3 
#---------------------------------------------------------------------------------------------------------------

stage3NRDLM <- dlm(
  
  FF = matrix(0,5,17),
  
  V = diag(0.00001, 5),
  
  GG =diag(0,17),
  
  W = diag(0,17),
  
  m0 = rep(0,17),
  
  C0 = diag(10000000,17)
  
  
)

# Matrix to parametrise VCV matrix W
R <- diag(0,17)


# build DLM
stage3buildNRDLM <- function(p){
  
  FF(stage3NRDLM)[1,1] <- 1
  FF(stage3NRDLM)[1,2] <- 1
  
  FF(stage3NRDLM)[2,11] <- 1
  FF(stage3NRDLM)[2,12] <- 1
  
  FF(stage3NRDLM)[3,14] <- 1
  
  FF(stage3NRDLM)[4,17] <- 1
  
  FF(stage3NRDLM)[5,7] <- 1
  
  
  GG(stage3NRDLM)[1,1]  <- 1 
  GG(stage3NRDLM)[1,5]  <- 1
  GG(stage3NRDLM)[2,2]  <- p[1] 
  GG(stage3NRDLM)[2,3]  <- p[2]

  GG(stage3NRDLM)[2,7]  <- p[3]/2
  GG(stage3NRDLM)[2,8]  <- p[3]/2
  GG(stage3NRDLM)[2,5]  <- -p[4]*p[3]/2
  GG(stage3NRDLM)[2,6]  <- -p[4]*p[3]/2
  
  GG(stage3NRDLM)[2,9]  <- -p[3]/2
  GG(stage3NRDLM)[2,10]  <- -p[3]/2
  
  
  GG(stage3NRDLM)[3,2]  <- 1 #ygap(-1)
  GG(stage3NRDLM)[4,3]  <- 1
  
  GG(stage3NRDLM)[5,5] <- 1  # mu
  GG(stage3NRDLM)[6,5] <- 1
  
  GG(stage3NRDLM)[7,7] <- 1   # r
  GG(stage3NRDLM)[8,7] <- 1
  
  GG(stage3NRDLM)[9,9] <- 1  # z
  GG(stage3NRDLM)[10,9] <- 1
  
  GG(stage3NRDLM)[11,11] <- 1  #u*
  
  GG(stage3NRDLM)[12,2]  <- p[5]*(0.4*p[1]+0.3)        
  GG(stage3NRDLM)[12,3]  <- p[5]*(0.4*p[2]+0.2)        
  GG(stage3NRDLM)[12,4]  <- p[5]*0.1       
  GG(stage3NRDLM)[12,5]  <- -p[5]*0.4*p[4]*p[3]/2       
  GG(stage3NRDLM)[12,6]  <- -p[5]*0.4*p[4]*p[3]/2       
  GG(stage3NRDLM)[12,7] <- p[5]*0.4*p[3]/2            
  GG(stage3NRDLM)[12,8] <- p[5]*0.4*p[3]/2  
  GG(stage3NRDLM)[12,9] <- -p[5]*0.4*p[3]/2            
  GG(stage3NRDLM)[12,10] <- -p[5]*0.4*p[3]/2  
  
  
  GG(stage3NRDLM)[13,12] <- 1           
  
  GG(stage3NRDLM)[14,13] <- p[6]
  GG(stage3NRDLM)[14,14] <- p[7]/3
  GG(stage3NRDLM)[14,15] <- p[7]/3
  GG(stage3NRDLM)[14,16] <- p[7]/3
  GG(stage3NRDLM)[14,17] <- (1-p[7])
  GG(stage3NRDLM)[15,14] <- 1
  GG(stage3NRDLM)[16,15] <- 1
  
  GG(stage3NRDLM)[17,17] <- 1
  
  
  # Variance covariance - RR'
  
  # R[1,1] <- (1+lambda.g)*p[8]
  # R[5,1] <- lambda.g*p[8]
  # R[2,2] <- p[9]
  # R[5,5] <-lambda.g*p[8]
  # R[9,9] <- (lambda.z*p[9])/p[3]
  # R[11,11] <- p[10]
  # R[12,12] <- p[11]+p[5]*p[9]
  # R[12,2] <- p[5]*p[9]
  # R[14,14] <- p[12]+(1-p[7])*p[13]
  # R[17,14] <- (1-p[7])*p[13]
  # R[17,17] <- p[13]
  # R[7,7] <- p[14]
  

   R[1,1] <- (1+lambda.g^2)*p[8]^2
   R[5,1] <- R[1,5] <- (lambda.g*p[8])^2
   R[2,2] <- p[9]^2
   R[5,5] <- (lambda.g*p[8])^2
   R[9,9] <- (lambda.z*p[9]/p[3])^2
   R[11,11] <- p[10]^2
   R[12,12] <- p[11]^2+(p[5]*p[9])^2
   R[12,2] <-  R[2,12] <-  (p[5]*p[9])^2
   R[14,14] <- p[12]^2+((1-p[7])*p[13])^2
   R[17,14] <- R[14,17] <-  ((1-p[7])*p[13])^2
   R[17,17] <- p[13]^2
   R[7,7] <- p[14]^2
  
    
  W(stage3NRDLM) <- R #R%*%t(R)
  
  # PRIORS
  m0(stage3NRDLM) <- c(NRdata$log.output[1],0,0,0,mean(diff(NRdata$log.output[1:4])),mean(diff(NRdata$log.output[1:5])),NRdata$real.r[2],NRdata$real.r[1],0,0,NRdata$unr[1],0,0,NRdata$Inflation[3],NRdata$Inflation[2],NRdata$Inflation[1], NRdata$Inflation.e[1])
  
  
  
  return(stage3NRDLM)
  
}

theta <- c(1.53,-0.54 ,-0.15, 1, -0.62, -0.32, 0.5, sqrt(0.38), sqrt(0.54), sqrt(0.15), sqrt(0.07),sqrt(0.79),0.5,0.5) # estimates from RBA  paper

LC <- rep(-Inf,14)
UC <- rep(Inf,14)
UC[3] <- -0.005  # Calibrated to give a neutral rate consistent with other studies 
UC[6] <- -0.32   # Parameter on Phillips curve is problematic
LC[c(7:14)] <- exp(-8)
UC[c(7:14)] <- exp(12)


stage3NRDLM.est <-  dlmMLE(y = cbind(NRdata$log.output,NRdata$unr,NRdata$Inflation, NRdata$Inflation.e, NRdata$real.r), parm = theta, build = stage3buildNRDLM, 
                           lower = LC, upper = UC,
                           control = list(trace = 6, REPORT = 5, maxit = 1000), debug = FALSE, method = "L-BFGS-B", hessian = TRUE)



#--------------------------------------------------------------------------------------------------------------------------
# filtered and smoothed estimates
#--------------------------------------------------------------------------------------------------------------------------

# Model is estimating a postive output gap / rnu gap. Most likely something going on with phillips curve. 

stage3NRDLMbuilt <- stage3buildNRDLM(stage3NRDLM.est$par)

stage3filtered <- dlmFilter(y =cbind(NRdata$log.output,NRdata$unr,NRdata$Inflation, NRdata$Inflation.e,NRdata$real.r), mod = stage3NRDLMbuilt)

stage3smoothed <- dlmSmooth(y = cbind(NRdata$log.output,NRdata$unr,NRdata$Inflation,NRdata$Inflation.e,NRdata$real.r), mod = stage3NRDLMbuilt)


## One-sided (filtered) R-star estimates
trend.filtered      <- stage3filtered$m[,4] * 4
z.filtered          <- stage3filtered$m[,9]
rstar.filtered      <- trend.filtered * stage3NRDLM.est$par[4] + z.filtered

cbind(NRdata$real.r,rstar.filtered[-1])

## One-sided (smoothed) R-star estimates
trend.smoothed     <- stage3smoothed$s[,5] * 4
z.smoothed          <- stage3smoothed$s[,9]
rstar.smoothed      <- trend.smoothed * stage3NRDLM.est$par[4] + z.smoothed

cbind(NRdata$real.r[-c(1:9)],rstar.smoothed[-c(1:10)]) %>% 
  matplot(type = "l")


cbind(stage3filtered$y[-c(2:9),1],stage3filtered$m[-c(1:9),1]) %>% matplot(type ="l")

cbind(stage3filtered$y[-c(2:9),1],stage3smoothed$s[-c(1:9),1]) %>% matplot(type ="l")

cbind(stage3filtered$m[-c(1:9),5],stage3smoothed$s[-c(1:9),5]) %>% matplot(type ="l")

#--------------------------------------------------------------------------------------------------------------------------
# Charts and tables
#--------------------------------------------------------------------------------------------------------------------------

# Chart one - potential output and actual output

# Chart two - output gap

# Chart three - nairu and unemployment rate

# Chart four - r* and real cash rate



