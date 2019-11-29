#-----------------------------------------------------------------------------------------
# Neutral rate estimation for Australia, following DR
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# load libraries
#-----------------------------------------------------------------------------------------

library(dplyr)
library(tidyverse)
library(lubridate)
library(zoo)
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

# We need output, unemployment, real interest rates, and inflation

NRdata <- tibble(log.output =    100*log(m$GDPE_r),
                 log.output.l1 = lag(log.output),
                 log.output.l2 = lag(log.output,2),
                 log.output.l3 = lag(log.output,3),
                 log.output.l4 = lag(log.output,4),
                 
                 
                 unr = m$rate_UNE,
                 unr.l1 = lag(unr),
                 d_unr = m$d_rate_UNE,
                 
                 
                 Inflation =    100*m$d4l_CPI, #m$rate_INFNFNE,
                 Inflation.h = 100*m$d4l_CPI,
                 Inflation.q = 100*m$dl_CPI,
                 Inflation.l1 = lag(Inflation),
                 Inflation.q.l2 = lag(Inflation.q,2),
                 Inflation.q.l3 = lag(Inflation.q,3),
                 Inflation.q.l4 = lag(Inflation.q,4),
                 
                 Inflation.q.e = rollapplyr(Inflation.q,list(-(4:1)),mean,fill=NA),
                 Inflation.e = rollapplyr(Inflation,list(-(12:1)),mean,fill=NA),
                 
                 nominal.r = m$rate_90,
                 
                 Date = m$Date) %>% 
  filter(Date >= "1980-03-01") %>% 
  mutate(real.r = nominal.r-Inflation.h,
         real.r.l1 = lag(real.r),
         real.r.l2 = lag(real.r,2) 
  )


#-----------------------------------------------------------------------------------------
# DLM - stage1 no interest rate
#-----------------------------------------------------------------------------------------

# Build model

NRDLM1 <- dlm(
  
  
  FF = matrix(c(1,0,0,0,1,0,0,0,0,0,0, #y
                0,0,0,0,0,0,1,1,0,0,0, #u
                0,0,0,0,0,0,0,0,0,1,0, #pi
                0,0,0,0,0,0,0,0,0,0,1 #pi
               ) 
              , 
              ncol = 11, byrow = TRUE),
  

  V = diag(0.00001, 4),
  
  GG = diag(0,11),
  
  W = diag(0,11),
  

  m0 = rep(0,11),
  
  C0 = diag(10e7,11)
)


R <- diag(1,11)

buildNRDLM1 <- function(p){
  
  
  # Y gap (and lags)
  GG(NRDLM1)[1,1]  <- p[1]
  GG(NRDLM1)[1,2]  <- p[2]

  GG(NRDLM1)[2,1]  <- 1
  GG(NRDLM1)[3,2]  <- 1
  GG(NRDLM1)[4,3]  <- 1
  
  # Y*
  GG(NRDLM1)[5,5]  <- 1
  GG(NRDLM1)[5,6]  <- 1
  GG(NRDLM1)[6,6]  <- 1
  
  # U*
  GG(NRDLM1)[7,7]  <- 1
  
  #U gap
  GG(NRDLM1)[8,1]  <- p[1]*p[3]
  GG(NRDLM1)[8,2]  <- p[2]*p[3]
  GG(NRDLM1)[9,8]  <- 1
  
  # PI
  GG(NRDLM1)[10,8] <-  p[4]
  GG(NRDLM1)[10,10] <-  p[5]
  GG(NRDLM1)[10,11] <- 1-p[5]
  GG(NRDLM1)[11,11] <- 1

  dd <- c(exp(p[6]),0,0,0,exp(p[7]),exp(p[8]),exp(p[9]),exp(p[10]),0,exp(p[11]),exp(p[12]))
  
  R[8,1] <- p[3]
  R[10,11] <- 1-p[5]
  
  W(NRDLM1) <- tcrossprod(R*rep(dd, each =11))
  
  
  
  m0(NRDLM1) <- c(0,0,0,0,NRdata$log.output[1],mean(diff(NRdata$log.output)),
                 NRdata$unr[1],0,0,NRdata$Inflation[1],NRdata$Inflation.e[1])
  
  return(NRDLM1)
  
}

theta <- c(1.5,-0.5,-0.3,0.3,0.7,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5) 

# Estimate model
NRDLM1.est <-  dlmMLE(y = cbind(NRdata$log.output,NRdata$unr,NRdata$Inflation,NRdata$Inflation.e), parm = theta, build = buildNRDLM1, lower =c(rep(-Inf,5),rep(-8,7),-Inf), upper= c(rep(Inf,5),rep(12,7),Inf), 
                      control = list(trace = 1, REPORT = 5, maxit = 1000), hessian = TRUE)

# Create table of estimates, standard erros, and t-stats
tibble(Var =c("a1", "a2", "O.1","b1","b2","sigma ygap", "sigma y*", "sigma g", "sigma u*","sigma ugap", "sigma pi", "sigma piexp"),
       Par =round(NRDLM1.est$par,3),
       SE = sqrt(diag(solve(NRDLM1.est$hessian))),
       `t-stat` = Par/SE) %>% 
  mutate(SE = ifelse(grepl("*sigma", .$Var), NA, SE),
         `t-stat`= ifelse(grepl("*sigma", .$Var), NA, Par/SE),
         Par = ifelse(grepl("*sigma", .$Var), exp(Par), Par),
  )


#--------------------------------------------------------------------------------------------------------------------------
# filtered and smoothed estimates
#--------------------------------------------------------------------------------------------------------------------------

NRDLMbuilt1 <- buildNRDLM1(NRDLM1.est$par)

filtered1 <- dlmFilter(y =cbind(NRdata$log.output,NRdata$unr,NRdata$Inflation,NRdata$Inflation.e), mod = NRDLMbuilt1)

smoothed1 <- dlmSmooth(y = cbind(NRdata$log.output,NRdata$unr,NRdata$Inflation,NRdata$Inflation.e), mod = NRDLMbuilt1)

variances <-  dlmSvd2var(smoothed1$U.S,smoothed1$D.S) %>%
  lapply(sqrt) %>%
  sapply(diag) %>%
  t() %>% 
  data.frame()


# Chart one output and potential
tibble(Potential = smoothed1$s[-1,5],
       Output = NRdata$log.output,
       Date = NRdata$Date) %>% 
  ggplot()+
  geom_line(aes(Date,Output), colour = tst_colors[1])+
  geom_line(aes(Date,Potential), colour = tst_colors[2])+
  tst_theme()+
  scale_colour_tst()+
  xlab("")+
  ylab("100 times natural log")+
  labs(colour = "")+
  theme(legend.position = c(0.8,0.5))+
  # Label for potential output
  annotate("text", x=ymd("1990-06-01") , y=1230 , label = "Potential", colour = tst_colors[2])+
  # Label for actual output
  annotate("text", x=ymd("2011-06-01") , y=1270 , label = "Output", colour = tst_colors[1])+
  # Chart title
  ggtitle("Output and potential", subtitle = "Stage 1 - without interest rates")


# Chart two output gap
tibble(`Output gap` = smoothed1$s[-1,1],
       Date = NRdata$Date,
       upper = smoothed1$s[-1,1]+(variances[-1,1]*qnorm(0.025,lower.tail = FALSE)) ,
       lower = smoothed1$s[-1,1]-(variances[-1,1]*qnorm(0.025,lower.tail = FALSE))
) %>% 
  ggplot()+
  geom_line(aes(Date, `Output gap`),  colour = tst_colors[2])+
  geom_ribbon(aes(Date, ymin =upper , ymax = lower), alpha = 0.2)+
  tst_theme()+
  scale_colour_tst()+
  xlab("")+
  ylab("Per cent")+
  # Label for probability intervals
  annotate("text", x=ymd("2010-06-01") , y= -3.5, label = "95% probability limit", colour = "dark grey")+
  # Label for output gap
  annotate("text", x=ymd("2015-06-01") , y=2.5 , label = "Output gap", colour = tst_colors[2])+
  # Chart title
  ggtitle("The output gap", subtitle = "Stage 1 - without interest rates")

  


# Chart three NAIRU
tibble(NAIRU = smoothed1$s[-1,7],
       `Unemployment rate` = NRdata$unr,
       Date = NRdata$Date,
       upper = smoothed1$s[-1,7]+(variances[-1,7]*qnorm(0.025,lower.tail = FALSE)) ,
       lower = smoothed1$s[-1,7]-(variances[-1,7]*qnorm(0.025,lower.tail = FALSE))
       ) %>% 
  ggplot()+
  geom_line(aes(Date, NAIRU),  colour = tst_colors[2])+
  geom_line(aes(Date, `Unemployment rate`), colour = tst_colors[1])+
  geom_ribbon(aes(Date, ymin =upper , ymax = lower), alpha = 0.2)+
  tst_theme()+
  scale_colour_tst()+
  xlab("")+
  ylab("Per cent")+
  annotate("text", x=ymd("2001-06-01") , y= 4.3, label = "95% probability limit", colour = "dark grey")+

  annotate("text", x=ymd("2012-06-01") , y=7 , label = "NAIRU", colour = tst_colors[2])+
  
  annotate("text", x=ymd("1990-06-01") , y=5.4 , label = "Unemployment rate", colour = tst_colors[1])+
  
  ggtitle("Unemployment rate and the NAIRU", subtitle = "Stage 1 - without interest rates")

# For stage 
lambda.g <- median.unbiased.estimator.stage1(smoothed1$s[,5])


#-----------------------------------------------------------------------------------------
# DLM - stage 2 (freely estimated variances)
#-----------------------------------------------------------------------------------------

NRDLM <- dlm(
  
              # 1 2 3 4 5 6 7 8 9 0 1 2  
  FF = matrix(c(1,0,0,1,0,0,0,0,0,0, # y
                0,0,0,0,0,0,0,0,0,0, # u
                0,0,0,0,0,0,0,0,0,0 # pi
                )
              , 
              ncol = 12, byrow = TRUE),
  
  
  V = diag(0.00001, 6),
  
  GG =matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0, # ygap      1
               0,0,0,0,0,0,0,0,0,0,0,0,0, # ygap-1    2
               0,0,0,0,0,0,0,0,0,0,0,0,0, # ygap-2    3
               0,0,0,0,0,0,0,0,0,0,0,0,0, # y*        4
               0,0,0,0,0,0,0,0,0,0,0,0,0, # g         5
               0,0,0,0,0,0,0,0,0,0,0,0,0, # a         6
               0,0,0,0,0,0,0,0,0,0,0,0,0, # u*        7
               0,0,0,0,0,0,0,0,0,0,0,0,0, # ugap      8
               0,0,0,0,0,0,0,0,0,0,0,0,0, # ugap-1    9
               0,0,0,0,0,0,0,0,0,0,0,0,0 # Inflation 10  
                 ),ncol = 13, byrow = TRUE),
  
  
JGG = matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0, # ygap      1
               0,0,0,0,0,0,0,0,0,0,0,0,0, # ygap-1    2
               0,0,0,0,0,0,0,0,0,0,0,0,0, # ygap-2    3
               0,0,0,0,0,0,0,0,0,0,0,0,0, # y*        4
               0,0,0,0,0,0,0,0,0,0,0,0,0, # g         5
               0,0,0,0,0,0,0,0,0,0,0,0,0, # a         6
               0,0,0,0,0,0,0,0,0,0,0,0,0, # u*        7
               0,0,0,0,0,0,0,0,0,0,0,0,0, # ugap      8
               0,0,0,0,0,0,0,0,0,0,0,0,0, # ugap-1    9
               0,0,0,0,0,0,0,0,0,0,0,0,0 # Inflation 10  
                                ),ncol = 13, byrow = TRUE),
  
  W = matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0, # ygap      1
               0,0,0,0,0,0,0,0,0,0,0,0,0, # ygap-1    2
               0,0,0,0,0,0,0,0,0,0,0,0,0, # ygap-2    3
               0,0,0,0,0,0,0,0,0,0,0,0,0, # y*        4
               0,0,0,0,0,0,0,0,0,0,0,0,0, # g         5
               0,0,0,0,0,0,0,0,0,0,0,0,0, # a         6
               0,0,0,0,0,0,0,0,0,0,0,0,0, # r-1       7
               0,0,0,0,0,0,0,0,0,0,0,0,0, # r-2       8
               0,0,0,0,0,0,0,0,0,0,0,0,0, # u*        9
               0,0,0,0,0,0,0,0,0,0,0,0,0, # ugap      10
               0,0,0,0,0,0,0,0,0,0,0,0,0, # ugap-1    11
               0,0,0,0,0,0,0,0,0,0,0,0,0, # Inflation 12  
               0,0,0,0,0,0,0,0,0,0,0,0,0  # Inlfation expectations
  ),ncol = 13, byrow = TRUE),
  

  m0 = rep(0,13),
  
  C0 = diag(10e7,13)
)


R <- diag(13)

buildNRDLM <- function(p){
  
  
  # Y gap (and lags)
  GG(NRDLM)[1,1]  <- p[1]
  GG(NRDLM)[1,2]  <- p[2]
  GG(NRDLM)[1,7]  <- p[3]/2
  GG(NRDLM)[1,8]  <- p[3]/2
  GG(NRDLM)[1,6]  <- p[16]
  GG(NRDLM)[1,5]  <- p[4]
  
  GG(NRDLM)[2,1]  <- 1
  GG(NRDLM)[3,2]  <- 1

  # Y*
  GG(NRDLM)[4,4]  <- 1
  GG(NRDLM)[4,5]  <- 1
  GG(NRDLM)[5,5]  <- 1 
  
  # a0
  GG(NRDLM)[6,6]  <- 1

  # r and lags
  GG(NRDLM)[7,7]  <- 1
  GG(NRDLM)[8,7]  <- 1

  # U*
  GG(NRDLM)[9,9]  <- 1
  
  #U gap
  GG(NRDLM)[10,1]  <- p[1]*p[5]
  GG(NRDLM)[10,2]  <- p[2]*p[5]
  GG(NRDLM)[11,10]  <- 1
  
  #pi 
  GG(NRDLM)[12,11]  <- p[6]
  GG(NRDLM)[12,12]  <- p[7]
  GG(NRDLM)[12,13]  <- p[8]
  
  #PIEXP
  GG(NRDLM)[13,13]  <- 1
 
  # Variance covariance - RR'
 
  R[1,10] <- p[5]
  
  R[12,13] <- p[8]
  
   
  
  dd <- c(exp(p[9]),0,0,exp(p[10]),exp(log(lambda.g)+p[10]),0,exp(p[11]),0,exp(p[12]),exp(p[13]),0,exp(p[14]),exp(p[15]))
  
  W(NRDLM) <- tcrossprod(R*rep(dd, each =13))
  
  
  
  m0(NRDLM) <- c(0,0,0,NRdata$log.output[1],mean(diff(NRdata$log.output)),1,
                 NRdata$real.r[2],NRdata$real.r[1],NRdata$unr[1],0,0,NRdata$Inflation.l1[1],NRdata$Inflation.e[1])
  
  
  return(NRDLM)
  
}

theta <- c(1.5,-0.5,-0.8,0.5, -0.5,-0.5,0.5,0.3,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,0.5,-0.5) 

# Estimate model

NRDLM.est <-  dlmMLE(y = cbind(NRdata$log.output,NRdata$unr,NRdata$Inflation,NRdata$real.r,NRdata$Inflation.e,1), parm = theta, build = buildNRDLM, lower =c(rep(-Inf,6),rep(-8,9), -Inf,-8), upper= c(rep(Inf,6),rep(12,9),Inf,12),
                                          control = list(trace = 1, REPORT = 5, maxit = 1000), hessian = TRUE)


# Create table of estimates, standard erros, and t-stats
tibble(Var =c("a1", "a2", "a3", "a4", "O.1","b1","b2", "b3","sigma ygap", "sigma y*","sigma r", "sigma u*","sigma ugap", "sigma pi", "sigma piexp","a0","sigma g"),
       Par =round(NRDLM.est$par,3),
       SE = sqrt(diag(solve(NRDLM.est$hessian))),
       `t-stat` = Par/SE) %>% 
  mutate(SE = ifelse(grepl("*sigma", .$Var), NA, SE),
         `t-stat`= ifelse(grepl("*sigma", .$Var), NA, Par/SE),
         Par = ifelse(grepl("*sigma", .$Var), exp(Par), Par),
  )




#--------------------------------------------------------------------------------------------------------------------------
# filtered and smoothed estimates
#--------------------------------------------------------------------------------------------------------------------------

NRDLMbuilt <- buildNRDLM(NRDLM.est$par)

filtered <- dlmFilter(y =cbind(NRdata$log.output,NRdata$unr,NRdata$Inflation,NRdata$real.r,NRdata$Inflation.e,1), mod = NRDLMbuilt)

smoothed <- dlmSmooth(y = cbind(NRdata$log.output,NRdata$unr,NRdata$Inflation,NRdata$real.r,NRdata$Inflation.e,1), mod = NRDLMbuilt)


cbind(filtered$y[-c(1:9),1],filtered$m[-c(1:10),4]) %>% matplot(type ="l")


# Two-sided (smoothed) estimates
trend.smoothed      <- smoothed$s[,5] * 4
potential.smoothed  <- smoothed$s[,4]
output.gap.smoothed <- smoothed$s[,1]

# Inputs for median.unbiased.estimator.stage2.R
y <- output.gap.smoothed[3:length(output.gap.smoothed)]
x <- cbind(output.gap.smoothed[2:(length(output.gap.smoothed)-1)],
           output.gap.smoothed[1:(length(output.gap.smoothed)-2)],
           (x.data[,8]+x.data[,9])/2,
           NRsmoothed2$xi.s[,4],
           rep(1,start))


lambda.z <- median.unbiased.estimator.stage2(y ,x)


#-----------------------------------------------------------------------------------------
# DLM - stage 3 (freely estimated variances)
#-----------------------------------------------------------------------------------------

NRDLM <- dlm(
  
              # 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6
  FF = matrix(c(1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0, # y
                0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0, # u
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0, # pi
                0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0, # r
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1) # piexp
              , 
              ncol = 16, byrow = TRUE),
  
  
  V = diag(0.00001, 5),
  
  GG =matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # ygap      1
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # ygap-1    2
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # ygap-2    3
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # y*        4
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # g         5
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # r*        6
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # z         7
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # r gap     8
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # r gap-1   9
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # r-gap-2   10
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # r-gap-3   11
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # u*        12
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # ugap      13
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # ugap-1    14
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # Inflation 15  
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0  # Inlfation expectations
  ),ncol = 16, byrow = TRUE),
  
  
  W = matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # ygap-1      1
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # ygap-2    2
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # ygap-3    3
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # y*        4
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # g         5
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # r*-1        6
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # z         7
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # r*-2      8
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # r*-3      9
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # r-1       10
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # r-2       11
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # u*        12
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # ugap      13
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # ugap-1    14
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # Inflation 15  
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0  # Inlfation expectations
  ),ncol = 16, byrow = TRUE),
  
  
  m0 = rep(0,16),
  
  C0 = diag(100000,16)
)


R <- diag(0,16)

buildNRDLM <- function(p){
  
  #C0(NRDLM)[7,7] <- 1
  # Y gap (and lags)
  GG(NRDLM)[1,1]  <- p[1]
  GG(NRDLM)[1,2]  <- p[2]
  GG(NRDLM)[1,8]  <- -p[3]
  GG(NRDLM)[1,10]  <- p[3]
  
  GG(NRDLM)[2,1]  <- 1
  GG(NRDLM)[3,2]  <- 1
  
  # Y*
  GG(NRDLM)[4,4]  <- 1
  GG(NRDLM)[4,5]  <- 1
  GG(NRDLM)[5,5]  <- 1
  
  # r*
  GG(NRDLM)[6,5]  <- p[8]
  GG(NRDLM)[6,7]  <- 1#p[17]
  
  # z
  GG(NRDLM)[7,7]  <-  1#p[9]
  
  # r star lags
  GG(NRDLM)[8,6]  <- 1
  GG(NRDLM)[9,8]  <- 1
  
  # r and lags
  GG(NRDLM)[10,10]  <- 1
  GG(NRDLM)[11,10]  <- 1
  
  # U*
  GG(NRDLM)[12,12]  <- 1
  
  #U gap
  GG(NRDLM)[13,1]  <- p[1]*p[4]
  GG(NRDLM)[13,2]  <- p[2]*p[4]
  GG(NRDLM)[13,8]  <- -p[3]*p[4] 
  GG(NRDLM)[13,10] <- p[3]*p[4]
  GG(NRDLM)[14,13]  <- 1
  
  #pi 
  GG(NRDLM)[15,13]  <- p[5]
  GG(NRDLM)[15,15]  <- p[6]
  GG(NRDLM)[15,16]  <- -p[6]
  
  #PIEXP
  GG(NRDLM)[16,16]  <- 1

  # Variance covariance - RR'
  R[1,1] <- exp(p[10])
  R[4,4] <- exp(p[11])
  R[5,5] <- exp(p[12])
  R[7,7] <- sqrt(0.6) #exp(p[13])
  R[10,10] <- exp(p[14])
  R[12,12] <- exp(p[15])
  R[13,13] <- exp(p[16])
  R[15,15] <- exp(p[17])
  R[16,16] <- exp(p[18])
  
  R[7,6] <- sqrt(0.6)             # r* - sum of z and 4* g
  R[6,5] <- p[8]*exp(p[12])  
  
  R[13,1] <- p[4]*exp(p[10])
  R[16,15] <- p[7]*exp(p[17])
  R[10,1] <- p[3]*exp(p[14])
#  R[13,10] <- p[3]*exp(p[14])
# R[13,8] <- -p[3]*exp(p[14])
#  R[8,1] <- -p[3]*exp(p[14])
  
  
  W(NRDLM) <- R%*%t(R) 
  
  
  
  m0(NRDLM) <- c(0,
                 0,
                 0,
                 NRdata$log.output[1],
                 mean(diff(NRdata$log.output)),
                 NRdata$real.r[3],
                 0,
                 NRdata$real.r[2],
                 NRdata$real.r[1],
                 NRdata$real.r[2],
                 NRdata$real.r[1],
                 NRdata$unr[1],
                 0,
                 0,
                 NRdata$Inflation.l1[1],
                 NRdata$Inflation.e[1])
  
  
  return(NRDLM)
  
}

theta <- c(1.5,-0.5,-1,-0.3, -1,0.7,0.3, 4,0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5) 

# Estimate model

NRDLM.est <-  dlmMLE(y = cbind(NRdata$log.output,NRdata$unr,NRdata$Inflation,NRdata$real.r,NRdata$Inflation.e), parm = theta, build = buildNRDLM, lower =c(rep(-Inf,9),rep(-8,9)), upper= c(rep(Inf,2),-0.03,Inf,-0.2,rep(Inf,4),rep(12,9)),
                     control = list(trace = 6, REPORT = 5, maxit = 500), hessian = TRUE)

delta <- diag(sqrt(diag(exp(NRDLM.est$par))%*%solve(NRDLM.est$hessian)%*%diag(exp(NRDLM.est$par))))

# Create table of estimates, standard erros, and t-stats
tibble(Var =c("a1", "a2", "a3", "O.1","b1","b2","b3","d","z1","sigma ygap", "sigma y*", "sigma g","sigma r", "sigma u*","sigma ugap", "sigma pi", "sigma piexp","sigmaz"),
       Par =round(NRDLM.est$par,3),
       SE = sqrt(diag(solve(NRDLM.est$hessian))),
       SE2 = delta,
       `t-stat` = Par/SE) %>% 
  mutate(SE = ifelse(grepl("*sigma", .$Var), delta, SE),
         `z-stat`= ifelse(grepl("*sigma", .$Var), NA, Par/SE),
         Par = ifelse(grepl("*sigma", .$Var), exp(Par), Par),
  )


#--------------------------------------------------------------------------------------------------------------------------
# filtered and smoothed estimates
#--------------------------------------------------------------------------------------------------------------------------

NRDLMbuilt <- buildNRDLM(NRDLM.est$par)

filtered <- dlmFilter(y =cbind(NRdata$log.output,NRdata$unr,NRdata$Inflation,NRdata$real.r,NRdata$Inflation.e), mod = NRDLMbuilt)

smoothed <- dlmSmooth(y = cbind(NRdata$log.output,NRdata$unr,NRdata$Inflation,NRdata$real.r,NRdata$Inflation.e), mod = NRDLMbuilt)


cbind(filtered$y[-c(1:9),4],filtered$m[-c(1:10),6]) %>% matplot(type ="l")




#--------------------------------------------------------------------------------------------------------------------------
# try making r exogenous
#--------------------------------------------------------------------------------------------------------------------------


NRDLM <- dlm(
  
                # 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6
  FF = matrix(c(1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0, # y
                0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0, # u
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0, # pi
                0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0, # r
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1) # piexp
              , 
              ncol = 16, byrow = TRUE),
  
  
  V = diag(0.00001, 5),
  
  GG =matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # ygap      1
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # ygap-1    2
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # ygap-2    3
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # y*        4
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # g         5
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # r*        6
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # z         7
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # r gap     8
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # r gap-1   9
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # r-gap-2   10
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # r-gap-3   11
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # u*        12
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # ugap      13
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # ugap-1    14
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # Inflation 15  
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0  # Inlfation expectations
  ),ncol = 16, byrow = TRUE),
  
  
  W = matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # ygap      1
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # ygap-1    2
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # ygap-2    3
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # y*        4
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # g         5
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # r*        6
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # z         7
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # r*-1      8
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # r*-2      9
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # r-1       10
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # r-2       11
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # u*        12
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # ugap      13
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # ugap-1    14
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # Inflation 15  
               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0  # Inlfation expectations
  ),ncol = 16, byrow = TRUE),
  
  
  m0 = rep(0,16),
  
  C0 = diag(100000,16)
)


R <- diag(0,16)

buildNRDLM <- function(p){
  
  #C0(NRDLM)[7,7] <- 1
  # Y gap (and lags)
  GG(NRDLM)[1,1]  <- p[1]
  GG(NRDLM)[1,2]  <- p[2]
  GG(NRDLM)[1,9]  <- p[3]
  
  GG(NRDLM)[2,1]  <- 1
  GG(NRDLM)[3,2]  <- 1
  
  # Y*
  GG(NRDLM)[4,4]  <- 1
  GG(NRDLM)[4,5]  <- 1
  GG(NRDLM)[5,5]  <- 1
  
  # r*
  GG(NRDLM)[6,5]  <- p[8]
  GG(NRDLM)[6,7]  <- p[17]
  
  # z
  GG(NRDLM)[7,7]  <-  1 #p[9]
  
  # r star lags
  GG(NRDLM)[8,6]  <- 1
  GG(NRDLM)[9,9]  <- 0
  GG(NRDLM)[9,10]  <- 1
  GG(NRDLM)[9,8]  <- -1
  
  # r and lags
  GG(NRDLM)[10,10]  <- 1
  GG(NRDLM)[11,10]  <- 1
  
  # U*
  GG(NRDLM)[12,12]  <- 1
  
  #U gap
  GG(NRDLM)[13,1]  <- p[1]*p[4]
  GG(NRDLM)[13,2]  <- p[2]*p[4]
  GG(NRDLM)[13,9]  <- (p[3])*p[4] 
  #GG(NRDLM)[13,11] <- (p[3])*p[4]
  GG(NRDLM)[13,2]  <- p[4]
  GG(NRDLM)[14,13]  <- 1
  
  #pi 
  GG(NRDLM)[15,14]  <- p[5]
  GG(NRDLM)[15,15]  <- p[6]
  GG(NRDLM)[15,16]  <- p[7]
  
  #PIEXP
  GG(NRDLM)[16,16]  <- 1
  
  # Variance covariance - RR'
  R[1,1] <- exp(p[10])
  R[4,4] <- exp(p[11])
  R[5,5] <- exp(p[12])
  R[7,7] <-  sqrt(0.6)#exp(p[13])
  R[10,10] <- exp(p[14])
  R[12,12] <- exp(p[15])
  R[13,13] <- exp(p[16])
  R[15,15] <- exp(p[17])
  R[16,16] <- exp(p[18])
  
  #  R[7,6] <- 1             # r* - sum of z and 4* g
  #  R[6,5] <- p[8]  
  
  #  R[13,1] <- p[4]
  #  R[15,16] <- p[7]
  #  R[13,10] <- p[3]
  #  R[10,1] <- p[3]
  
  
  W(NRDLM) <- R%*%t(R) 
  
  
  
  m0(NRDLM) <- c(0,
                 0,
                 0,
                 NRdata$log.output[1],
                 mean(diff(NRdata$log.output)),
                 NRdata$real.r[3],
                 0,
                 NRdata$real.r[2],
                 0,
                 NRdata$real.r[2],
                 NRdata$real.r[1],
                 NRdata$unr[1],
                 0,
                 0,
                 NRdata$Inflation.l1[1],
                 NRdata$Inflation.e[1])
  
  
  return(NRDLM)
  
}

theta <- c(1.5,-0.5,-0.8,-0.3, 0.3,0.7,0.3, 4,0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5) 

# Estimate model

NRDLM.est <-  dlmMLE(y = cbind(NRdata$log.output,NRdata$unr,NRdata$Inflation,NRdata$real.r,NRdata$Inflation.e), parm = theta, build = buildNRDLM, lower =c(rep(-Inf,9),rep(-8,9)), upper= c(rep(Inf,2),-0.03,rep(Inf,6),rep(12,9)),
                     control = list(trace = 6, REPORT = 5, maxit = 500), hessian = TRUE)

delta <- diag(sqrt(diag(exp(NRDLM.est$par))%*%solve(NRDLM.est$hessian)%*%diag(exp(NRDLM.est$par))))

# Create table of estimates, standard erros, and t-stats
tibble(Var =c("a1", "a2", "a3", "O.1","b1","b2","b3","d","z1","sigma ygap", "sigma y*", "sigma g","sigma r", "sigma u*","sigma ugap", "sigma pi", "sigma piexp","sigmaz"),
       Par =round(NRDLM.est$par,3),
       SE = sqrt(diag(solve(NRDLM.est$hessian))),
       SE2 = delta,
       `t-stat` = Par/SE) %>% 
  mutate(SE = ifelse(grepl("*sigma", .$Var), delta, SE),
         `z-stat`= ifelse(grepl("*sigma", .$Var), NA, Par/SE),
         Par = ifelse(grepl("*sigma", .$Var), exp(Par), Par),
  )


#--------------------------------------------------------------------------------------------------------------------------
# filtered and smoothed estimates
#--------------------------------------------------------------------------------------------------------------------------

NRDLMbuilt <- buildNRDLM(NRDLM.est$par)

filtered <- dlmFilter(y =cbind(NRdata$log.output,NRdata$unr,NRdata$Inflation,NRdata$real.r,NRdata$Inflation.e), mod = NRDLMbuilt)

smoothed <- dlmSmooth(y = cbind(NRdata$log.output,NRdata$unr,NRdata$Inflation,NRdata$real.r,NRdata$Inflation.e), mod = NRDLMbuilt)


cbind(filtered$y[-c(1:9),4],filtered$m[-c(1:10),6]) %>% matplot(type ="l")



#--------------------------------------------------------------------------------------------------------------------------
# try making r exogenous 2
#--------------------------------------------------------------------------------------------------------------------------


NRDLM <- dlm(
  
              # 1 2 3 4 5 6 7 8 9 0 
  FF = matrix(c(1,0,0,0,0,0,0,0,0,0,0, # y
                0,0,0,0,0,0,0,0,0,0,0, # u
                0,0,0,0,0,0,0,0,0,0,0 # pi
                ) # piexp
              , 
              ncol = 16, byrow = TRUE),
  
  JFF = matrix(c(1,0,0,0,0,0,0,0,0,0,0, # y
                0,0,0,0,0,0,0,0,0,0,0, # u
                0,0,0,0,0,0,0,0,0,0,0# pi
  ) # piexp
  , 
  ncol = 16, byrow = TRUE),
  
  
  V = diag(0.00001, 5),
  
  GG =matrix(c(0,0,0,0,0,0,0,0,0,0,0, # a1      1
               0,0,0,0,0,0,0,0,0,0,0, # y*   2
               0,0,0,0,0,0,0,0,0,0,0, # g    3
               0,0,0,0,0,0,0,0,0,0,0, # g-1        4
               0,0,0,0,0,0,0,0,0,0,0, # g-2         5
               0,0,0,0,0,0,0,0,0,0,0, # y-1        6
               0,0,0,0,0,0,0,0,0,0,0, # y-2        7
               0,0,0,0,0,0,0,0,0,0,0, # z-1     8
               0,0,0,0,0,0,0,0,0,0,0, # z-2   9
               0,0,0,0,0,0,0,0,0,0,0,  # a2
               0,0,0,0,0,0,0,0,0,0,0  # u*   10
                 ),ncol = 11, byrow = TRUE),
  
  
  W = matrix(c(0,0,0,0,0,0,0,0,0,0, # a1      1
               0,0,0,0,0,0,0,0,0,0, # y*   2
               0,0,0,0,0,0,0,0,0,0, # g    3
               0,0,0,0,0,0,0,0,0,0, # g-1        4
               0,0,0,0,0,0,0,0,0,0, # g-2         5
               0,0,0,0,0,0,0,0,0,0, # y-1        6
               0,0,0,0,0,0,0,0,0,0, # y-2        7
               0,0,0,0,0,0,0,0,0,0, # z-1     8
               0,0,0,0,0,0,0,0,0,0, # z-2   9
               0,0,0,0,0,0,0,0,0,0 # u*   10
  ),ncol = 10, byrow = TRUE),
  
  
  m0 = rep(0,10),
  
  C0 = diag(100000,10)
)


R <- diag(0,10)

buildNRDLM <- function(p){
  
  #C0(NRDLM)[7,7] <- 1
  # Y gap (and lags)
  GG(NRDLM)[1,1]  <- p[1]
  GG(NRDLM)[1,2]  <- p[2]
  GG(NRDLM)[1,9]  <- p[3]
  
  GG(NRDLM)[2,1]  <- 1
  GG(NRDLM)[3,2]  <- 1
  
  # Y*
  GG(NRDLM)[4,4]  <- 1
  GG(NRDLM)[4,5]  <- 1
  GG(NRDLM)[5,5]  <- 1
  
  # r*
  GG(NRDLM)[6,5]  <- p[8]
  GG(NRDLM)[6,7]  <- p[17]
  
  # z
  GG(NRDLM)[7,7]  <-  1 #p[9]
  
  # r star lags
  GG(NRDLM)[8,6]  <- 1
  GG(NRDLM)[9,9]  <- 0
  GG(NRDLM)[9,10]  <- 1
  GG(NRDLM)[9,8]  <- -1
  
  # r and lags
  GG(NRDLM)[10,10]  <- 1
  GG(NRDLM)[11,10]  <- 1
  
  # U*
  GG(NRDLM)[12,12]  <- 1
  
  #U gap
  GG(NRDLM)[13,1]  <- p[1]*p[4]
  GG(NRDLM)[13,2]  <- p[2]*p[4]
  GG(NRDLM)[13,9]  <- (p[3])*p[4] 
  #GG(NRDLM)[13,11] <- (p[3])*p[4]
  GG(NRDLM)[13,2]  <- p[4]
  GG(NRDLM)[14,13]  <- 1
  
  #pi 
  GG(NRDLM)[15,14]  <- p[5]
  GG(NRDLM)[15,15]  <- p[6]
  GG(NRDLM)[15,16]  <- p[7]
  
  #PIEXP
  GG(NRDLM)[16,16]  <- 1
  
  # Variance covariance - RR'
  R[1,1] <- exp(p[10])
  R[4,4] <- exp(p[11])
  R[5,5] <- exp(p[12])
  R[7,7] <-  sqrt(0.6)#exp(p[13])
  R[10,10] <- exp(p[14])
  R[12,12] <- exp(p[15])
  R[13,13] <- exp(p[16])
  R[15,15] <- exp(p[17])
  R[16,16] <- exp(p[18])
  
  #  R[7,6] <- 1             # r* - sum of z and 4* g
  #  R[6,5] <- p[8]  
  
  #  R[13,1] <- p[4]
  #  R[15,16] <- p[7]
  #  R[13,10] <- p[3]
  #  R[10,1] <- p[3]
  
  
  W(NRDLM) <- R%*%t(R) 
  
  
  
  m0(NRDLM) <- c(0,
                 0,
                 0,
                 NRdata$log.output[1],
                 mean(diff(NRdata$log.output)),
                 NRdata$real.r[3],
                 0,
                 NRdata$real.r[2],
                 0,
                 NRdata$real.r[2],
                 NRdata$real.r[1],
                 NRdata$unr[1],
                 0,
                 0,
                 NRdata$Inflation.l1[1],
                 NRdata$Inflation.e[1])
  
  
  return(NRDLM)
  
}

theta <- c(1.5,-0.5,-0.8,-0.3, 0.3,0.7,0.3, 4,0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5) 

# Estimate model

NRDLM.est <-  dlmMLE(y = cbind(NRdata$log.output,NRdata$unr,NRdata$Inflation,NRdata$real.r,NRdata$Inflation.e), parm = theta, build = buildNRDLM, lower =c(rep(-Inf,9),rep(-8,9)), upper= c(rep(Inf,2),-0.03,rep(Inf,6),rep(12,9)),
                     control = list(trace = 6, REPORT = 5, maxit = 500), hessian = TRUE)

delta <- diag(sqrt(diag(exp(NRDLM.est$par))%*%solve(NRDLM.est$hessian)%*%diag(exp(NRDLM.est$par))))

# Create table of estimates, standard erros, and t-stats
tibble(Var =c("a1", "a2", "a3", "O.1","b1","b2","b3","d","z1","sigma ygap", "sigma y*", "sigma g","sigma r", "sigma u*","sigma ugap", "sigma pi", "sigma piexp","sigmaz"),
       Par =round(NRDLM.est$par,3),
       SE = sqrt(diag(solve(NRDLM.est$hessian))),
       SE2 = delta,
       `t-stat` = Par/SE) %>% 
  mutate(SE = ifelse(grepl("*sigma", .$Var), delta, SE),
         `z-stat`= ifelse(grepl("*sigma", .$Var), NA, Par/SE),
         Par = ifelse(grepl("*sigma", .$Var), exp(Par), Par),
  )


#--------------------------------------------------------------------------------------------------------------------------
# filtered and smoothed estimates
#--------------------------------------------------------------------------------------------------------------------------

NRDLMbuilt <- buildNRDLM(NRDLM.est$par)

filtered <- dlmFilter(y =cbind(NRdata$log.output,NRdata$unr,NRdata$Inflation,NRdata$real.r,NRdata$Inflation.e), mod = NRDLMbuilt)

smoothed <- dlmSmooth(y = cbind(NRdata$log.output,NRdata$unr,NRdata$Inflation,NRdata$real.r,NRdata$Inflation.e), mod = NRDLMbuilt)


cbind(filtered$y[-c(1:9),4],filtered$m[-c(1:10),6]) %>% matplot(type ="l")

