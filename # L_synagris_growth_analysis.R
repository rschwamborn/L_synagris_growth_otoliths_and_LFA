
### Analysis of Lutjanus synagris body growth
### using novel bootstrapped length-frequency 
### and length-at-age (otolith) analyses

# Ralf Schwamborn 
# March 20, 2022

###


# This script contains analyses conducted for the manuscript:

# "Comparing the accuracy and precision 
# of novel bootstrapped length-frequency 
# and length-at-age (otolith) analyses, 
# with a case study of lane snapper (Lutjanus synagris) 
# from the SW Atlantic"

# Authors of the manuscript:

# Ralf Schwamborn, 
# Matheus Oliveira Freitas,
# Rodrigo Leao de Moura
# and
# Alexandre Aschenbrenner




## 0. Clean memory ------------------------------------------------------------

# gc()  
# gc(reset=T)
# rm(list = ls()) 
#Ctrl+Shift+F10 clean memory

opar <- par() # save plot parameters

## 0. Load packages -----------------------------------------------------------

 library(parallel)
 library(TropFishR)
 library(ks)
 library(rfishbase)
  library(caroline)
 library(beanplot)

library(devtools)
install_github("rschwamborn/fishboot")
library(fishboot)


# 0. Function definition ---------------------------------------------------------

#### Define the function "lenage_boot" ----------------------------------
# lenage_boot
# Bootstrapped length-at-age analysis 
# input: lengths and ages (e.g., from otolith readings) 
# output: bootstrap posteriors of VGBF parameters
# part of the "fishboot" R package 

lenage_boot <- function(input.data, nboot = 200) {
  
  data.len.age <- input.data
  
  B = nboot ## number of bootstraps
  res = data.frame(Linf = numeric(B) , K = numeric(B), t0 = numeric(B)) ## vector to hold results
  n = length(data.len.age$length) # number of data pairs
  
  
  for(b in 1:B){
    
    
    tryCatch({   
      
      seed <- round(as.numeric(Sys.time())+runif(1, min = 0, max = 1e4)) # maximize stochasticity
      set.seed(seed)
      
      i = sample(x = 1:n, size = n, replace = TRUE) ## sample indices
      
      bootsam.dataxy <- data.len.age[i,] ## get data
      bootsam.dataxy <- list(age=  round(bootsam.dataxy$age,0) , length= bootsam.dataxy$length)
      
      
      output <- growth_length_age(param = bootsam.dataxy, method = "LSM",
                                  Linf_init = (max(bootsam.dataxy$length)), CI = FALSE, age_plot=NULL)
      
      output$Linf
      
      ## store results
      res$Linf[b]  <-   output$Linf 
      res$K[b]  <- output$K  
      res$t0[b]<- output$t0
      
    }, error=function(e){})
    
  }
  
  ret <- list()
  
  ret$bootRaw <- res
  
  ret$seed <- seed
  
  class(ret) <- "lfqBoot"
  
  
  
  return(ret)
  
  
}


## 1. Set working directory and load LFD data --------------------------------

 
 # setwd("D:/Documentos/Papers_Pesca_Ralph_Methods_R_Elephan")

# setwd("C:/Users/Ralf1/Desktop/Diversos_Work/Papers/00000000 - Ale - Paper Otolitos vs LFA Lutjanus")


# setwd ("C:/Users/User/Desktop/Ale_Paper_L_synagris/")
 
 
 # data.3 <- read.csv(file = "CopyPlanilha_LT_Synagris_Completa_n_2568b.csv", sep = ",")

 # data.3 <- read.csv(file = "Total_lengths_L_synagris_n_2568.csv", sep = ",")

urlfile<- 'https://raw.githubusercontent.com/rschwamborn/L_synagris_growth_otoliths_and_LFA/main/Total_lengths_L_synagris_n_2568.csv'
 
 data.3 <-read.csv(url(urlfile), sep = ",")
 
 head(data.3)
  
# data.3
 
 summary(data.3)
  
 
 
 
 ## 2. Traditional  LFD Analysis (non-bootstrapped) ####
 
 
 data.3$Data <- as.Date(data.3$Date, tryFormats = c("%d-%m-%Y"))
 
 data.3$Comp <- as.numeric(data.3$Total_Length_cm) 
 
 summary(data.3$Data )
 

 length(data.3$Comp ) # n = 2564 
 
 summary(data.3$Comp) # Total Lengths
 
 library(TropFishR)
 
 
 # ?lfqCreate
 
 # create LFD data
 
 lfq_dat.syn  <- lfqCreate(data.3,Lname = "Comp", Dname = "Data", aggregate_dates = TRUE,
                           length_unit = "cm", bin_size = 2, plot=TRUE)
 
 # test different Bin sizes
 
 lfq_dat.syn1  <- lfqCreate(data.3,Lname = "Comp", Dname = "Data", aggregate_dates = TRUE,
                            length_unit = "cm", bin_size = 1, plot=TRUE)#, main = "Bin size = 1 cm")
 
 lfq_dat.syn2  <- lfqCreate(data.3,Lname = "Comp", Dname = "Data", aggregate_dates = TRUE,
                           length_unit = "cm", bin_size = 2, plot=TRUE)#,, main = "Bin size = 2 cm")
 
 lfq_dat.syn3  <- lfqCreate(data.3,Lname = "Comp", Dname = "Data", aggregate_dates = TRUE,
                           length_unit = "cm", bin_size = 3, plot=TRUE)#,, main = "Bin size = 3 cm")
 
 lfq_dat.syn4  <- lfqCreate(data.3,Lname = "Comp", Dname = "Data", aggregate_dates = TRUE,
                           length_unit = "cm", bin_size = 4, plot=TRUE)#,, main = "Bin size = 4 cm")
 lfq_dat.syn5  <- lfqCreate(data.3,Lname = "Comp", Dname = "Data", aggregate_dates = TRUE,
                            length_unit = "cm", bin_size = 5, plot=TRUE)#,, main = "Bin size = 5 cm")
 lfq_dat.syn6  <- lfqCreate(data.3,Lname = "Comp", Dname = "Data", aggregate_dates = TRUE,
                            length_unit = "cm", bin_size = 6, plot=TRUE)#,, main = "Bin size = 6 cm")
 lfq_dat.syn8  <- lfqCreate(data.3,Lname = "Comp", Dname = "Data", aggregate_dates = TRUE,
                            length_unit = "cm", bin_size = 8, plot=TRUE)#,, main = "Bin size = 8 cm")
 
 
 
 sum(lfq_dat.syn$catch) # OK
 
 lfqbin <- lfqRestructure(lfq_dat.syn, MA = 5,addl.sqrt = FALSE)
 opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
 plot(lfqbin, Fname = "catch", date.axis = "modern")
 plot(lfqbin, Fname = "rcounts", date.axis = "modern")
 par(opar)
 
 
 #View(lfq_dat.syn$midLengths)
 sum(lfq_dat.syn$catch)
 
 
 # 2. Traditional LFA (PW-plot and K-scan) ------------------------------------
 
 
 #Powell Wetherall plot
 res_PW <- powell_wetherall(param = lfq_dat.syn,
                            catch_columns = 1:ncol(lfq_dat.syn$catch),
                            reg_int = c(10,30))
 
 # show results
 paste("Linf =",round(res_PW$Linf_est), "?", round(res_PW$se_Linf))
 #> [1] "Linf = 56 +- 2"
 
 # ELEFAN with K-Scan -> Shows the K up and down with fixed Linf -> No seasonal growth
 res_KScan <- ELEFAN(lfq_dat.syn, Linf_fix = res_PW$Linf_est,
                     MA=7, addl.sqrt = TRUE, hide.progressbar = TRUE, plot = FALSE)
 
 
 # show results
 res_KScan$par; res_KScan$Rn_max
 
 # $Linf
 # [1] 56.24408
 # 
 # $K
 # [1] 3.13
 # 
 # $t_anchor
 # [1] 0.6054715
 # 
 # $C
 # [1] 0
 # 
 # $ts
 # [1] 0
 # 
 # $phiL
 # [1] 3.995698
 # 
 
 
 # 3. TropFishR, modern search algorithm ELEFAN_SA  --------------------------
 
 # ELEFAN_SA ------------------------
 
# iterations using agemax = 18 years (maximum age from otoliths)

 set.seed(3)

 res_SA <- ELEFAN_SA(lfqbin,   SA_time = 60*2, SA_temp = 6e5,
                     agemax = 18,
                     MA = 5, seasonalised = TRUE, addl.sqrt = FALSE,
                     init_par = list(Linf = 50, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                     low_par = list(Linf = 30, K = 0.01, t_anchor = 0, C = 0, ts = 0),
                     up_par = list(Linf = 100, K = 1, t_anchor = 1, C = 1, ts = 1))
 # show results
 res_SA$par; res_SA$Rn_max


 ## resuts of ELEFAN_SA
 # $Linf
 # [1] 72.10377
 # 
 # $K
 # [1] 0.1472316
 # 
 # $t_anchor
 # [1] 0.7950251
 # 
 # $C
 # [1] 0.8873571
 # 
 # $ts
 # [1] 0.4050936
 # 
 # $phiL
 # [1] 2.883917
 # 
 # [1] 0.2242769
 # 
 #plot curve
 
 # plot(res_SA, draw = FALSE)
 # lfqFitCurves(res_SA, col="blue", par=res_SA$par, draw=TRUE)$ESP
 
 # Best fit (Rn= 0.567)
 # 
 # $Linf
 # [1] 410.59 mm
 # 
 # $K
 # [1] 0.01350 
 # 
 # $t_anchor
 # [1] 0.7921
 # 
 # $C
 # [1] 0.850
 # 
 # $ts
 # [1] 0.277
 # 
 # $phiL
 # [1] 3.36
 # 
 # Rn
 # [1] 0.567
 
 
 
 ###
 ## now restrain with otlith results, Linf = 60, K = 0.2 ####
 # 
 # res_SA_restr <- ELEFAN_SA(lfqbin,   SA_time = 60*2, SA_temp = 6e5,
 #                           agemax = 18,
 #                           MA = 5, seasonalised = TRUE, addl.sqrt = FALSE,
 #                           init_par = list(Linf = 60, K = 0.2, t_anchor = 0.5, C=0.5, ts = 0.5),
 #                           low_par = list(Linf = 55, K = 0.15, t_anchor = 0, C = 0, ts = 0),
 #                           up_par = list(Linf = 65, K = 0.25, t_anchor = 1, C = 1, ts = 1))
 # # show results
 # res_SA_restr$par; res_SA_restr$Rn_max
 # 
 #  Show Results ###
 # $Linf 
 # [1] 59.78423
 # 
 # $K
 # [1] 0.1907338
 # 
 # $t_anchor
 # [1] 0.9023946
 # 
 # $C
 # [1] 0.7749233
 # 
 # $ts
 # [1] 0.5693673
 # 
 # $phiL
 # [1] 2.833601
 # 
 # [1] 0.2210604
 
 #plot curve
 # 
 # plot(res_SA_restr, draw = FALSE)
 # lfqFitCurves(res_SA_restr, col="blue", par=res_SA$par, draw=TRUE)$ESP
 # 
 
 #####
 # 
 # 
 # Lmax.3 <- max(data.3$Comp.) # Lmax = 55.8
 # 
 # min_Lmaxfor_search <-  Lmax.3 * 0.75 # 41.9 cm
 # max_Lmaxfor_search <-  Lmax.3 * 2  # 111.6 cm
 # 
 # library(TropFishR)
 # 
 # # set.seed(6)
 # set.seed(as.numeric(Sys.time()))
 # 
 # seed.num.dummy <- (as.numeric(Sys.time()) + runif(1,10,1000))
 # 
 # 
 # 
 
 # # results   of ELEFAN_GA (MA = 7)
 
 ### run 1 (MA= 7)
 # $Linf
 # [1] 65.16758
 # 
 # $K
 # [1] 0.187528
 # 
 # $t_anchor
 # [1] 0.2814865
 # 
 # $C
 # [1] 0.5534256
 # 
 # $ts
 # [1] 0.336149
 # 
 # $phiL
 # [1] 2.901129
 # 
 # >  output5$ASP
 # [1] 32.76644
 # >  output5$Rn_max
 # [1] 0.201711
 # 
 
 #### run 2 (MA = 7)
 
 # $Linf
 # [1] 62.42689
 # 
 # $K
 # [1] 0.240235
 # 
 # $t_anchor
 # [1] 0.6398544
 # 
 # $C
 # [1] 0.4458948
 # 
 # $ts
 # [1] 0.3472186
 # 
 # $phiL
 # [1] 2.97138
 # 
 # >  output5$ASP
 # [1] 32.76644
 # >  output5$Rn_max
 # [1] 0.1737833
 # > 
 # 
 
### run 3 (MA = 7)
#  $Linf
#  [1] 73.51738
#  
#  $K
#  [1] 0.1663429
#  
#  $t_anchor
#  [1] 0.3699263
#  
#  $C
#  [1] 0.6721292
#  
#  $ts
#  [1] 0.3080939
#  
#  $phiL
#  [1] 2.953784
#  
#  >  output5$ASP
#  [1] 32.76644
#  >  output5$Rn_max
#  [1] 0.2326851
#  > 
 
 
 
 
 
 
 ## 4. Bootstrapped LFD Analysis (fishboot) ----------------------------------
 
 
 ##  bootstrapped LFD Analysis (fishboot)  with MA = 5,7,9,11,13  
 
 
 library(fishboot)
 #?ELEFAN_GA_boot()
 
 # 
 # set.seed(as.numeric(Sys.time()))
 # 
 # seed.num.dummy <- (as.numeric(Sys.time()) + runif(1,10,1000))
 # 
 # 
 # t1 <- Sys.time()
 # 
 # output_E_MA5_200l <- ELEFAN_GA_boot(lfqbin, nresamp = 200, seasonalised = TRUE,
 #                      low_par = list(Linf = 41.9, K = 0.01, t_anchor = 0, C = 0, ts= 0),
 #                      up_par = list(Linf = 111.6, K = 0.6, t_anchor = 1, C = 1, ts = 1),
 #                      popSize = 40, maxiter = 30, run = 10,
 #                      MA = 5,  seed = seed.num.dummy)
 # 
 # t2 <- Sys.time()
 # 
 # t2-t1 
 # 
 
 ##### Analyse bootstrap results #####
 
 # Experiment A
 # # Experim A = MA7, Binsize 2 cm, seasonal (C not fixed)
 # duration non-bootstrap: 56 secs, 38 secs, 60 secs, per run, with PC_six_core_Ralf_sala_Museu
 # duration GA_boot: 27 min for  100 runs, with PC_six_core_Ralf_sala_Museu, run 7a, file: output_MA7_100a.csv
 # duration GA_boot: 53 min for  200 runs, with PC_six_core_Ralf_sala_Museu, run 7b, file: output_MA7_200b.csv
 # duration GA_boot: 1.3 hs for  300 runs, with PC_six_core_Ralf_sala_Museu, run 7c, file: output_MA7_300c.csv
 # average duration: aprox. 25-30 mins for 100 runs 

 
 #### Combine experiment results into one big posterior (nboot = 600) ####
 
 # # MA 7
 # output_A_MA7_100a.csv <-read.csv("output_A_MA7_100a.csv")
 # output_A_MA7_200b.csv <-read.csv("output_A_MA7_200b.csv")
 # output_A_MA7_300c.csv <-read.csv("output_A_MA7_300c.csv")
 # 
 # EXP_A_MA7 <- rbind(output_A_MA7_200b.csv, output_A_MA7_300c.csv)
 # EXP_A_MA7_500 <- rbind(output_A_MA7_200b.csv, output_A_MA7_300c.csv)
 # EXP_A_MA7_600 <- rbind(output_A_MA7_100a.csv, output_A_MA7_200b.csv, output_A_MA7_300c.csv)
 # EXP_A_MA7_400 <- EXP_A_MA7_600[1:400,]
 # 
 # head(EXP_A_MA7) 
 # 
 # 
 # # View(EXP_A_MA7_400)
 # 
 # # MA 9
 # 
 #  output_B_MA9_200f <-read.csv("output_B_MA9_200f.csv")
 # output_B_MA9_400e <-read.csv("output_B_MA9_400e.csv")
 # 
 # EXP_B_MA9 <- rbind(output_B_MA9_200f, output_B_MA9_400e)
 # EXP_B_MA9_600 <- EXP_B_MA9
 # EXP_B_MA9_400 <- EXP_B_MA9_600[1:400,]
 # 
 #  
 # head(EXP_B_MA9) 
 # 
 # View(EXP_B_MA9_600)
 
 ## MA 11
 
 # output_C_MA11_200h <-read.csv("output_C_MA11_200h.csv")
 # output_C_MA11_400g <-read.csv("output_C_MA11_400g.csv")
 # 
 # EXP_C_MA11 <- rbind(output_C_MA11_200h, output_C_MA11_400g)
 # EXP_C_MA11_600 <- EXP_C_MA11
 # EXP_C_MA11_400 <- EXP_C_MA11_600[1:400,]
 # 
 # 
 # head(EXP_C_MA11_400) 
 # head(EXP_C_MA11_600) 
 # 
 # View(EXP_C_MA11_400)
 
 
 ## MA 13
 # 
 # output_D_MA13_200i <-read.csv("output_D_MA13_200i.csv")
 # output_D_MA13_400j <-read.csv("output_D_MA13_400j.csv")
 # 
 # EXP_D_MA13 <- rbind(output_D_MA13_200i, output_D_MA13_400j)
 # EXP_D_MA13_600 <- EXP_D_MA13
 # EXP_D_MA13_400 <- EXP_D_MA13_600[1:400,]
 # 
 # 
 # head(EXP_D_MA13_400) 
 # head(EXP_D_MA13_600) 
 # 
 # View(EXP_D_MA13_600)
 
 
 # # MA 5
 
 # setwd("C:/Users/Ralf1/Desktop/Diversos_Work/Papers/00000000 - Ale - Paper Otolitos vs LFA Lutjanus/Data_inputs_outputs")
 
 
 # output_E_MA5_200k <-read.csv("output_E_MA5_200k.csv")
 # output_E_MA5_200l <-read.csv("output_E_MA5_200l.csv")
 # 
 # EXP_E_MA5_400 <- rbind(output_E_MA5_200k, output_E_MA5_200l)
 # 
 # write.csv (EXP_E_MA5_400, file = "output_E_MA5_400.csv") 
 
 
 # EXP_E_MA5_400 <- read.csv("output_E_MA5_400.csv")
 # EXP_E_MA5_400 <- <- read.csv(url("http://some.where.net/data/foo.csv"))
 
 urlfile2 <- 'https://raw.githubusercontent.com/rschwamborn/L_synagris_growth_otoliths_and_LFA/main/output_E_MA5_400.csv'
 
 EXP_E_MA5_400 <- read.csv(url(urlfile2))
 
 

 head(EXP_E_MA5_400)

 # View(EXP_E_MA5_400)

 summary(EXP_E_MA5_400)
 
 # MEDIAN values (MA = 5)
 # 
 # Linf 70.46
 # 
 # K  0.1531
 # 
 # t_anchor  0.47
 # 
 # C 0.5270
 # 
 # ts 0.50096
 # 
 # phiL 2.902 
 
 
 #plot curve
 
 
 #  bootstrapped LFA - Medians and 95% Confidence Intervals  (MA = 5)  ---------------
 
 quantile(EXP_E_MA5_400$Linf, c(0.5, 0.025, 0.975))
 quantile(EXP_E_MA5_400$K, c(0.5, 0.025, 0.975))
 quantile(EXP_E_MA5_400$phiL, c(0.5, 0.025, 0.975))
 quantile(EXP_E_MA5_400$t_anchor, c(0.5, 0.025, 0.975))
 quantile(EXP_E_MA5_400$C, c(0.5, 0.025, 0.975))
 quantile(EXP_E_MA5_400$ts, c(0.5, 0.025, 0.975))
 
 #  LAA - Medians and 95% Confidence Intervals  (otoliths)  ---------------
 
# setwd("C:/Users/RALF/Desktop/Ralf_diversos 2018_19_iii/Papers/00 - Ale - Paper Otolitos vs LFA Lutjanus/Data_inputs_outputs")
 
 # EXP_Otol_len_age_400 <-read.csv("EXP_Otol_len_age_400.csv")
 
 urlfile3 <- 'https://raw.githubusercontent.com/rschwamborn/L_synagris_growth_otoliths_and_LFA/main/EXP_Otol_len_age_400.csv'
 
 EXP_Otol_len_age_400 <- read.csv(url(urlfile3))
 
 head(EXP_Otol_len_age_400)
 summary(EXP_Otol_len_age_400)
 
 
 # Phi_p = log10(K) + 2 log10(Linf)
 EXP_Otol_len_age_400$Phi_vec <- log10(EXP_Otol_len_age_400$K) + 2 * log10(EXP_Otol_len_age_400$Linf)
 quantile(EXP_Otol_len_age_400$Phi_vec, c(0.5, 0.025, 0.975))
 
 
 quantile(EXP_Otol_len_age_400$Linf, c(0.5, 0.025, 0.975))
 quantile(EXP_Otol_len_age_400$K, c(0.5, 0.025, 0.975))
 quantile(EXP_Otol_len_age_400$Phi_vec, c(0.5, 0.025, 0.975))
 
 
 
 
 
 
 
 res_MA5_OK <- res_SA

 res_MA5_OK$par$Linf <- 70.46
 
 res_MA5_OK$par$K <- 0.1531
 
 res_MA5_OK$par$t_anchor <- 0.47
 
 res_MA5_OK$par$C <- 0.5270
 
 res_MA5_OK$par$ts <- 0.50096
  
 res_MA5_OK$par$phiL <- 2.902 
  
  
 plot(res_MA5_OK, draw = FALSE)
 lfqFitCurves(res_MA5_OK, col="blue", par=res_MA5_OK$par, draw=TRUE)$ESP
 
 res_MA5_OK$par
 
 
 ##### Add columns with MA and BS (bin size)  and make one huge file #### 
 
 # EXP_A_MA7_400$MA <- rep(7, 400);   EXP_A_MA7_400$BS <- rep(2, 400) 
 # EXP_B_MA9_400$MA <- rep(9, 400);   EXP_B_MA9_400$BS <- rep(2, 400) 
 # EXP_C_MA11_400$MA <- rep(11, 400); EXP_C_MA11_400$BS <- rep(2, 400) 
 # EXP_D_MA13_400$MA <- rep(13, 400); EXP_D_MA13_400$BS <- rep(2, 400) 
 # EXP_E_MA5_400$MA <- rep(5, 400);   EXP_E_MA5_400$BS <- rep(2, 400) 
 # 
 # EXP_ABCDE_MA5791113_400 <- rbind (EXP_E_MA5_400, EXP_A_MA7_400,EXP_B_MA9_400,EXP_C_MA11_400,
 #                                   EXP_D_MA13_400)
 
 # View(EXP_ABCDE_MA5791113_400)
 
 
 # EXP_A_MA7_600$MA <- rep(7, 600);   EXP_A_MA7_600$BS <- rep(2, 600) 
 # EXP_B_MA9_600$MA <- rep(9, 600);   EXP_B_MA9_600$BS <- rep(2, 600) 
 # EXP_C_MA11_600$MA <- rep(11, 600); EXP_C_MA11_600$BS <- rep(2, 600) 
 # EXP_D_MA13_600$MA <- rep(13, 600); EXP_D_MA13_600$BS <- rep(2, 600) 
 # EXP_E_MA5_600$MA <- rep(5, 600);   EXP_E_MA5_600$BS <- rep(2, 600) 
 # 
 # EXP_ABCDE_MA5791113_600 <- rbind (EXP_E_MA5_600, EXP_A_MA7_600,EXP_B_MA9_600,EXP_D_MA11_600,
 #                                   EXP_D_MA13_600)
 
 
 #### write and read summary files#### 
 
 # getwd()
 
 # write.csv(EXP_E_MA5_400, file = "EXP_E_MA5_400.csv")
 # write.csv(EXP_A_MA7_400, file = "EXP_A_MA7_400.csv")
 # write.csv(EXP_B_MA9_400, file = "EXP_B_MA9_400.csv")
 # write.csv(EXP_C_MA11_400, file = "EXP_C_MA11_400.csv")
 # write.csv(EXP_D_MA13_400, file = "EXP_D_MA13_400.csv")
 # 
 # write.csv(EXP_ABCDE_MA5791113_400, file = "EXP_ABCDE_MA5791113_400.csv")
 # 
 
 # EXP_E_MA5_400 <-read.csv("EXP_E_MA5_400.csv")
 # EXP_A_MA7_400 <-read.csv("EXP_A_MA7_400.csv")
 # EXP_B_MA9_400 <-read.csv("EXP_B_MA9_400.csv")
 # EXP_C_MA11_400 <-read.csv("EXP_C_MA11_400.csv")
 # EXP_D_MA13_400 <-read.csv("EXP_D_MA13_400.csv")
 # 
 #EXP_ABCDE_MA5791113_400 <-read.csv("EXP_ABCDE_MA5791113_400.csv")
 
 urlfile4 <- 'https://raw.githubusercontent.com/rschwamborn/L_synagris_growth_otoliths_and_LFA/main/EXP_ABCDE_MA5791113_400.csv'
 
 EXP_ABCDE_MA5791113_400 <- read.csv(url(urlfile4))
 
 head(EXP_ABCDE_MA5791113_400)
 summary(EXP_ABCDE_MA5791113_400)
 
#####  Add otolith data (posteriors) #####
  
#  compare Conf interval widths -------------

 head(EXP_ABCDE_MA5791113_400)
 
 attach(EXP_ABCDE_MA5791113_400)
library(beanplot)  

 
 # Linf_Beanplot_MA5to13_400
 
 beanplot(Linf ~ as.factor(EXP_ABCDE_MA5791113_400$MA),
          what=c(0,1,1,0),
          boxwex = 0.9,
 col = c("grey", "lightgrey") , cut = 0.001, beanlines = "median", 
          border = c("darkgrey"), xlab = "Moving average, MA (BS = 2), lines: medians", 
 ylab = "Linf", ylim = c(40, 120))

 ## check 95% conf ints
 
 # 
 # quantile(EXP_B_MA9_400$Linf, probs = c(0.025,0.5,0.975)) # 95% CIs
 # abline(h = c(44.6, 105.51))
 # 
 
 # K_Beanplot_MA5to13_400
 
 beanplot(EXP_ABCDE_MA5791113_400$K ~ as.factor(EXP_ABCDE_MA5791113_400$MA),
          what=c(0,1,1,0),
          col = c("grey", "lightgrey") , cut = 0.95, beanlines = "median", 
          border = c("darkgrey"), xlab = "Moving average, MA (BS = 2)", 
          ylab = "K", ylim = c(0.02, 1))
 
 # C_Beanplot_MA5to13_400
 
 beanplot(EXP_ABCDE_MA5791113_400$C ~ as.factor(EXP_ABCDE_MA5791113_400$MA),
          what=c(0,1,1,0),
          col = c("grey", "lightgrey") , cut = 0.95, beanlines = "median", 
          border = c("darkgrey"), xlab = "Moving average, MA (BS = 2)", 
          ylab = "C", ylim = c(0.02, 1))
 
 
 ## now a new vioplot with 95 % quantiles as dotted box ("violins" in "caroline" package ######
 # 
 # # install.packages("caroline")
 # library(caroline)
 # #? violins
 # 
 # opar <- par()
 # 
 # # Linf_violinsplot_MA5to13_400
 # 
 # violins(list(MA5=EXP_E_MA5_400$Linf, MA7=EXP_A_MA7_400$Linf,MA9=EXP_B_MA9_400$Linf,
 #              MA11=EXP_C_MA11_400$Linf,MA13=EXP_D_MA13_400$Linf), 
 #            ylim=c(40,120),
 #         col=c('purple','lightblue','lightgreen','red','orange'),
 #         stats=TRUE, 
 #         quantiles = c(0.025, 0.975),
 #         connect = FALSE,
 #         drawRect = TRUE,
 #         deciles = F,
 #         CImed = F,
 #         ylab = "Linf (cm)")
 # 
 # quantile(EXP_B_MA9_400$Linf, probs = c(0.025,0.5,0.975)) # 95% CIs
 # abline(h = c(44.6, 105.51))
 # 
 # # K_violinsplot_MA5to13_400
 # 
 # 
 # violins(list(MA5=EXP_E_MA5_400$K, MA7=EXP_A_MA7_400$K,MA9=EXP_B_MA9_400$K,
 #              MA11=EXP_C_MA11_400$K,MA13=EXP_D_MA13_400$K), 
 #         ylim=c(0.04,0.65),
 #         col=c('purple','lightblue','lightgreen','red','orange'),
 #         stats=TRUE, 
 #         quantiles = c(0.025, 0.975),
 #         connect = FALSE,
 #         drawRect = TRUE,
 #         deciles = F,
 #         CImed = F,
 # #         ylab = "K (y-1)")
 # # 
 # # 
 # # # Now with otolith posteriors in one graph
 # # 
 # #  
 # setwd("C:/Users/RALF/Desktop/Ralf_diversos 2018_19_iii/Papers/00 - Ale - Paper Otolitos vs LFA Lutjanus/Data_inputs_outputs")
 # 
 # EXP_Otol_len_age_400 <-read.csv("EXP_Otol_len_age_400.csv")
 # 
 # #View(EXP_Otol_len_age_400) 
 # 
 # # clean (remove Linf = zero)
 # 
 # EXP_Otol_len_age_400<-EXP_Otol_len_age_400[!(EXP_Otol_len_age_400$Linf==0 ),]
 # 
 # #View(EXP_Otol_len_age_400)
 # 
 # length(EXP_Otol_len_age_400$Linf) # 399 data
 # 
 # # Linf_violinsplot_MA5to13_400_Otoliths
 # 
 # violins(list(Otoliths =EXP_Otol_len_age_400$Linf, MA5=EXP_E_MA5_400$Linf, MA7=EXP_A_MA7_400$Linf,MA9=EXP_B_MA9_400$Linf,
 #              MA11=EXP_C_MA11_400$Linf,MA13=EXP_D_MA13_400$Linf), 
 #         ylim=c(40,120),
 #         col=c('green','purple','lightblue','lightgreen','red','orange'),
 #         stats=TRUE, 
 #         quantiles = c(0.025, 0.975),
 #         connect = FALSE,
 #         drawRect = TRUE,
 #         deciles = F,
 #         CImed = F,
 #         ylab = "Linf (cm)")
 # 
 # # quantile(EXP_B_MA9_400$Linf, probs = c(0.025,0.5,0.975)) # 95% CIs
 # # abline(h = c(44.6, 105.51))
 # 
 # # K_violinsplot_MA5to13_400_Otoliths
 # 
 # 
 # violins(list(Otoliths =EXP_Otol_len_age_400$K, MA5=EXP_E_MA5_400$K, MA7=EXP_A_MA7_400$K,MA9=EXP_B_MA9_400$K,
 #              MA11=EXP_C_MA11_400$K,MA13=EXP_D_MA13_400$K), 
 #         ylim=c(0.04,0.65),
 #         col=c('green','purple','lightblue','lightgreen','red','orange'),
 #         stats=TRUE, 
 #         quantiles = c(0.025, 0.975),
 #         connect = FALSE,
 #         drawRect = TRUE,
 #         deciles = F,
 #         CImed = F,
 #         ylab = "K (y-1)")
 # 
 # 
 # 
 

 # Experiment B
 # # Experim B = MA9, Binsize 2 cm, seasonal (C not fixed)
 # duration GA_boot: 1.5 hs for  400 runs, with PC_six_core_Ralf_sala_Museu, file: output_B_MA9_400e.csv
 # duration GA_boot:  47 mins  200 runs, with PC_six_core_Ralf_sala_Museu, file:output_B_MA9_200f.csv
 # average duration: aprox. 20-25 mins for 100 runs 
 

 
 # Experiment C
 # # Experim C = MA11, Binsize 2 cm, seasonal (C not fixed)
 # duration GA_boot: 1.5 hs  for  400 runs, with PC_six_core_Ralf_sala_Museu, run 7g, file: output_C_MA11_400g.csv
 # duration GA_boot: 48 min.  for  200 runs, with PC_six_core_Ralf_sala_Museu, run 7h, file: output_C_MA11_200h.csv
 # average duration: aprox. 20-25 mins for 100 runs 
   
 
 # Experiment D
 # # Experim D = MA13, Binsize 2 cm, seasonal (C not fixed)
 # duration GA_boot: 41 Min.  for  200 runs, with PC_six_core_Ralf_sala_Museu, run 7i, file: output_D_MA13_200i.csv
 # duration GA_boot: 1.28 h  for  400 runs, with PC_six_core_Ralf_sala_Museu, run 7j, file: output_D_MA13_400j.csv
 # average duration: aprox. 20-23 mins for 100 runs 
 
 
 # Experiment E
 # # Experim D = MA5, Binsize 2 cm, seasonal (C not fixed)
 # duration GA_boot: 1 h  for  200 runs, with PC_six_core_Ralf_sala_Museu, run 7i, file: output_E_MA5_200k.csv
 # duration GA_boot: xx min.  for  400 runs, with PC_six_core_Ralf_sala_Museu, run 7j, file: output_E_MA5_200l.csv
 # average duration: aprox. 30 mins for 100 runs 
 
 
 
 
 
 
 
  
 #outputs
 
 # res7a <- output_E_MA5_200l$bootRaw
 
 #res7a <- EXP_A_MA7_400
 # 
 # hist(res7a$Linf)
 # summary(res7a$Linf)
 # quantile(res7a$Linf, probs = c(0.025,0.5,0.975)) # 95% CIs
 # 
 # hist(res7a$K)
 # summary(res7a$K)
 # quantile(res7a$K, probs = c(0.025,0.5,0.975)) # 95% CIs
 # 
 # hist(res7a$t_anchor)
 # summary(res7a$t_anchor)
 # quantile(res7a$t_anchor, probs = c(0.025,0.5,0.975)) # 95% CIs
 # 
 # hist(res7a$C)
 # summary(res7a$C)
 # quantile(res7a$C, probs = c(0.025,0.5,0.975)) # 95% CIs
 # 
 # hist(res7a$ts)
 # summary(res7a$ts)
 # quantile(res7a$ts, probs = c(0.025,0.5,0.975)) # 95% CIs
 # 
 # hist(res7a$phiL)
 # summary(res7a$phiL)
 # quantile(res7a$phiL, probs = c(0.025,0.5,0.975)) # 95% CIs
 # 
 # 
 ####  write.csv() -> Save results do Hard Disk ####
 
 # 
 # # write.csv(res7a, file = "output_E_MA5_200l.csv")
 # 
 #  write.csv(res7a, file = "EXP_A_MA7.csv")
 # 
 # 
 ########
 ########
 ########
 ########
 
 
 
 
 
 # 5. Bootstrapped length_at_age (from otoliths)   ----------------------------
 
### Function Len.age_boot  V.06
## otolith length_at_age bootstrap V.05
# Copyrigth: Ralf Schwamborn, 2020
# includes Len.age_boot function

# Analysis of L. synagris otolith data

#####  Bootstrap the VBGF fit to length-at-age data  #####

#### load  data ####

 # setwd("C:/Users/Ralf1/Desktop/Diversos_Work/Papers/00000000 - Ale - Paper Otolitos vs LFA Lutjanus")
 
#Lsynagris_lenage <- read.csv("Idade_Otolitos_Lsynagris.csv")

urlfile5 <- 'https://raw.githubusercontent.com/rschwamborn/L_synagris_growth_otoliths_and_LFA/main/Idade_Otolitos_Lsynagris.csv'

Lsynagris_lenage <- read.csv(url(urlfile5))

head(Lsynagris_lenage)
summary(Lsynagris_lenage)

 
#   View(Lsynagris_lenage)

names(Lsynagris_lenage) <- c("age", "length"  )

length(Lsynagris_lenage$length)

summary(Lsynagris_lenage)

attach(Lsynagris_lenage)
plot(length ~ age)

## adjust a VBGF function

# non linear least squares method , TropFishR::growth_length_age,method = "LSM"

library(TropFishR)

# output2 <- growth_length_age(param = Lsynagris_lenage, method = "LSM",Linf_init = 30, CI = TRUE, age_plot=NULL)
# 
# summary(output2$mod)

# ?growth_length_age

output3 <- growth_length_age(param = Lsynagris_lenage, method = "LSM",Linf_init = (max(length)), CI = TRUE, age_plot=NULL)

summary(output3$mod)

output3$Linf
output3$K
output3$t0

# ?growth_length_age


# growth_length_age(Lsynagris_lenage, method = "GullandHolt")

# Bertalaffy plot
#growth_length_age(Lsynagris_lenage, method = "BertalanffyPlot", Linf_est = (max(length)))

# non linear least squares method

output <- growth_length_age(Lsynagris_lenage , method = "LSM",
                              CI = TRUE, age_plot=NULL,
                            Linf_init = (max(length)))
summary(output$mod)



## now bootsrap the VBGF with length-at-age  ###
#### DEFINE Function lenage_boot ####

lenage_boot <- function(input.data, nboot = 200) {

data.len.age <- input.data

B = nboot ## number of bootstraps
res = data.frame(Linf = numeric(B) , K = numeric(B), t0 = numeric(B)) ## vector to hold results
n = length(data.len.age$length) # number of data pairs


for(b in 1:B){
  
  
  tryCatch({   
  
    seed <- round(as.numeric(Sys.time())+runif(1, min = 0, max = 1e4)) # maximize stochasticity
    set.seed(seed)
    
  i = sample(x = 1:n, size = n, replace = TRUE) ## sample indices
  
  bootsam.dataxy <- data.len.age[i,] ## get data
  bootsam.dataxy <- list(age=  round(bootsam.dataxy$age,0) , length= bootsam.dataxy$length)
  
  
  output <- growth_length_age(param = bootsam.dataxy, method = "LSM",
                              Linf_init = (max(bootsam.dataxy$length)), CI = FALSE, age_plot=NULL)
  
  output$Linf
  
  ## store results
  res$Linf[b]  <-   output$Linf 
  res$K[b]  <- output$K  
  res$t0[b]<- output$t0
  
  }, error=function(e){})
  
}

ret <- list()

ret$bootRaw <- res

ret$seed <- seed

class(ret) <- "lfqBoot"


return(ret)


}



resok <- lenage_boot(Lsynagris_lenage, nboot = 2000)

# View(resok$bootRaw)

res <- resok$bootRaw


hist(res$Linf)
summary(res$Linf)
quantile(res$Linf, probs = c(0.025,0.5,0.975)) # 95% CIs

hist(res$K)
summary(res$K)
quantile(res$K, probs = c(0.025,0.5,0.975)) # 95% CIs

hist(res$t0)
summary(res$t0)
quantile(res$t0, probs = c(0.025,0.5,0.975)) # 95% CIs


####################################################
# Bootstrap results
# > res <- resok$bootRaw
# > hist(res$Linf)
# > summary(res$Linf)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   58.33   59.99   59.90   61.87   76.95 
# > quantile(res$Linf, probs = c(0.025,0.5,0.975)) # 95% CIs
# 2.5%      50%    97.5% 
# 55.53008 59.99258 68.04329 
# > 
#   > hist(res$K)
# > summary(res$K)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.1821  0.1962  0.1955  0.2113  0.3227 
# > quantile(res$K, probs = c(0.025,0.5,0.975)) # 95% CIs
# 2.5%       50%     97.5% 
# 0.1453072 0.1961775 0.2453651 
# > 
#   > hist(res$t0)
# > summary(res$t0)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.663164 -0.141050 -0.019470 -0.007274  0.120815  1.000757 
# > quantile(res$t0, probs = c(0.025,0.5,0.975)) # 95% CIs
# 2.5%         50%       97.5% 
# -0.36941694 -0.01946974  0.40829093 
########################################################################


output2 <- growth_length_age(param = Lsynagris_lenage, method = "LSM",Linf_init = 30, CI = TRUE, age_plot=NULL)

summary(output2$mod)


# Works OK!
# Conclusion: Bootstrap 95% CI are much larger (CI width more than double) than nls (TropFishR) estimates

######
# plot LFA (otolith) results, nboot = 400 ###########

res <- EXP_Otol_len_age_400 

data.len.age <- Lsynagris_lenage

attach(Lsynagris_lenage)




output3$Linf
output3$K
output3$t0


t.seq2 <- seq(0,max(age) , by = 0.1)
Linf.2 <-  round((output3$Linf),2)
K.2 <- round((output3$K),3)
t0.2 <- output3$t0
L.seq.2 <- VBGF(t = t.seq2, list(Linf=Linf.2, K=K.2, C= 0, t0 = t0.2))  # t.zero = 0 ???


plot(t.seq2, L.seq.2, t="l", col = "blue", xlim = c(0,max(age)), ylim = c(0,max(length)*1.2),
     xlab = "age (y)", ylab = "L (cm)")
points(data.len.age$length ~ data.len.age$age, col = "navyblue")


###
# plot with grey curve swarms
###


plot(t.seq2, L.seq.2, t="l", col = "green", xlim = c(0,max(age)), ylim = c(0,max(length)*1.2),
     xlab = "age (y)", ylab = "L (cm)")
points(data.len.age$length ~ data.len.age$age, col = "navyblue")


# library(TropFishR)
# calculate Lt from t, with t0
# t.seq.p = seq(0,max(agerec.2),0.1)
# Lt.p <- VBGF(list(Linf=res$bootRaw$Linf[1], 
#                   K= res$bootRaw$K[1], 
#                   t0= res$bootRaw$t_zero_i[1]),
#              t=t.seq.p)
# lines(t.seq.p, Lt.p, col = rgb(0, 0, 0, 0.25), lwd = 0.5 ) # transparent lines
# 


res <- EXP_Otol_len_age_400


library(scales)

for (w in  1:length(res$K)) {

  t.seq.p = seq(0,max(data.len.age$age),0.1)
  Lt.p <- VBGF(list(Linf=res$Linf[w],
                    K= res$K[w],
                    t0= res$t0[w]),
               t=t.seq.p)
  # lines(t.seq.p, Lt.p, col = rgb(0, 0, 0, 0.2), lwd = 0.3 ) # transparent lines

  lines(t.seq.p, Lt.p, col = alpha("grey", 0.1), lwd = 0.3 ) # transparent lines

}

points(data.len.age$length ~ data.len.age$age, col = "navyblue")
#points(data.len.age$length ~ data.len.age$age, col = "green")

### nice plot

plot(t.seq2, L.seq.2, t="l", col = "green", xlim = c(0,max(age)), ylim = c(0,max(length)*1.2),
     xlab = "age (y)", ylab = "L (cm)")
lines(t.seq2, L.seq.2, t="l", col = "red")
library(scales)
points( data.len.age$length ~ data.len.age$age,
     col = alpha(blues9, 0.4), pch=16, cex = 1) 



##### 
### calculate the 95 % confid. envelope 
#####

# make a  matrix for  L(t) and t, with 'nruns' columns and 't.seq' rows.

L.mat <- matrix(data = NA, nrow = (length(t.seq.p)), ncol = (length(res$K)))

# fill the matrix with VBGF curves (length data for each age 't') from bootstrapping

L.mat[1,2] <- 0

# L.mat[r.i,c.w] <- 0

Lt.9 <- VBGF(list(Linf=res$Linf[9], 
                  K= res$K[9], 
                  t0= res$t0[9]),
             t=t.seq.p)

L.mat[,9] <- round(Lt.9,2)

# fill the matrix with data (curve swarm length data, one curve per column)

for (w in 1:(ncol(L.mat))) {
  
  Lt <- VBGF(list(Linf=res$Linf[w], 
                  K= res$K[w], 
                  t0= res$t0[w]),
             t=t.seq.p)
  
  
  L.mat[,w] <- round(Lt,2)
  
}

# calculate 95% quantiles for each age 't'

L.t.i <- L.mat[20,]
lo.ci.i <-  quantile(L.t.i, probs = c(0.025))
up.ci.i <- quantile(L.t.i, probs = c(0.975))

lo.ci.vec <- as.numeric(apply(L.mat, 1, quantile, probs = 0.025))
up.ci.vec <- as.numeric(apply(L.mat, 1, quantile, probs = 0.975))

# insert CIs into the plot

plot(t.seq2, L.seq.2, t="l", col = "red", xlim = c(0,max(age)), ylim = c(0,max(length)*1.2),
     xlab = "age (y)", ylab = "L (cm)")
points(data.len.age$length ~ data.len.age$age, col = "navyblue")

lines(t.seq2, lo.ci.vec, lty = 3, col = "navyblue", lwd = 1.5)
lines(t.seq2, up.ci.vec, lty = 3, col = "navyblue", lwd = 1.5)


# Kimura plot ---------------------------

plot(res$K~res$Linf , ylim = c(0,0.5), xlim = c(30,100), col = alpha("blue", 0.2))

#fishboot style

library(fishboot)

res2 <- list(bootRaw = res)
res2$bootRaw$Linf

LinfK_scatterhist(res2, xlim = c(50,80))

#univariate_density(res2)


######
# plot LFQ MA 5 results, nboot = 400 ###########



res <-  EXP_E_MA5_400

#plot curve

res_SA

res_MA5_OK <-  res_SA

res_MA5_OK$bootRaw <- res

res_MA5_OK$par$Linf <- 70.46

res_MA5_OK$par$K <- 0.1531

res_MA5_OK$par$t_anchor <- 0.47

res_MA5_OK$par$C <- 0.5270

res_MA5_OK$par$ts <- 0.50096

res_MA5_OK$par$phiL <- 2.902 


par <- opar

library(TropFishR)
plot(res_MA5_OK, draw = TRUE)
lfqFitCurves(res_MA5_OK, col="blue", par=res_MA5_OK$par, draw=TRUE)

res_MA5_OK$par




t.seq2 <- seq(0,max(age) , by = 0.1)
Linf.2 <-  round((output3$Linf),2)
K.2 <- round((output3$K),3)
t0.2 <- output3$t0
L.seq.2 <- VBGF(t = t.seq2, list(Linf=Linf.2, K=K.2, C= 0, t0 = t0.2))  # t.zero = 0 ???


plot(t.seq2, L.seq.2, t="l", col = "blue", xlim = c(0,max(age)), ylim = c(0,max(length)*1.2),
     xlab = "age (y)", ylab = "L (cm)")
points(data.len.age$length ~ data.len.age$age, col = "navyblue")


###
# plot with grey curve swarms
###


plot(t.seq2, L.seq.2, t="l", col = "green", xlim = c(0,max(age)), ylim = c(0,max(length)*1.2),
     xlab = "age (y)", ylab = "L (cm)")
points(data.len.age$length ~ data.len.age$age, col = "navyblue")


# library(TropFishR)
# calculate Lt from t, with t0
# t.seq.p = seq(0,max(agerec.2),0.1)
# Lt.p <- VBGF(list(Linf=res$bootRaw$Linf[1], 
#                   K= res$bootRaw$K[1], 
#                   t0= res$bootRaw$t_zero_i[1]),
#              t=t.seq.p)
# lines(t.seq.p, Lt.p, col = rgb(0, 0, 0, 0.25), lwd = 0.5 ) # transparent lines

res <- EXP_Otol_len_age_400

library(scales)

for (w in  1:length(res$K)) {
  
  t.seq.p = seq(0,max(data.len.age$age),0.1)
  Lt.p <- VBGF(list(Linf=res$Linf[w], 
                    K= res$K[w], 
                    t0= res$t0[w]),
               t=t.seq.p)
  # lines(t.seq.p, Lt.p, col = rgb(0, 0, 0, 0.2), lwd = 0.3 ) # transparent lines
  
  lines(t.seq.p, Lt.p, col = alpha("grey", 0.1), lwd = 0.3 ) # transparent lines
  
}

#points(data.len.age$length ~ data.len.age$age, col = "navyblue")
points(data.len.age$length ~ data.len.age$age, col = "green")

### nice plot

plot(t.seq2, L.seq.2, t="l", col = "green", xlim = c(0,max(age)), ylim = c(0,max(length)*1.2),
     xlab = "age (y)", ylab = "L (cm)")
lines(t.seq2, L.seq.2, t="l", col = "red")
library(scales)
points( data.len.age$length ~ data.len.age$age,
        col = alpha(blues9, 0.4), pch=16, cex = 1) 

dev.off()

# reset_par() # resets params ... all in one window



##### 
### calculate the 95 % confid. envelope 
#####

# make a  matrix for  L(t) and t, with 'nruns' columns and 't.seq' rows.

L.mat <- matrix(data = NA, nrow = (length(t.seq.p)), ncol = (length(res$K)))

# fill the matrix with VBGF curves (length data for each age 't') from bootstrapping

L.mat[1,2] <- 0

# L.mat[r.i,c.w] <- 0

Lt.9 <- VBGF(list(Linf=res$Linf[9], 
                  K= res$K[9], 
                  t0= res$t0[9]),
             t=t.seq.p)

L.mat[,9] <- round(Lt.9,2)

# fill the matrix with data (curve swarm length data, one curve per column)

for (w in 1:(ncol(L.mat))) {
  
  Lt <- VBGF(list(Linf=res$Linf[w], 
                  K= res$K[w], 
                  t0= res$t0[w]),
             t=t.seq.p)
  
  
  L.mat[,w] <- round(Lt,2)
  
}

# calculate 95% quantiles for each age 't'

L.t.i <- L.mat[20,]
lo.ci.i <-  quantile(L.t.i, probs = c(0.025))
up.ci.i <- quantile(L.t.i, probs = c(0.975))

lo.ci.vec <- as.numeric(apply(L.mat, 1, quantile, probs = 0.025))
up.ci.vec <- as.numeric(apply(L.mat, 1, quantile, probs = 0.975))

# insert CIs into the plot

plot(t.seq2, L.seq.2, t="l", col = "red", xlim = c(0,max(age)), ylim = c(0,max(length)*1.2),
     xlab = "age (y)", ylab = "L (cm)")
points(data.len.age$length ~ data.len.age$age, col = "navyblue")

lines(t.seq2, lo.ci.vec, lty = 3, col = "navyblue", lwd = 1.5)
lines(t.seq2, up.ci.vec, lty = 3, col = "navyblue", lwd = 1.5)



# Kimura plot

plot(res$K~res$Linf , ylim = c(0,0.5), xlim = c(30,100), col = alpha("blue", 0.2))

#fishboot style

library(fishboot)

res2 <- list(bootRaw = res)
res2$bootRaw$Linf

LinfK_scatterhist(res2, xlim = c(50,80))

#univariate_density(res2)


### 6. Tests and comparisons --------------------------------------------------


# 6.1 CI width   ------------

# test used to compare 95%CI widths: interquantile range test

# Source: 
# The interquantile range test, V01
# Ralf Schwamborn
# 25 de July 2019

# https://rpubs.com/rschwamborn/515419

# Defines the function 'interquant_r.test'  # ------
# tests for differences between 95% Inter-Quantile Range of standardized distributions    
# performs a non-parametric Harrell-Davis quantile test
# uses the function Qanova within the R package WRS2 (Mair and Wilcox, 2017).

library(WRS2)

interquant_r.test <- function(dat.A, dat.B, n.boot) {
  
  
  dat.A.stand <-   dat.A - quantile(dat.A, 0.025 )                             
  
  dat.B.stand <-   dat.B - quantile(dat.B, 0.025 )                             
  
  
  dat.frameAB <- data.frame(dat.A.stand, dat.B.stand)
  dat.frameAB.stack <-  stack(dat.frameAB)
  
  
  res <- Qanova(values ~ ind, dat.frameAB.stack, q =  0.975, nboot = n.boot)
  
  res.p.value <-   res$p.value[1,2]
  
  
  if (res.p.value == 0)
    
    paste("p-value < 0.0001")
  
  else
    
    paste ("p-value = ", res.p.value )
  
}

# Example: three datsets 

A <- rnorm (2000,mean = 1, sd = 5 ) # standard deviation = 5
B <- rnorm (2000,mean = 3000, sd = 5.02 ) # standard deviation =  basically the same as A
C <- rnorm (2000,mean = 1000, sd = 50 ) # much larger standard deviation

interquant_r.test(A,C, 500) # performs the interquant_r.test between "A" and "B", with 500 runs 

## [1] "p-value =  0.104"

# the large "p" value (p> 0.05) means that 95% interquantile ranges are not signifficantly different between datasets "A' and "B".

interquant_r.test(A,C, 500) # performs the interquant_r.test between "A" and "C", with 500 runs 

## [1] "p-value < 0.0001"

# the small "p" value (p < 0.05) means that 95% interquantile ranges are si

### Apply the test to the present data ----------

## Otoliths vs MA 5 -----------------

# Linf
interquant_r.test(EXP_E_MA5_400$Linf,EXP_Otol_len_age_400$Linf, 500) # performs the interquant_r.test between "A" and "C", with 500 runs 
# "p-value < 0.0001"

# K
interquant_r.test(EXP_E_MA5_400$K,EXP_Otol_len_age_400$K, 500) # performs the interquant_r.test between "A" and "C", with 500 runs 
# "p-value < 0.0001"


# Phi'
#interquant_r.test(EXP_E_MA5_400$phiL,EXP_Otol_len_age_400$Phi_vec, 500) # performs the interquant_r.test between "A" and "C", with 500 runs 


EXP_Otol_len_age_400[392,] # contains one -Inf value

EXP_Otol_len_age_400.b <- EXP_Otol_len_age_400

EXP_Otol_len_age_400.b[- grep("-Inf", EXP_Otol_len_age_400.b$Phi_vec),]

EXP_Otol_len_age_400.b[ grep("-Inf", EXP_Otol_len_age_400.b$Phi_vec, invert = TRUE) , ]

length(EXP_Otol_len_age_400.b$K)
length(EXP_Otol_len_age_400$K)

EXP_Otol_len_age_400.b[392,]  # contains a  "-Inf" value  in  row 392

EXP_Otol_len_age_400.b = EXP_Otol_len_age_400.b[-392,] 
# delete the bad row (otoliths)


EXP_E_MA5_400.b <- EXP_E_MA5_400 # remove one row at random (LFA)
rowx<- (round( runif(n=1,0,400))) # remove one row at random
EXP_E_MA5_400.b= EXP_E_MA5_400.b[-rowx,] # remove one row at random


# Phi'
interquant_r.test(EXP_E_MA5_400.b$phiL,EXP_Otol_len_age_400.b$Phi_vec, 500) # performs the interquant_r.test between "A" and "C", with 500 runs 
"p-value < 0.0001"


# Compare MA widths ------------------

# compare MA = 5 , 7, 9 11, 13, (pcrit with Bonferroni correction!)
# Multiple comparisons: use pcrit with Bonferroni correction! 
# pBonf = pcrit / B   B = number of comparisons 
# pBonf = pcrit / B

# number of data sets : 5
# number of possible comparisons : 10
# pBonf = pcrit / B   B = number of comparisons 
 pBonf = 0.05 / 10 #   B = number of comparisons 
pBonf # 0.005 pcrit is   0.005


#  MA 7
# output_A_MA7_100a.csv <-read.csv("output_A_MA7_100a.csv")
# output_A_MA7_200b.csv <-read.csv("output_A_MA7_200b.csv")
# output_A_MA7_300c.csv <-read.csv("output_A_MA7_300c.csv")
# #
# EXP_A_MA7 <- rbind(output_A_MA7_200b.csv, output_A_MA7_300c.csv)
# EXP_A_MA7_500 <- rbind(output_A_MA7_200b.csv, output_A_MA7_300c.csv)
# EXP_A_MA7_600 <- rbind(output_A_MA7_100a.csv, output_A_MA7_200b.csv, output_A_MA7_300c.csv)
# EXP_A_MA7_400 <- EXP_A_MA7_600[1:400,]
# 
# # head(EXP_A_MA7)
# #
#
#  View(EXP_A_MA7_400)
#
#  MA 9
#
#  output_B_MA9_200f <-read.csv("output_B_MA9_200f.csv")
# output_B_MA9_400e <-read.csv("output_B_MA9_400e.csv")
# 
# EXP_B_MA9 <- rbind(output_B_MA9_200f, output_B_MA9_400e)
# EXP_B_MA9_600 <- EXP_B_MA9
# EXP_B_MA9_400 <- EXP_B_MA9_600[1:400,]


#head(EXP_B_MA9)
#
# View(EXP_B_MA9_600)

# MA 11
# 
# output_C_MA11_200h <-read.csv("output_C_MA11_200h.csv")
# output_C_MA11_400g <-read.csv("output_C_MA11_400g.csv")
# 
# EXP_C_MA11 <- rbind(output_C_MA11_200h, output_C_MA11_400g)
# EXP_C_MA11_600 <- EXP_C_MA11
# EXP_C_MA11_400 <- EXP_C_MA11_600[1:400,]


# head(EXP_C_MA11_400)
# head(EXP_C_MA11_600)

# View(EXP_C_MA11_400)


# MA 13
# #
# output_D_MA13_200i <-read.csv("output_D_MA13_200i.csv")
# output_D_MA13_400j <-read.csv("output_D_MA13_400j.csv")
# 
# EXP_D_MA13 <- rbind(output_D_MA13_200i, output_D_MA13_400j)
# EXP_D_MA13_600 <- EXP_D_MA13
# EXP_D_MA13_400 <- EXP_D_MA13_600[1:400,]
# 
# 
# head(EXP_D_MA13_400)
# head(EXP_D_MA13_600)

# View(EXP_D_MA13_600)


# MA 5

# 
# output_E_MA5_200k <-read.csv("output_E_MA5_200k.csv")
# output_E_MA5_200l <-read.csv("output_E_MA5_200l.csv")
# 
# EXP_E_MA5_400 <- rbind(output_E_MA5_200k, output_E_MA5_200l)


head(EXP_E_MA5_400)

#View(EXP_E_MA5_400)

summary(EXP_E_MA5_400)


# 
# # MA5 vs MA13
# 
# # Linf
# interquant_r.test(EXP_E_MA5_400$Linf,EXP_D_MA13_400$Linf, 1000) # performs the interquant_r.test between "A" and "C", with 500 runs 
# # "p-value < 0.0001"
# 
# # K
# interquant_r.test(EXP_E_MA5_400$K,EXP_D_MA13_400$K, 1000) # performs the interquant_r.test between "A" and "C", with 500 runs 
# # "p-value < 0.0001"
# 
# 
# # Phi'
# interquant_r.test(EXP_E_MA5_400$phiL,EXP_D_MA13_400$phiL, 1000) # performs the interquant_r.test between "A" and "C", with 500 runs 
# # "p-value < 0.0001"
# 
# 
# # 
# # MA5 vs MA9
# 
# # Linf
# interquant_r.test(EXP_E_MA5_400$Linf,EXP_B_MA9_400$Linf, 1000) # performs the interquant_r.test between "A" and "C", with 500 runs 
# # "p-value < 0.0001"
# 
# # K
# interquant_r.test(EXP_E_MA5_400$K,EXP_B_MA9_400$K, 1000) # performs the interquant_r.test between "A" and "C", with 500 runs 
# # "p-value < 0.0001"
# 
# 
# # Phi'
# interquant_r.test(EXP_E_MA5_400$phiL,EXP_B_MA9_400$phiL, 1000) # performs the interquant_r.test between "A" and "C", with 500 runs 
# # ""p-value < 0.0001"
# 
# 
# 
# # MA5 vs MA7
# 
# # Linf
# interquant_r.test(EXP_E_MA5_400$Linf,EXP_A_MA7_400$Linf, 1000) # performs the interquant_r.test between "A" and "C", with 500 runs 
# # "p-value =  0.352" not sign.
# 
# # K
# interquant_r.test(EXP_E_MA5_400$K,EXP_A_MA7_400$K, 1000) # performs the interquant_r.test between "A" and "C", with 500 runs 
# # "p-value =  0.008" not sign
# 
# 
# # Phi'
# interquant_r.test(EXP_E_MA5_400$phiL,EXP_A_MA7_400$phiL, 1000) # performs the interquant_r.test between "A" and "C", with 500 runs 
# # "p-value =  0.272" n.s.
# 
# 
# 
# 
# 
# 
# # 6.2 Two-sample Median test ---------------
# 
# # Comparison of median values of growth parameter estimates
# # (K, Linf and Phi')
# 
# # test used:
# 
# #### The Bootstrapped two-sample-test  -------------------
# # Version 1.2 beta
# ## R. Schwamborn, March 2022
# 
# # Source: 
# 
# # https://rpubs.com/rschwamborn/880212
# 
# 
# ##  Short description ------------
# 
# #  Calculates a "p" value for the difference in medians. 
# #  The "Bootstrapped two-sample-test" gives a non-parametric two-sided "p" value (we do not assume normal distribution), based on the bootstrap posteriors for the estimate (median).
# #  Inputs:  Two bootstrap posteriors from  samples: "A" and "B" (see function below).
# 
# ## Function --------------
# 
# #   boot.p.two.sample.diff.posteriors -------------------
# 
# # Provides a two-sided "p" value for the difference in mean (or medians) for two samples 
# # Important: smaller sample  first
# # Uses two  posterior distributions (for means or for medians)  


boot.p.two.sample.diff.posteriors <- function (boot.post.A, boot.post.B) {
  
  # difference between posteriors
  
  boot.statisticsA <- boot.post.A # posterior for sample A
  boot.statisticsB <- boot.post.B # posterior for sample B
  
  
  boot.statisticsDiff <- boot.statisticsB - boot.statisticsA
  
  
  # calculate "p" value  from bootstrap diff.  
  
  N.below.zero <- length(boot.statisticsDiff [boot.statisticsDiff< 0]) # 7 are below 0 
  
  N.total <- length(boot.statisticsDiff ) # all 
  
  ( one.sided.p_value <- N.below.zero / N.total) # one-sided "p" value
  
  (two.sided.p_value <- 2* (N.below.zero / N.total)) # two-sided "p" value
  
  c("two.sided.p_value", two.sided.p_value)  # output
  
}

# 
# # K
# boot.p.two.sample.diff.posteriors(EXP_E_MA5_400$K, 
#                                   EXP_Otol_len_age_400$K
#                                   ) # 
# # "two.sided.p_value" "0.425"  # not signifficant! 
# 
# # Linf
# boot.p.two.sample.diff.posteriors(EXP_Otol_len_age_400$Linf,
#                                   EXP_E_MA5_400$Linf) # 
# # "two.sided.p_value" "0.215"  # not signifficant!
# 
# # Phi'
# boot.p.two.sample.diff.posteriors(EXP_Otol_len_age_400.b$Phi_vec,
#                                   EXP_E_MA5_400.b$phiL) # 
# 
# # "two.sided.p_value" "0.626"  # not signifficant!
# 
# 

# 6.3 diverse other tests

### permutation tests for "p" ####
####
### permutate the length-at-age-fit of the  VBGF #####
# source: https://stats.stackexchange.com/questions/316483/manually-bootstrapping-linear-regression-in-r
### problem: bad, noisy data  give many "zero" values... 
# possible solution: note error as "zero"...




attach(Lsynagris_lenage)
plot(length ~ age)


## adjust a VBGF function

# non linear least squares method , TropFishR::growth_length_age,method = "LSM"

library(TropFishR)

output2 <- growth_length_age(param = Lsynagris_lenage, method = "LSM",Linf_init = 30, CI = TRUE, age_plot=NULL)

summary(output2$mod)

# ?growth_length_age

output3 <- growth_length_age(param = Lsynagris_lenage, method = "LSM",Linf_init = (max(length)), CI = TRUE, age_plot=NULL)

summary(output3$mod)

output3$Linf
output3$K
output3$t0

# ?growth_length_age

## now permutate the VBGF with length-at-age

data.len.age <- data.frame(Lsynagris_lenage)

P = 2000 ## number of permutations
res = data.frame(Linf = numeric(P) , K = numeric(P), t0 = numeric(P)) ## vector to hold results
n = length(data.len.age$length) # number of data pairs


attach(data.len.age)
#plot(length ~ age) # plot background points

for(b in 1:P){
  
  
  tryCatch({   
  
  seed <- round(as.numeric(Sys.time())+runif(1, min = 0, max = 1e4)) # maximize stochasticity
  set.seed(seed)
  
  i = sample(x = 1:n, size = n, replace = FALSE) ## sample indices
  
  bootsam.dataxy <- data.len.age[i,] ## get data
  bootsam.dataxy <- list(age=  round(bootsam.dataxy$age,0) , length= data.len.age$length)
  
  
  output <- growth_length_age(param = bootsam.dataxy, method = "LSM",
                              Linf_init = (max(bootsam.dataxy$length)), CI = FALSE, age_plot=NULL)
  
  output$Linf
  
  ## store results
  res$Linf[b]  <-   output$Linf 
  res$K[b]  <- output$K  
  res$t0[b]<- output$t0
  
  }, error=function(e){})
  
}



# View(res)

hist(res$Linf)
summary(res$Linf)
quantile(res$Linf, probs = c(0.025,0.5,0.975))

hist(res$K)
summary(res$K)
quantile(res$K, probs = c(0.025,0.5,0.975))

hist(res$t0)
summary(res$t0)
quantile(res$t0, probs = c(0.025,0.5,0.975))

######
# plot results ###########

attach(data.len.age)

output3$Linf
output3$K
output3$t0

# hist(res$Linf)
# hist(res$Linf, xlim = c((min(res$Linf)),(output3$Linf*1) ))
quantile(res$Linf, probs = c(0.95))
pcrit.Linf <- quantile(res$Linf, probs = c(0.95))
abline(v = pcrit.Linf, col = "red")
abline (v = output3$Linf, col = "blue")


###
# calculate "p" from absolute "n" of observations for K, Linf, t0, etc. #####

# calculate "p" for Linf
#View (res)

Linf.obs <- output3$Linf  

n.larger <- length(res$Linf[res$Linf> Linf.obs])
n.smaller <- length(res$Linf[res$Linf< Linf.obs])
n.Perm <- P

p.one.sided.Linf <- n.larger /n.Perm # 0.02 -> OK
p.one.sided.Linf

# now calculate "p" for K 
# View (res)

K.obs <- output3$K

n.larger <- length(res$K[res$K> K.obs])
n.smaller <- length(res$K[res$K< K.obs])
n.Perm <- P

p.one.sided.K <- n.larger /n.Perm # p = 0.02 -> OK
p.one.sided.K


# now calculate "p" for t_zero 
# View (res)

t0.obs <- output3$t0

if (t0.obs > 0) {
n.larger <- length(res$t0[ (res$t0)  > (t0.obs)])
n.smaller <- length(res$t0[  (res$t0)  < (t0.obs) ])
n.Perm <- P
p.one.sided.t0 <- n.larger /n.Perm # p = 0.02 -> OK
p.one.sided.t0 }

if (t0.obs < 0)  { n.larger <- length(res$t0[ (res$t0)  < (t0.obs)])
n.smaller <- length(res$t0[  (res$t0)  > (t0.obs) ])
n.Perm <- P
p.one.sided.t0 <- n.larger /n.Perm 
p.one.sided.t0 }

p.one.sided.t0 # p = 0.02 -> OK




###  7. FishBase serach and comparisons ---------------------------------------

# Lutjanidae growth - from Fishbase ----------------

### rfishbase search --------------------
##  113 Lutjanidae species in FishBase
##  57 Lutjanidae species with growth data (317 data)
##  
##  30 useful growth for species "Lutjanus synagris"

# Species examples  from rfishbase -----------------

# source: https://cran.r-project.org/web/packages/rfishbase/readme/README.html

# remotes::install_github("ropensci/rfishbase")
# 
# library("rfishbase")
# library("dplyr") # convenient but not required


library(rfishbase)

# Fast  growing species  ---------------

# Fast  growth: K above 0.2
# Very fast  growth: K above 0.4


# example of very fast growth: Tilapia


P_char <- popchar("Oreochromis niloticus") 

P_char$Lmax

pop1 <- popgrowth("Oreochromis niloticus",
                  server = NULL) # server = getOption("FISHBASE_API", "fishbase")



median(pop1$K) # median for Tilapia is K = 2.9 # very fast growth!
quantile(pop1$K, c(0.025, 0.975))
# 95% of the K values for Tilapa  in FishBase 
# are between K = 0.21 and K = 13 (fast growth) # fast growth: above 0.2 

# Phi_p = log10(K) + 2 log10(Linf)
Phi_vec.pop1 = log10(pop1$K) + 2 * log10(pop1$Loo)
quantile(Phi_vec.pop1, c(0.025, 0.975))
# 95% of the Phi' values  in FishBase 
#Phi' are between Phi' = 2.25 and Phi'= 3.58

## Lutjanidae (all) --------------


famlist <- species_list(Family = "Lutjanidae")

famlist # all Lutjanidae species in fishbase (113 species)

popLutfam <- popgrowth(famlist,
                       server = NULL) # server = getOption("FISHBASE_API", "fishbase")


popLutfam$K

median(popLutfam$K, na.rm = TRUE) # 0.22 = median for Lutjanidae is K = 0.22 # moderate growth!
# 0.22 = "moderate growth"
quantile(popLutfam$K, c(0.025, 0.975), na.rm = TRUE)
# 95% of the K values 
# are between K = 0.0886 and K = 0.82 (fast growth) # fast growth: above 0.2 
# 2.5%  97.5% 
#   0.0886 0.8212  for Lutjanidae  in FishBase
# Wide range for K, from very slow to very fast! 

# Phi_p = log10(K) + 2 log10(Linf)

Phi_vec.popLutfam = log10(popLutfam$K) + 2 * log10(popLutfam$Loo)
median(Phi_vec.popLutfam, na.rm = TRUE) # median = 2.87 (moderate perf.)
quantile(Phi_vec.popLutfam, c(0.025, 0.975), na.rm = TRUE)
# 95% of the Phi' values  in FishBase 
#Phi' are between Phi' = 2.135 (low perf.) and Phi'= 3.404 (high perf.)


popLutfam$Loo

median(popLutfam$Loo, na.rm = TRUE) # 0.22 = median for Lutjanidae is Loo = 0.22 # moderate growth!
# 60  = "median  Linf for Lutjanids"
quantile(popLutfam$Loo, c(0.025, 0.975), na.rm = TRUE)
# 95% of the Loo values 
# are between Loo = 21.7 and Loo = 109.8 
# 2.5%  97.5% 
#   21.65 109.80   for Lutjanidae  in FishBase
# Wide range for Loo, from very small  to very large 

# How many useful data:
length(popLutfam$Loo) # 373 data 
d <- popLutfam$Loo
d <- d[!is.na(d)]
d
length(d) # 317 useful growth for Lutjanidae  in FishBase  (113 species) 



# L. synagris only (Ariaco) -------------

# new data (28 feb) - directly from FishBase (web), with conversion for FL and Sl, to TL  

#Lsyn_FB_WEB <- read.csv("C:/Users/RALF/Desktop/Ralf_diversos 2018_19_iii/Papers/00 - Ale - Paper Otolitos vs LFA Lutjanus/Lsynagris_FishBaseOKII.csv")


urlfile5 <- 'https://raw.githubusercontent.com/rschwamborn/L_synagris_growth_otoliths_and_LFA/main/Lsynagris_FishBaseOKII.csv'

Lsyn_FB_WEB <- read.csv(url(urlfile5))

head(Lsyn_FB_WEB)
summary(Lsyn_FB_WEB)



# View(Lsyn_FB_WEB) # 28 data sets


attach(Lsyn_FB_WEB)

median(Lsyn_FB_WEB$K) # median K for Lane snapper  is 0.225
quantile(Lsyn_FB_WEB$K, c(0.025, 0.975))# 0.0955  to 0.3646 
summary(Lsyn_FB_WEB$K)


median(Lsyn_FB_WEB$Linf_OK) # median for Lane snapper  is  50.18 cm
quantile(Lsyn_FB_WEB$Linf_OK, c(0.025, 0.975)) # 32.81 to 70.86
summary(Lsyn_FB_WEB$Linf_OK)


median(Lsyn_FB_WEB$Phi_p) # median for Lane snapper  is  50.18 cm
quantile(Lsyn_FB_WEB$Phi_p, c(0.025, 0.975)) # 32.81 to 70.86
summary(Lsyn_FB_WEB$Phi_p)



# Compare  FishBase-web and rfishbase ---------------------

# FishBase web - K - L. synagris only
sort(Lsyn_FB_WEB$K)


# Lutjanus synagris,  Lane snapper (Ariaco)
popLs <- popgrowth("Lutjanus synagris",
                   server = NULL) # server = getOption("FISHBASE_API", "fishbase")

popLs$Species

popLs$K

# rfishbase  - K - L. synagris only
sort(popLs$K)



# How many useful data:
length(popLs$Loo) # 30 data  for Lutjanus synagris
d <- popLs$Loo
d <- d[!is.na(d)]
d
length(d) # 30 useful growth for Lutjanus synagris,  Lane snapper  in FishBase  (1 species) 

median(popLs$K) # median for Lane snapperntile(popLs$K, c(0.025, 0.975))
# 2.5%    97.5% 
#   0.097025 0.432125 
# 95% of the K values for Lane snapper  in FishBase 
# are between K = 0.097 (very slow growth) and K = 0.432 (very fast growth) 
median( popLs$K) #  median = 0.2305, moderate growth
# 0.2305

# Phi_p = log10(K) + 2 log10(Linf)
Phi_vec.popLs = log10(popLs$K) + 2 * log10(popLs$Loo)
quantile(Phi_vec.popLs, c(0.025, 0.975))
# 95% of the Phi' values  in FishBase 
#Phi' are between Phi' = 2.218758 (low perf.) and Phi' = 3.043 (moderate perf., , close to 2.8) 
median(Phi_vec.popLs) # 2.649 (moderate perf., close to 2.8)

median(popLs$Loo) # median for Lane snapper  is K = 2.9 # very fast growth!
quantile(popLs$Loo, c(0.025, 0.975))
# 2.5%    97.5% 
#  31.0625 71.00 
# 95% of the Linf values for Lane snapper  in FishBase 
# are between K = 31 (very small) and Linf = 71 (very large Linf) 
median( popLs$Loo) #  median = 43.8 cm, small Linf
# 0.2305


P_char <- popchar("Lutjanus synagris") 

P_char$Lmax # Lmax in FishBase ranges between 18.2 cm and 43 cm

Lmax_lane_snapper <- c(41.0, 36.0 , 38.0 ,43.0 , 18.2)
median(Lmax_lane_snapper)
# 38 cm = median Lmax !
quantile(Lmax_lane_snapper, c(0.025, 0.975))
# 2.5%    97.5% 
#  19.98 42.80 



# All fishbase growth data -----------------


#    fb_tbl("popgrowth") # 11,861 growth data 

#### Plot Figure for paper ---------------
# plot growth data ------------

plot(popLutfam$K ~popLutfam$Loo,  xlim = c(0,170), ylim = c(0,1.4),
     col = "lightgrey", pch = 18, cex = 0.8, xaxt="n", yaxt="n")
points(popLs$K ~popLs$Loo, col = "blue", type = "p" , pch = 16)
points(popLs$K ~popLs$Loo, col = "darkblue", type = "p" , pch = 1)



# Changing x axis ticks
xtick<-seq(0, 170, by=10)
axis(side=1, at=xtick, labels = TRUE)
# text(x=xtick,  par("usr")[3], 
#     labels = xtick, srt = 45, pos = 1, xpd = TRUE)
# Changing y axis
ytick<-seq(0, 1.4, by=0.2)
axis(side=2, at=ytick, labels = TRUE)
#text(par("usr")[1], ytick,  
#     labels = ytick, srt = 45, pos = 2, xpd = TRUE)


# fishboot style with Phi' isopleths------------
library(fishboot)

#setwd("C:/Users/RALF/Desktop/Ralf_diversos 2018_19_iii/Papers/00 - Ale - Paper Otolitos vs LFA Lutjanus/Data_inputs_outputs")
# EXP_Otol_len_age_400 <-read.csv("EXP_Otol_len_age_400.csv")


urlfile6 <- 'https://raw.githubusercontent.com/rschwamborn/L_synagris_growth_otoliths_and_LFA/main/EXP_Otol_len_age_400.csv'

EXP_Otol_len_age_400 <- read.csv(url(urlfile6))

head(EXP_Otol_len_age_400)
summary(EXP_Otol_len_age_400)



res <- EXP_Otol_len_age_400
res2 <- list(bootRaw = res)
res2$bootRaw$Linf


# contours and isopleths only for plot - otoliths
LinfK_scatterhist(res2, pt.col = "white",
                  probs = c(95),
                  xlim = c(0,170), ylim = c(0,1.4) ,
                  phi.contour = TRUE, phi.contour.col = "lightgrey")

# contours and isopleths only for plot - LFA - MA5
#EXP_E_MA5_400 <-read.csv("EXP_E_MA5_400.csv")
#EXP_E_MA5_400 <- read.csv(url("http://some.where.net/data/foo.csv"))

summary(EXP_E_MA5_400)

res <-  EXP_E_MA5_400

res2 <- list(bootRaw = res)
res2$bootRaw$Linf

# contours and isopleths only for plot - otoliths
LinfK_scatterhist(res2, pt.col = "white",
                  probs = c(95),
                  xlim = c(0,170), ylim = c(0,1.4) ,
                  phi.contour = TRUE, phi.contour.col = "lightgrey")










