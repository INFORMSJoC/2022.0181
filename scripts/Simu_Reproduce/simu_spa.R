## data describe
# correspond to simulation data set with secure(X,B,E) + r4 outlier

############ simulatipn 

rm(list = ls())
setwd("~/scripts/Simu_Reproduce")
# load package
library(purrr) # for cross3
library(ncvreg) # for Mest computing
library(snowfall) # parallel programming
library(quadprog)
library(parallel)

##############    update to 22-4-17 (FARMA-1st) ################################
# initialization
path_result <- "Results-Sim22_5_1_context/"
if (!file_test("-d", path_result)) {
  dir.create(path_result)
}

# function of mask, swamp, detection
MSD = function(C_est,n,nout){
  mask = 0
  swamp = 0
  detection = 0
  for (i in 1:n) {
    if(i<=nout)
    {
      if( norm(C_est[i,],"2") == 0 )
      {mask = mask + 1}
    }else{
      if( norm(C_est[i,],"2") != 0 )
      {swamp = swamp + 1}
    }
  }
  
  if(mask == 0 & swamp == 0) detection = 1
  
  Mmask = mask/nout #the fraction of undetected outliers
  Mswamp = swamp/(n-nout) #the fraction of good points labeled as outliers
  Mdetection = detection #with no masking and no swamping.
  
  return(list(mask = Mmask, swamp = Mswamp, detection = Mdetection))
}





### step1 simulate 100 data set
runtimes = 100

cal_loss_once <- function(ii) 
{
  # varying alpha
  set = "settingSpa"
  n = 400
  q = 200
  nrank = 3
  rho_X = 0.5
  rho_E = 0.5
  nout = 5
  vout = NULL
  voutsd = 2
  nlev = 0
  vlev = NULL
  vlevsd = NULL
  #sparsity = c(8,9,9)
  snr = 0.75
  p = 500
  
  # load settings
  for (kk in 1:length(settings)) {
    if(ii < kk*100+1)
    {idx = kk
    break}
  }
  sparsity= settings[[idx]]
  
  
  #
  method_list <- c("ROP", "ROPr")
  
  # output
  B_methods <- array(0, dim = length(method_list))
  C_methods <- array(0, dim = length(method_list))
  BC_methods <- array(0, dim = length(method_list))
  rank_methods = array(0, dim = length(method_list))
  time_methods = array(0, dim = length(method_list))
  mask_methods = array(0, dim = length(method_list))
  swamp_methods = array(0, dim = length(method_list))
  detection_methods = array(0, dim = length(method_list))
  
  # output files
  file_name_B <- paste(set,idx, "B.txt", sep = "_")
  file_path_B <- paste(path_result, file_name_B, sep = "")
  
  file_name_C <- paste(set,idx, "C.txt", sep = "_")
  file_path_C <- paste(path_result, file_name_C, sep = "")
  
  file_name_BC <- paste(set,idx, "BC.txt", sep = "_")
  file_path_BC <- paste(path_result, file_name_BC, sep = "")
  
  file_name_rank <- paste(set,idx, "rank.txt", sep = "_")
  file_path_rank <- paste(path_result, file_name_rank, sep = "")
  
  file_name_time <- paste(set,idx, "time.txt", sep = "_")
  file_path_time <- paste(path_result, file_name_time, sep = "")
  
  file_name_mask <- paste(set,idx, "mask.txt", sep = "_")
  file_path_mask <- paste(path_result, file_name_mask, sep = "")
  
  file_name_swamp <- paste(set,idx, "swamp.txt", sep = "_")
  file_path_swamp <- paste(path_result, file_name_swamp, sep = "")
  
  file_name_detection <- paste(set,idx, "detection.txt", sep = "_")
  file_path_detection <- paste(path_result, file_name_detection, sep = "")
  
  # debug
  file_name_debug <- paste(set,idx, "debug.txt", sep = "_")
  file_path_debug <- paste(path_result, file_name_debug, sep = "")
  cat(ii, "\t", file = file_path_debug, append = TRUE)
  
  # library or source some useful packages and funcitons
  source("RS1.R")
  library("rrpack")
  source("setup_secure.R")
  source("r4.1.R")  # r4 something wrong with its original code
  source("RCGL.R")
  source("secure.R")
  
  # import data
  seed = 2022+ii
  simdata = setsecure(n = n,
                      p = p,
                      q = q,
                      s = sparsity,
                      nrank = nrank,
                      rho_X = rho_X,
                      rho_E = rho_E,
                      snr = snr,
                      nout = nout,
                      vout = vout,
                      voutsd = voutsd,
                      nlev = nlev,
                      vlev = vlev,
                      vlevsd = vlevsd) 
  
  Y <- simdata$Y
  X <- simdata$X
  B = simdata$B
  C = simdata$C
  
  trueEST = rbind(B,C)
  r = qr(trueEST)$rank
  Y.mean = X %*% B + C
  
  
  
  # ROP
  start1 = Sys.time()
  re1 = RS1(Y,X,r = (r+4),outlier = TRUE)
  stop1 = Sys.time()
  timecost1 = stop1 - start1
  rank.est1 = re1 $rank
  B.est1 = re1 $B
  C.est1 = re1 $C
  BC.est1 = re1 $coef
  temp1 = MSD(C.est1,n,nout)
  
  BC_methods[1] = norm(Y.mean - X %*% B.est1 - C.est1,"F")/norm(Y.mean,"F")
  B_methods[1] = norm(B - B.est1,"F")/norm(B,"F")
  C_methods[1] = norm(C - C.est1,"F")/norm(C,"F")
  rank_methods[1] = rank.est1
  time_methods[1] = timecost1
  mask_methods[1] = temp1$mask
  swamp_methods[1] = temp1$swamp
  detection_methods[1] = temp1$detection
  
  # ROPr  
  start2 = Sys.time()
  Y1 = Y - C.est1
  re2 = RS1(Y1,X,rank.est1,outlier = FALSE)
  stop2 = Sys.time()
  timecost2 = stop2 - start2 + timecost1
  B.est2 = re2$B
  rank.est2 = re2$rank
  
  BC_methods[2] = norm(Y.mean - X %*% B.est2 - C.est1,"F")/norm(Y.mean,"F")
  B_methods[2] = norm(B - B.est2,"F")/norm(B,"F")
  C_methods[2] = norm(C - C.est1,"F")/norm(C,"F")
  rank_methods[2] = rank.est2
  time_methods[2] = timecost2
  mask_methods[2] = temp1$mask
  swamp_methods[2] = temp1$swamp
  detection_methods[2] = temp1$detection
  
  
  # output results
  Sys.sleep(runif(1, min=0, max=1))
  cat(BC_methods, "\n", file = file_path_BC, append = TRUE)
  cat(B_methods, "\n", file = file_path_B, append = TRUE)
  cat(C_methods, "\n", file = file_path_C, append = TRUE)
  cat(rank_methods, "\n", file = file_path_rank, append = TRUE)
  cat(time_methods, "\n", file = file_path_time, append = TRUE)
  cat(mask_methods, "\n", file = file_path_mask, append = TRUE)
  cat(swamp_methods, "\n", file = file_path_swamp, append = TRUE)
  cat(detection_methods, "\n", file = file_path_detection, append = TRUE)
  
  # debug
  cat("done\t", ii, "\n", file = file_path_debug, append = TRUE)
  # end debug
}

# Setting: 
num_rep <- 100 # number of duplication, default:100
sparsity_list = list(c(8,9,9),c(16,18,18),c(24,27,27),c(32,36,36))
#p_list = c(500,800,1000,1200,1500,1800,2000)
settings <- sparsity_list


# snowfall
sfInit(parallel = TRUE, cpus = detectCores() - 1)
sfExport("path_result", "MSD","settings")
sfLapply(c(1:(length(settings) * num_rep)), cal_loss_once)
sfStop()

