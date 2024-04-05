rm(list = ls())
# Table 3
# prepossing  
library(readr)
# note the file is too big to upload in the github, please follow the instrcution in README.md file in "d:/Downloads/2022.0181/data/" to get the file.
User_info <- read_csv("d:/Downloads/2022.0181/data/top_category_new_data.csv") 
C = User_info[which(User_info$category_id ==4145813),] # use different Category ID: 4756105,  2355072, 4145813,3607361,982926
User_used = C[,c(1,2,3,5)]
User_used = cbind(User_used,-1)
head(User_used)

cal_loss_once = function(i) {  
  library(RhpcBLASctl)
  blas_set_num_threads(1)
  temp = unclass(as.POSIXlt(User_used$time[i]))
  if(temp$mon == 10){
    timesC  = (temp$mday-25)*24+temp$hour
  }else if(temp$mon == 11){
    timesC  = 144 + (temp$mday-1)*24+temp$hour
  }
  return(timesC)
}


# snowfall 
library(purrr) # for cross3
library(ncvreg) # for Mest computing
library(snowfall) # parallel programming
library(quadprog)
library(parallel)
sfInit(parallel = TRUE, cpus = detectCores() - 1)
sfExport("User_used")
parallel_output = snowfall::sfSapply(1:nrow(User_used), cal_loss_once) 
sfStop()
User_used[,5]=parallel_output



#####
colnames(User_used)= c("user_id" ,    "item_id"  ,   "category_id", "time",        "timestamp") 
User_used = User_used[,-4]
#write.csv(User_used,"usedinformation.csv")

#timestamps
timestamps = length(rle(sort(User_used$timestamp))$values)
User_update = User_used[-which(User_used$timestamp<0),]
timestamps = length(rle(sort(User_update$timestamp))$values)
timestamps = c(1:timestamps)

# information
length(rle(sort(User_update$item_id))$values) 
length(rle(sort(User_update$user_id))$values) 
nrow(User_update)