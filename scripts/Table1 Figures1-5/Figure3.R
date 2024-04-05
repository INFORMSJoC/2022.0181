rm(list = ls())
setwd("d:/Downloads/2022.0181/scripts/Table1 Figures1-5")
# Figure 3
# read results
path_result <- "Results-Sim22_5_1_context/"
method_list <- c("ROP", "ROPr")
# varying voutsd
set = "settingalpha"
voutsd_list = c(4,6,8,10,20)
RE_B = array(dim = c(2,length(voutsd_list),100)) 
RE_C  = array(dim = c(2,length(voutsd_list),100)) 
RE_time  = array(dim = c(2,length(voutsd_list),100)) 
i = 1
for (voutsd in voutsd_list) {
  print(voutsd)
  
  file_name_B <- paste(set,voutsd, "B.txt", sep = "_")
  file_path_B <- paste(path_result, file_name_B, sep = "")
  rt = read.table(file_path_B)
  RE_B[1,i,] = as.numeric(round(rt[,1], 3))
  RE_B[2,i,] = as.numeric(round(rt[,2], 3))
  
  file_name_C <- paste(set,voutsd, "C.txt", sep = "_")
  file_path_C <- paste(path_result, file_name_C, sep = "")
  rt = read.table(file_path_C)
  RE_C[1,i,] = as.numeric(round(rt[,1], 3))
  RE_C[2,i,] = as.numeric(round(rt[,2], 3))
  
  file_name_time <- paste(set,voutsd, "time.txt", sep = "_")
  file_path_time <- paste(path_result, file_name_time, sep = "")
  rt = read.table(file_path_time)
  RE_time[1,i,] = as.numeric(round(rt[,1], 3))
  RE_time[2,i,] = as.numeric(round(rt[,2], 3))
  
  i = i+1
}

##### draw plot   
library(ggplot2)
library(latex2exp)

da1=c(RE_B[2,1,],RE_B[2,2,],RE_B[2,3,],RE_B[2,4,],RE_B[2,5,])
noutl=factor(rep( voutsd_list,each=100))
data1 = data.frame(da1, noutl )
P1 = ggplot( data1,aes(x = noutl, y = da1)) +
  geom_boxplot(alpha=0.7) +
  coord_cartesian(ylim = c(0,0.1))+
  scale_y_continuous(name=TeX("$Err(\\hat{B})$"))+
  scale_x_discrete(name =expression(alpha))

da2=c(RE_C[2,1,],RE_C[2,2,],RE_C[2,3,],RE_C[2,4,],RE_C[2,5,])
nout2=factor(rep( voutsd_list,each=100))
data2 = data.frame(da2, nout2)
P2 = ggplot( data2,aes(x = nout2, y = da2)) +
  geom_boxplot(alpha=0.7) +
  coord_cartesian(ylim = c(0,0.2))+
  scale_y_continuous(name=TeX("$Err(\\hat{C})$"))+
  scale_x_discrete(name =expression(alpha))

############
library(patchwork)
P1+P2