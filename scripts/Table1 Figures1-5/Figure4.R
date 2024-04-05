rm(list = ls())
setwd("d:/Downloads/2022.0181/scripts/Table1 Figures1-5")
# Figure 4
path_result <- "Results-Sim22_5_1_context/"
# read results
method_list <- c("ROP", "ROPr")
# varying voutsd
set = "settingSpa"
sparsity_list = c(c(8,9,9),c(16,18,18),c(24,27,27),c(32,36,36))
RE_B = array(dim = c(2,length(sparsity_list),100)) 
RE_C  = array(dim = c(2,length(sparsity_list),100)) 
i = 1
for (i in 1:4) {
  print(i)
  
  file_name_B <- paste(set,i, "B.txt", sep = "_")
  file_path_B <- paste(path_result, file_name_B, sep = "")
  rt = read.table(file_path_B)
  RE_B[1,i,] = as.numeric(round(rt[,1], 3))
  RE_B[2,i,] = as.numeric(round(rt[,2], 3))
  
  file_name_C <- paste(set,i, "C.txt", sep = "_")
  file_path_C <- paste(path_result, file_name_C, sep = "")
  rt = read.table(file_path_C)
  RE_C[1,i,] = as.numeric(round(rt[,1], 3))
  RE_C[2,i,] = as.numeric(round(rt[,2], 3))
  
  i = i+1
}

##### draw plot   
library(ggplot2)
library(latex2exp)

da1=c(RE_B[2,1,],RE_B[2,2,],RE_B[2,3,],RE_B[2,4,])
noutl=factor(rep( c("(8,9,9)","(16,18,18)","(24,27,27)","(32,36,36)"),each=100))
data1 = data.frame(da1, noutl )
P1 = ggplot( data1,aes(x = factor(noutl,levels = unique(noutl)), y = da1)) +
  geom_boxplot(alpha=0.7) +
  coord_cartesian(ylim = c(0,0.1))+
  scale_y_continuous(name=TeX("$Err(\\hat{B})$"))+
  scale_x_discrete(name =  TeX("$s_u$"))

da2=c(RE_C[2,1,],RE_C[2,2,],RE_C[2,3,],RE_C[2,4,])
nout2=factor(rep( c("(8,9,9)","(16,18,18)","(24,27,27)","(32,36,36)"),each=100))
data2 = data.frame(da2, nout2)
P2 = ggplot( data2,aes(x = factor(nout2,levels = unique(nout2)), y = da2)) +
  geom_boxplot(alpha=0.7) +
  coord_cartesian(ylim = c(0,0.5))+
  scale_y_continuous(name=TeX("$Err(\\hat{C})$"))+
  scale_x_discrete(name =TeX("$s_u$"))

############
library(patchwork)
P1+P2

