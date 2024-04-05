################# Figure 6 
# 
rm(list = ls())
setwd("d:/Downloads/2022.0181/scripts/Table2 Figure6")
# load package
library(purrr) # for cross3
library(ncvreg) # for Mest computing
library(snowfall) # parallel programming
library(quadprog)
library(parallel)

X <- read.csv("d:/Downloads/2022.0181/data/file1used.csv")
Y <- read.csv("d:/Downloads/2022.0181/data/file2used.csv")

index = which(Y$Pathwayname %in% c("Plastoquinonebiosynthesis",
                                   "Phytosterolbiosynthesis",
                                   "Carotenoidbiosynthesis",
                                   "PorphyrinChlorophyllmetabolism"))
Y0 = Y[index ,]

Xused = t(X[,7:124]) 
Yused = t(Y0[,5:122])
n=nrow(Yused);q=ncol(Yused);p=ncol(Xused)
 
 
# ROP
source("RS2.R") 
start1 = Sys.time()
re1 = RS2(Yused,Xused,rmax = NULL,r = 5,outlier = TRUE)
stop1 = Sys.time()
timecost1 = stop1 - start1
rank.est1 = re1 $rank
B.est1 = re1 $B
C.est1 = re1 $C
BC.est1 = re1 $coef

order(apply(abs(C.est1), 1, sum),decreasing=TRUE)[1:10]
apply(abs(C.est1[order(apply(abs(C.est1), 1, sum),decreasing=TRUE)[1:10],]), 1, sum)


# ROPr
Y1 = Yused-C.est1
re2 = RS2(Y1,Xused,rmax = NULL,r = 5,outlier = FALSE)
B.est2 = re2$B
rank.est2 = re2$rank

 

# refit 
SVDY = svd(Y1,nu=3,nv=3)
factors = SVDY$u
values = SVDY$d[1:3]
values.stand = (values)/sqrt(q)
loadings = SVDY$v*sqrt(q)


path =c(rep("Carotenoid",11),rep("Phytosterol",25),rep("Plastoquinone",2),rep("Chlorophyll",24))

load = data.frame(sample = c(1:62),factor1 = loadings[,1],factor2 = loadings[,2],factor3 = loadings[,3],pathways = path)



xlab = c("Carotenoid","Phytosterol","Plastoquinone","Chlorophyll")

library(ggplot2)
P1 = ggplot(load, aes(x=sample, y=factor1, colour=pathways, shape = pathways)) + 
  geom_point()+ 
  geom_text(aes(label = sample, vjust = 1.1, hjust = -0.5, angle = 0), show_guide = FALSE,check_overlap = TRUE)+
  theme(axis.text.x = element_blank())+
  theme(panel.grid =element_blank())+
  coord_cartesian(ylim=c(-5,5))+
  geom_hline(aes(yintercept = 0),colour = "grey",linetype = "dashed")+
  geom_hline(aes(yintercept = 1),colour = "black",linetype = "dashed")+
  geom_hline(aes(yintercept = -1),colour = "black",linetype = "dashed")+
  geom_vline(aes(xintercept = 11),colour = "black",linetype = "dotted")+
  geom_vline(aes(xintercept = 36),colour = "black",linetype = "dotted")+
  geom_vline(aes(xintercept = 38),colour = "black",linetype = "dotted")+
  labs(x = "Genes from four pathways", y = "Factor coefficient",title = "First factor")+
  theme(legend.justification = c(1,0),legend.position = c(1,0))




P2 = ggplot(load, aes(x=sample, y=factor2, colour=pathways, shape = pathways)) + 
  geom_point()+ 
  geom_text(aes(label = sample,vjust = 0.5, hjust = -0.5,  angle = 0), show_guide = FALSE,check_overlap = TRUE)+
  theme(axis.text.x = element_blank())+
  theme(panel.grid =element_blank())+
  coord_cartesian(ylim=c(-5,5))+
  geom_hline(aes(yintercept = 0),colour = "grey",linetype = "dashed")+
  geom_hline(aes(yintercept = 1),colour = "black",linetype = "dashed")+
  geom_hline(aes(yintercept = -1),colour = "black",linetype = "dashed")+
  geom_vline(aes(xintercept = 11),colour = "black",linetype = "dotted")+
  geom_vline(aes(xintercept = 36),colour = "black",linetype = "dotted")+
  geom_vline(aes(xintercept = 38),colour = "black",linetype = "dotted")+
  labs(x = "Genes from four pathways", y = "  ",title = "Second factor")+
  theme(legend.justification = c(1,0),legend.position = c(1,0))

P3 = ggplot(load, aes(x=sample, y=factor3, colour=pathways, shape = pathways)) + 
  geom_point()+ 
  geom_text(aes(label = sample, vjust = 1.1, hjust = -0.5, angle = 0), show_guide = FALSE,check_overlap = TRUE)+
  theme(axis.text.x = element_blank())+
  theme(panel.grid =element_blank())+
  coord_cartesian(ylim=c(-5,5))+
  geom_hline(aes(yintercept = 0),colour = "grey",linetype = "dashed")+
  geom_hline(aes(yintercept = 1),colour = "black",linetype = "dashed")+
  geom_hline(aes(yintercept = -1),colour = "black",linetype = "dashed")+
  geom_vline(aes(xintercept = 11),colour = "black",linetype = "dotted")+
  geom_vline(aes(xintercept = 36),colour = "black",linetype = "dotted")+
  geom_vline(aes(xintercept = 38),colour = "black",linetype = "dotted")+
  labs(x = "Genes from four pathways", y = "  ",title = "Third factor")+
  theme(legend.justification = c(1,0),legend.position = c(1,0))


###
library("grid")
library("gridExtra")
library("lemon")
grid_arrange_shared_legend(P1,P2,P3,ncol = 3,position='bottom')

# refit: B

Plast = which(Y0$Pathwayname %in% c("Plastoquinonebiosynthesis"))
B.plast  = B.est2[,Plast]
order(abs(B.plast[,1]),decreasing=TRUE)[1:10]
order(abs(B.plast[,2]),decreasing=TRUE)[1:10]
which ((abs(B.plast[,1])>1))
which ((abs(B.plast[,2])>1))
abs(B.plast[,1])[order(abs(B.plast[,1]),decreasing=TRUE)[1:10]]


Carote = which(Y0$Pathwayname %in% c("Carotenoidbiosynthesis"))
B.Carote  = B.est2[,Carote]
for (i in 1:11) {
  print(i)
  print(which ((abs(B.Carote[,i])>1)))
  i=i+1
}


Chloro = which(Y0$Pathwayname %in% c("PorphyrinChlorophyllmetabolism"))
B.Chloro  = B.est2[,Chloro]
for (i in 1:length(Chloro)) {
  print(i)
  print(which ((abs(B.Chloro[,i])>1)))
  i=i+1
}


Phyto = which(Y0$Pathwayname %in% c("Phytosterolbiosynthesis"))
B.Phyto  = B.est2[,Phyto]
for (i in 1:length(Phyto)) {
  print(i)
  print(which ((abs(B.Phyto[,i])>1)))
  i=i+1
}


 