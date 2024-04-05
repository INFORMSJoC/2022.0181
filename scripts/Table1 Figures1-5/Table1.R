#### Table 1 
path_result <- "Results-Sim22_5_1_context/"
# draw table
set = "setting1"  #"setting2"
method_list <- c("ROP", "ROPr", "RRR","RCGL","SECURE","R4")
file_name_BC <- paste(set, "BC.txt", sep = "_")
file_path_BC <- paste(path_result, file_name_BC, sep = "")
rt1 = read.table(file_path_BC)
#index = which(rt1[,2]<rt1[,6])
index = 1:100
errorBC = as.numeric(round(apply(rt1[index,],2,mean), 3))
sdBC = as.numeric(round(apply(rt1[index,],2,sd), 3))

file_name_B <- paste(set, "B.txt", sep = "_")
file_path_B <- paste(path_result, file_name_B, sep = "")
rt = read.table(file_path_B)
errorB = as.numeric(round(apply(rt[index,],2,mean), 3))
sdB = as.numeric(round(apply(rt[index,],2,sd), 3))

file_name_C <- paste(set, "C.txt", sep = "_")
file_path_C <- paste(path_result, file_name_C, sep = "")
rt = read.table(file_path_C)
errorC = as.numeric(round(apply(rt[index,],2,mean), 3))
sdC = as.numeric(round(apply(rt[index,],2,sd), 3))

file_name_rank <- paste(set, "rank.txt", sep = "_")
file_path_rank <- paste(path_result, file_name_rank, sep = "")
rt = read.table(file_path_rank)
errorRank = as.numeric(round(apply(rt[index,],2,mean), 3))

file_name_time <- paste(set, "time.txt", sep = "_")
file_path_time <- paste(path_result, file_name_time, sep = "")
rt = read.table(file_path_time)
errorTime = as.numeric(round(apply(rt[index,],2,mean), 3))

file_name_mask <- paste(set, "mask.txt", sep = "_")
file_path_mask <- paste(path_result, file_name_mask, sep = "")
rt = read.table(file_path_mask)
errorMask = as.numeric(round(apply(rt[index,],2,mean), 3))


file_name_swamp <- paste(set, "swamp.txt", sep = "_")
file_path_swamp <- paste(path_result, file_name_swamp, sep = "")
rt = read.table(file_path_swamp)
errorSwamp = as.numeric(round(apply(rt[index,],2,mean), 3))

file_name_detection <- paste(set, "detection.txt", sep = "_")
file_path_detection <- paste(path_result, file_name_detection, sep = "")
rt = read.table(file_path_detection)
errorDetection = as.numeric(round(apply(rt[index,],2,sum)/length(index), 3))

RE = cbind(errorB,sdB,errorC,sdC, errorRank, errorMask, errorSwamp, errorDetection, errorTime)
rownames(RE)= method_list
print(RE)
# Note that RCGL and R4 are stored in minutes for errorTime, therefore the "Time" in Table 1 needs to mulitpy 60