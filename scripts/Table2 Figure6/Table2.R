######### To produce Table 2
path_result <- "Results-Sim23_7_5_Arabidopsis/"
### read tables  PE 
file_name_PE_rop <- paste("PE_rop.txt", sep = "_")
file_path_PE_rop <- paste(path_result, file_name_PE_rop, sep = "")
PE_rop = data.table::fread(file_path_PE_rop)


file_name_PE_ropr <- paste("PE_ropr.txt", sep = "_")
file_path_PE_ropr <- paste(path_result, file_name_PE_ropr, sep = "")
PE_ropr = data.table::fread(file_path_PE_ropr)


file_name_PE_r4 <- paste("PE_r4.txt", sep = "_")
file_path_PE_r4 <- paste(path_result, file_name_PE_r4, sep = "")
PE_r4 = data.table::fread(file_path_PE_r4)


file_name_PE_rrr <- paste("PE_rrr3.txt", sep = "_")
file_path_PE_rrr <- paste(path_result, file_name_PE_rrr, sep = "")
PE_rrr = data.table::fread(file_path_PE_rrr)

file_name_PE_sss <- paste("PE_sss.txt", sep = "_")
file_path_PE_sss <- paste(path_result, file_name_PE_sss, sep = "")
PE_sss = data.table::fread(file_path_PE_sss)

file_name_PE_rcgl <- paste("PE_rcgl.txt", sep = "_")
file_path_PE_rcgl <- paste(path_result, file_name_PE_rcgl, sep = "")
PE_rcgl = data.table::fread(file_path_PE_rcgl)

file_name_PE_secure <- paste("PE_secure.txt", sep = "_")
file_path_PE_secure <- paste(path_result, file_name_PE_secure, sep = "")
PE_secure = data.table::fread(file_path_PE_secure) 


round(apply(PE_rop, 2, mean),3) 
round(apply(PE_ropr, 2, mean),3) 
round(apply(PE_r4, 2, mean)  ,3)
round(apply(PE_rrr, 2, mean),3)
round(apply(PE_sss, 2, mean),3)
round(apply(PE_rcgl, 2, mean),3)
round(apply(PE_secure, 2, mean),3)

