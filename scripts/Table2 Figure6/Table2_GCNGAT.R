# For Table 2
# GCN
gnn_re = read.csv("gnn_results_gene_10_linear_unscaled.csv")
Re_mat = matrix(nrow = 6, ncol = 2)
Re_mat[,1] = c(4,6,8,16,32,64)
i = 1
for(dim in c(4,6,8,16,32,64)){
  Re_mat[i ,2] = mean(gnn_re[which(gnn_re[,2] == dim),3])
  i = i+1
}
round(min(Re_mat[,2]),3)

# GAT
GAT_re = read.csv("GAT_results_gene_10.csv")
Re_mat = matrix(nrow = 6, ncol = 2)
Re_mat[,1] = c(4,6,8,16,32,64)
i = 1
for(dim in c(4,6,8,16,32,64)){
  Re_mat[i ,2] = mean(GAT_re[which(GAT_re[,2] == dim),3])
  i = i+1
}
round(min(Re_mat[,2]),3)
