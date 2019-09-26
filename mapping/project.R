source("readin_YWX.R")
source("run_projection.R")

args <- commandArgs(TRUE)
FM_fn    = args[1]
inputprefix   = args[2]
LDr2 = args[3]
startidx = args[4]

## read in data
data  = readin_data(FM_fn, paste0(inputprefix, '_', LDr2), Z_score = T, startidx)
Genes = data[[1]]
SNPs  = data[[2]]
eqtl_mat   = data[[3]]
factor_mat = data[[4]]
comp = data[[5]]


pve_mat <- project_eqtl_to_flash_factor(as.matrix(eqtl_mat), factor_mat, full = T)$full_variance_explained
membership <- get_membership_from_pve_mat(pve_mat)
membership <- membership[which(!is.na(membership))]


name_to_idx = min(as.numeric(names(membership))) - 1
assign_mx = matrix(0, ncol = dim(factor_mat)[2], nrow = dim(eqtl_mat)[1])
for(i in seq(1, length(membership))){
  assign_mx[as.numeric(names(membership))[i] - name_to_idx, as.vector(membership)[i]] = 1
}

assign_mx = as.data.frame(assign_mx)
rownames(assign_mx) = paste(Genes, SNPs)

## save results
outdir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/LL/'

FM_fn = gsub('.mat', '', FM_fn)
if(startidx != "NONE"){
  FM_fn = paste0(FM_fn, '_startidx', startidx)
}

outFn = paste0(outdir, FM_fn, paste0('_LD', LDr2), '_Loadings_projection.txt')
print(outFn)
write.table(assign_mx, outFn, sep='\t')


