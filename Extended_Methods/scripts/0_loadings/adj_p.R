library(data.table)

args <- commandArgs(TRUE)
FDR_alpha = as.numeric(args[1])
FMfn = args[2]

LLdir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/Revision_Zscore/LL/'

### read in 
pV_fn = paste0(LLdir, FMfn, '_Loadings_pvalue.txt')
bT_fn = paste0(LLdir, FMfn, '_Loadings_beta.txt')

pValues = fread(pV_fn, sep='\t')
pValues = data.frame(pValues, row.names=1)
pValues[is.na(pValues)] = -1

Betas = fread(bT_fn, sep='\t')
Betas = data.frame(Betas, row.names=1)

### BH-correction

# ignore those with missing pvalue

p = unlist(pValues)
idx = which(p!=-1)

adj_p = rep(-1, length(p))
adj_p[idx] = p.adjust(p[idx], method = 'BH')
adj_p = data.frame(matrix(adj_p, ncol = dim(pValues)[2]))
rownames(adj_p) = rownames(pValues)
stopifnot(sum(which(adj_p == -1) != which(pValues==-1)) == 0)

sig_pairs = (adj_p >= 0) * (adj_p < FDR_alpha)
adj_beta = Betas * sig_pairs
adj_beta[is.na(adj_beta)] = 0
rownames(adj_beta) = rownames(Betas)

### save results
pV_fn = paste0(LLdir, FMfn, '_Loadings_pvalue_BH.txt')
bT_fn = paste0(LLdir, FMfn, '_Loadings_beta_BH_alpha',FDR_alpha,'.txt')

write.table(adj_p, pV_fn, sep='\t')
write.table(adj_beta, bT_fn, sep='\t')



