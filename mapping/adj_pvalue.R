adj_pvalue <- function(Betas, pValues, FM_fn, outdir, FDR_alpha = 0.05){

	### save result
	dir.create(outdir, showWarnings = FALSE)

    	outFn = paste0(outdir, 'mapping_', FM_fn, '_Loadings_pvalue.txt')
    	print(outFn)
    	write.table(pValues, outFn, sep='\t')

    	outFn = paste0(outdir, 'mapping_', FM_fn, '_Loadings_beta.txt')
    	print(outFn)
    	write.table(Betas, outFn, sep='\t')


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
        outFn = paste0(outdir, 'mapping_', FM_fn, '_Loadings_pvalue_BH.txt')
        print(outFn)
        write.table(adj_p, outFn, sep='\t')

        outFn = paste0(outdir, 'mapping_', FM_fn, '_Loadings_beta_alpha',FDR_alpha,'.txt')
        print(outFn)
        write.table(adj_beta, outFn, sep='\t')

}

