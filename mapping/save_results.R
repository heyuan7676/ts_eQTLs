save_results <- function(Betas, pValues, FM_fn, prefix = "", startidx = "NONE"){

    ### save result
    outdir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/LL/'

    FM_fn = gsub('.mat', '', FM_fn)
    if(startidx != "NONE"){
        FM_fn = paste0(FM_fn, '_startidx', startidx)
    }

    outFn = paste0(outdir, FM_fn, prefix, '_Loadings_pvalue.txt')
    print(outFn)
    write.table(pValues, outFn, sep='\t')

    outFn = paste0(outdir, FM_fn, prefix, '_Loadings_beta.txt')
    print(outFn)
    write.table(Betas, outFn, sep='\t')

}

