save_results <- function(Betas, pValues, FM_fn, startidx = "NONE"){

    ### save result
    outdir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/Revision_Zscore/LL/'

    FM_fn = gsub('.mat', '', FM_fn)

    if(startidx != "NONE"){
        FM_fn = paste0(FM_fn, '_startidx', startidx)
    }


    outFn = paste0(outdir, FM_fn, '_Loadings_pvalue.txt')
    print(outFn)
    write.table(pValues, outFn, sep='\t')

    outFn = paste0(outdir, FM_fn, '_Loadings_beta.txt')
    print(outFn)
    write.table(Betas, outFn, sep='\t')

}

