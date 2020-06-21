from GLOBAL_VAR import *

r2 = '1'
FDR = 0.05

fmDir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/FL/coph'
ll_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/LL'
prefix       = 'v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete_filteredSNPs.LDblocks_%s' % str(r2)
LDprefix     = '_LD1'

FMfn = 'SparseMF_coph_%s_topPair_K30_a11_l110' % prefix.replace(r2, '0.2')
LMfn= '%s%s_Loadings_beta_BH_corrected_alpha%s' % (FMfn, LDprefix, str(FDR))


shared_text = 0
[X, Comp_tissues, tissues] =  readin_X()
