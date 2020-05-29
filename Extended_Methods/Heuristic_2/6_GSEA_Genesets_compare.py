from GLOBAL_VAR import *

### Test number of ts-eGenes

# 

Tested_Method_N = {}
for geneGroup in ['topGene4', 'topGene5', 'topGene6']:
    Tested_Method_N[geneGroup] = []
    for tis in range(28):
        tp = pd.read_csv('%s/%s_%s_group%d.txt' % (pairdir, LMfn, geneGroup, tis), sep='\t')
        Tested_Method_N[geneGroup].append(len(tp))


# spMF

r2 = '1'
FDR = 0.05

fmDir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/FL/coph'
ll_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/LL'
prefix       = 'v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete_filteredSNPs.LDblocks_%s' % str(r2)
LDprefix     = '_LD1'

FMfn = 'SparseMF_coph_%s_topPair_K30_a11_l110' % prefix.replace(r2, '0.2')
LMfn= '%s%s_Loadings_beta_BH_corrected_alpha%s' % (FMfn, LDprefix, str(FDR))

spMF_N = {}

for geneGroup in ['topGene5', 'topGene6', 'topGene7']:
    spMF_N[geneGroup] = []
    for k in range(1, 23):
        tp = pd.read_csv('%s/%s_%s_group%d.txt' % (pairdir, LMfn, geneGroup, k), sep='\t')
        spMF_N[geneGroup].append(len(tp))
       


# flashr
LMfn = 'flashr_Loadings_beta_BH_alpha0.05_corrected'

flashr_N = {}

for geneGroup in ['topGene5']:
    flashr_N[geneGroup] = []
    for k in range(1, 23):
        tp = pd.read_csv('%s/%s_%s_group%d.txt' % (pairdir, LMfn, geneGroup, k), sep='\t')
        flashr_N[geneGroup].append(len(tp)) 
        
x = Tested_Method_N['topGene5']
y = spMF_N['topGene5']
z = flashr_N['topGene5']
print(np.mean(x), np.mean(y), np.mean(z))

