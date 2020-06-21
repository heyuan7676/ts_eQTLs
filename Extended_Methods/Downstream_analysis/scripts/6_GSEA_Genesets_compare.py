from GLOBAL_VAR import *
from METHOD_VAR import *

### Test number of ts-eGenes
Tested_Method_N = {}
for geneGroup in ['topGene4','topGene6', 'topGene5', 'topGene30']:
    Tested_Method_N[geneGroup] = []
    for tis in range(1,23):
        tp = pd.read_csv('%s/%s_%s_group%d.txt' % (pairdir, LMfn, geneGroup, tis), sep='\t')
        Tested_Method_N[geneGroup].append(len(tp))


# spMF

old_pairdir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/downstream/pairSets/'

r2 = '1'
FDR = 0.05

FMfn = 'SparseMF_coph_%s_topPair_K30_a11_l110' % prefix.replace(r2, '0.2')
LMfn= '%s%s_Loadings_beta_BH_corrected_alpha%s' % (FMfn, LDprefix, str(FDR))

spMF_N = {}

for geneGroup in ['topGene5', 'topGene22']:
    spMF_N[geneGroup] = []
    for k in range(1, 23):
        tp = pd.read_csv('%s/%s_%s_group%d.txt' % (old_pairdir, LMfn, geneGroup, k), sep='\t')
        spMF_N[geneGroup].append(len(tp))
      
 
x = Tested_Method_N['topGene4']
y = spMF_N['topGene5']
print(np.mean(x), np.mean(y))


x = Tested_Method_N['topGene30']
y = spMF_N['topGene22']
print(np.mean(x), np.mean(y))

