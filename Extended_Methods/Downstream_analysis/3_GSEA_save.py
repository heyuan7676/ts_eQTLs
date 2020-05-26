from GLOBAL_VAR import *
from METHOD_VAR import *
from statsmodels.stats import multitest

fn=open(gsea_fn, 'r')

gsea_pathways = {}
gsea_pathways_number = {}

for line in fn.readlines():
    gsea_pathways[line.split('\t')[0]] = line.split('\t')[2:]
    x  = np.ceil(int(len(line.split('\t')[2:]))/10)
    if  x in gsea_pathways_number.keys():
        gsea_pathways_number[x].append(line.split('\t')[0])
    else:
        gsea_pathways_number[x] = [line.split('\t')[0]]

fn.close()




### weighted GO terms
go_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/annotations/TFmotif/weighted_GO_terms/'
fn = 'SetTestWeights_C5.BP_RNA.and.protein.txt'
dat = pd.read_csv('%s/%s' % (go_dir, fn), sep='\t')

dat = dat.apply(lambda x: (x - np.min(x)) / (np.max(x) - np.min(x)), axis=0)

dat.columns = ['Adipose_Subcutaneous', 'Adrenal_Gland', 'appendix', 'bone marrow', 
             'Breast_Mammary_Tissue', 'Brain_Cortex', 'Uterus', 'Colon_Sigmoid',
             'duodenum', 'endometrium', 'epididymis', 'Esophagus_Mucosa', 'fallopian tube',
      'gallbladder', 'Heart_Left_Ventricle', 'Kidney_Cortex', 'Liver', 'Lung', 
            'lymph node', 'ovary', 'Pancreas', 'parathyroid gland', 'placenta', 'prostate',
      'rectum', 'salivary gland', 'seminal vesicle', 'Muscle_Skeletal', 'Skin_Not_Sun_Exposed_Suprapubic',
             'Small_Intestine_Terminal_Ileum', 'smooth muscle', 'Spleen', 'Stomcah', 
             'Testis', 'Thyroid', 'tonsil', 'urinary bladder']

dat['Esophagus_Muscularis'] = dat['Esophagus_Mucosa']
dat['Adipose_Visceral_Omentum'] = dat['Adipose_Subcutaneous']
dat['Brain_Cerebellar_Hemisphere'] = dat['Brain_Cortex']
dat['Brain_Caudate_basal_ganglia'] = dat['Brain_Cortex']

dat['Pituitary'] = dat['Brain_Cortex']
dat['Whole_Blood'] = dat['bone marrow']


dat = dat[np.intersect1d(dat.columns, tissues)]



def readin_gsea_result(geneGroup, HIT, lower = 20, higher = 500):
    gsea_dat = pd.DataFrame()

    for g in range(len(Comp_tissues)):
        gsea_fn = '%s/%s_%s_group%d_enrichment.txt' % (gsea_dir, LMfn,geneGroup, g)
        try:
            gsea_dat_tis = pd.read_csv(gsea_fn, sep='\t', index_col = 0)
            if np.sum(gsea_dat_tis['test_inset']!=0) == 0:
                continue
        except:
            continue

        gsea_dat_tis['group'] = g
        gsea_dat = gsea_dat.append(gsea_dat_tis)
    
    gsea_dat = gsea_dat[gsea_dat['test_inset'] + gsea_dat['background_inset'] > HIT]
    gsea_dat['N_genes'] = [len(gsea_pathways[g.replace(' ','')]) for g in np.array(gsea_dat.index)]
    gsea_dat = gsea_dat[gsea_dat['N_genes'] < higher]
    gsea_dat = gsea_dat[gsea_dat['N_genes'] > lower]
    
    gsea_dat['BH_p'] = multitest.multipletests(gsea_dat['pvalue'], method = 'fdr_bh')[1]
    gsea_dat['Bonf_p'] = multitest.multipletests(gsea_dat['pvalue'])[1]
    gsea_dat = gsea_dat.sort_values('pvalue')
    

    return gsea_dat




### save
HIT = 10
LOWER = 20
UPPER = 500
alpha = 0.05

topGeneN = int(sys.argv[1])

mf_genes = readin_gsea_result('topGene%d' % topGeneN,HIT, lower=LOWER, higher=UPPER)
mf_genes['BH_p'] = multitest.multipletests(mf_genes['pvalue'], method = 'fdr_bh')[1]
suffix = '%s_stringent' % FMfn
mf_genes.to_csv('%s/Fig3_GSEA_%s.txt' % (plot_save_dir, suffix), sep='\t', index=True)
x = mf_genes[mf_genes['group'] != 0]
print('Stringent definition: Total enriched = %d, in ts-eGenes = %d' % (len(set(mf_genes[mf_genes['BH_p'] < alpha].index)), len(set(x[x['BH_p'] < alpha].index))))



mf_genes = readin_gsea_result('topGene30',HIT, lower=LOWER, higher=UPPER)
mf_genes['BH_p'] = multitest.multipletests(mf_genes['pvalue'], method = 'fdr_bh')[1]
suffix = '%s' % FMfn
mf_genes.to_csv('%s/Fig3_GSEA_%s.txt' % (plot_save_dir, suffix), sep='\t', index=True)
x = mf_genes[mf_genes['group'] != 0]
print('Total enriched = %d, in ts-eGenes = %d' % (len(set(mf_genes[mf_genes['BH_p'] < alpha].index)), len(set(x[x['BH_p'] < alpha].index))))

