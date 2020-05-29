from GLOBAL_VAR import *

group = -1
pair_Fn = '%s/%s_outlierPairs_group%d.txt' % (pairdir, LMfn, group)
N0 = len(pd.read_csv(pair_Fn, header= None, sep='\t', usecols=[0]))


factor_names_thr = ['Adipose', 'Adrenal_Gland', 'Artery',
                'Brain', 'LCL', 'Cells_Cultured_fibroblasts', 
                'Colon', 'Esophagus', 'Heart',
                'Kidney_Cortex', 'Liver', 'Lung',
                'Minor_Salivary_Gland', 'Muscle_Skeletal', 'Nerve_Tibial',
                'Ovary', 'Pancreas', 'Pituitary', 'Prostate',
                'Skin', 'Small_Intestine_Terminal_Ileum', 'Spleen',
                'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina', 'Whole_Blood']

pairs_N = []
for group in range(28):
    pair_Fn = '%s/%s_outlierPairs_group%d.txt' % (pairdir, LMfn, group)
    N = len(pd.read_csv(pair_Fn, header= None, sep='\t', usecols=[0]))
    for tis in [x for x in tissues if factor_names_thr[group] in x]:
        pairs_N.append([N, tis, N0])

pairs_N = pd.DataFrame(pairs_N)
pairs_N.columns = ['ts_N', 'tissue', 'shared_N']
pairs_N.to_csv('Fig2_sig_prop_%s.txt' % FMfn, sep='\t', index = False)
