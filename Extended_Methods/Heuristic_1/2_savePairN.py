from GLOBAL_VAR import *

group = 0
pair_Fn = '%s/%s_outlierPairs_group%d.txt' % (pairdir, LMfn, group)
N0 = len(pd.read_csv(pair_Fn, header= None, sep='\t', usecols=[0]))


pairs_N = []
for group in range(1, 23):
    pair_Fn = '%s/%s_outlierPairs_group%d.txt' % (pairdir, LMfn, group)
    N = len(pd.read_csv(pair_Fn, header= None, sep='\t', usecols=[0]))
    pairs_N.append([N, N0])

pairs_N = pd.DataFrame(pairs_N)
pairs_N.columns = ['ts_N', 'shared_N']
pairs_N.to_csv('Fig2_sig_prop_%s.txt' % FMfn, sep='\t', index = False)
