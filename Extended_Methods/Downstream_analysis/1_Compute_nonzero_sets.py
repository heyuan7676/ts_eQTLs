from GLOBAL_VAR import *
from METHOD_VAR import *
from collections import Counter

def readin_geneannotation():
    ### read in gene annotation
    gene_annotation = {}
    gene_annotation_fn = open('/work-zfs/abattle4/lab_data/annotation/gencode.v26/gencode.v26.annotation.gene.txt','r')
    for l in gene_annotation_fn.readlines():
        ensemID = l.rstrip().split('\t')[0].replace('"','')
        geneName = l.rstrip().split('\t')[5]
        gene_annotation[ensemID] = geneName
    return gene_annotation



def readin_LDblocks():
    ld_pairs = dict()
    for ith in range(39):
        data = np.load('%s/npz_ith/%s_ith_%d.npz' % (inputdir, prefix, ith))
        ld_pairs_ith = data['arr_0'].item()
        ld_pairs.update(ld_pairs_ith)
    return ld_pairs



def pairSet():

    B = pd.read_csv('%s/%s.txt' % (ll_dir, LMfn), sep='\t', index_col = [0, 1], nrows = None)
    [N, K] = B.shape
    B.columns = range(K)

    ### read in the Loading matrix
    pairs = np.array([' '.join((p[0], p[1])) for p in np.array(B.index)])
    snpsets = {}    

    ### map back to all SNPs
    ld_pairs        = readin_LDblocks()
    for k in range(K):
        good_sig_pair_idx = np.where(B[k]!=0)[0]
        outlierpairs = pairs[good_sig_pair_idx]

        keyPairs = np.intersect1d(list(ld_pairs.keys()), outlierpairs)
        alloutlierPairs = [ld_pairs[p] for p in keyPairs]
        alloutlierPairs = np.unique([a for b in alloutlierPairs for a in b] + list(outlierpairs))
        outliergenes = [x.split(' ')[0] for x in alloutlierPairs]
	snpsets[k]   = [x.split(' ')[1] for x in alloutlierPairs]

        print('Number of outlier in loading %d: %d -> %d' % (k, len(outlierpairs), len(alloutlierPairs)))

        f = open('%s/%s_outlierPairs_group%d.txt' % (pairdir, LMfn, k), 'w')
        for p in alloutlierPairs:
            f.write('%s\t%s\n' % (p.split(' ')[0], p.split(' ')[1]))
        f.close()



if __name__ == '__main__':
    K = len(Comp_tissues)
    gene_annotation = readin_geneannotation()
    pairSet()

