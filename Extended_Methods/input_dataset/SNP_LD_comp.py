import sys
sys.path.append('/home-4/yhe23@jhu.edu/work/yuan/tools/python_lib/lib/python2.7/site-packages')
sys.path.append('/home-4/yhe23@jhu.edu/work/yuan/tools/python_lib/lib/python2.7/site-packages/lib/python2.7/site-packages')

import pdb
import networkx as nx
import pickle

r2_thr = str(sys.argv[1])

LD_block_dir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/datasets/GTEx_datasets/genotypes'
fn = '%s/plink_r2_%s.ld' % (LD_block_dir, r2_thr)

LD_fn = open(fn, 'rb')
graph=nx.read_edgelist(LD_fn)
LD_fn.close()
print 'Edges added'

outfn = '%s.graphs.pickle' % fn 
with open(outfn, 'wb') as output:
	pickle.dump(graph, output, protocol=pickle.HIGHEST_PROTOCOL)

print 'Saved the graph'
