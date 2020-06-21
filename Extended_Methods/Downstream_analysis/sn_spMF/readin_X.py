import os
import sys
#sys.path.append('/home-4/yhe23@jhu.edu/work/yuan/tools/python_lib/lib/python2.7/site-packages')

import pandas as pd
import numpy as np
from scipy.io import loadmat

mfDir = '/work-zfs/abattle4/heyuan/tissue_spec_eQTL_v8/FL/coph'

def prototype_f(X, tissues):
    P_V, F_tissues  = [], []
    X = X / np.max(X, 1)[np.newaxis].transpose()
    [K, T] = X.shape
    for t in range(1, T):
        proto_vector = np.zeros(T)
        proto_vector[:t] = 1.0
        P_V.append(proto_vector)
    for x in X:
        dis = [np.linalg.norm(v - np.sort(x)[::-1]) for v in P_V]
        tissue_in_f =  np.argsort(x)[::-1][:np.argmin(dis)+1]
        F_tissues.append(tissues[tissue_in_f])
    return F_tissues




def tissues_factor(X, tissues):
    Comp_tissues = []
    for x in X:
        Comp_tissues.append(tissues[np.where(x>0.01)])
    return Comp_tissues



def readin_X(FLfn, save_clean_X = False):
    '''
    X dim: D x T 
    '''
    data= loadmat('%s/%s.mat' % (mfDir, FLfn))
    X = data['Xall'][0][0]
    tissues = pd.read_csv('tissues.txt', sep='\t', header=None)
    tissues = np.array(tissues[0])
    
    scale_X = np.median(np.max(X, axis=0))
    X = np.array([scale_X/np.max(x) * x for x in X.transpose()])
    assert np.sum((np.unique(np.max(X, axis=1)) - np.max(X)) > 1e-10) == 0
    
    Comp_tissues = tissues_factor(X, tissues)
    
    if save_clean_X:
        clean_X = np.zeros(X.shape)
        for k in range(len(Comp_tissues)):
            c = Comp_tissues[k]
            remain_tissue_idx = np.array([list(tissues).index(t) for t in c])
            clean_X[k, remain_tissue_idx]  = X[k, remain_tissue_idx]
            
            clean_X[0, :] = np.ones(T)
            data['Xall'][0,0] = clean_X_2.transpose()
            savemat('%s/%s_cleanX_2.mat' % (mfDir, FLfn), data)

    return [X, Comp_tissues, tissues]
