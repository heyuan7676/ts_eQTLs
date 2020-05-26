from GLOBAL_VAR import *
from METHOD_VAR import *
from sklearn.linear_model import LinearRegression



def readin_YW_BP(prefix, Nrows_to_read = None):
    Y = pd.read_csv('%s/%s_slope.txt' % (inputdir, prefix),  sep='\t', nrows = Nrows_to_read, index_col=[0,1])
    Pairs = np.array(Y.index)
    genes = [p[0] for p in Pairs]
    SNPs  = [p[1] for p in Pairs]
    Y.index = pd.MultiIndex.from_arrays([genes, SNPs])

    W = pd.read_csv('%s/%s_se.txt' % (inputdir, prefix),  sep='\t', index_col = [0,1], nrows = Nrows_to_read)
    W = 1/W
    Pairs = np.array(Y.index)
    genes = [p[0] for p in Pairs]
    SNPs  = [p[1] for p in Pairs]
    W.index = pd.MultiIndex.from_arrays([genes, SNPs])

    B = pd.read_csv('%s/%s.txt' % (ll_dir, LMfn), sep='\t', nrows = Nrows_to_read, index_col = [0,1])
    Y = Y.loc[B.index]    
    W = W.loc[B.index]

    print 'Number of pairs with >= 1 non-zero loading: %.5f' % (np.sum(np.sum(B!=0, axis=1) > 0) / float(len(B)))
    return [Y, W, B]





def compute_R2_one_pair(y, w, X, b):
        if np.sum(b!=0) == 0:
            return 0
        [T, K] = X.shape
        yw = y * w
        idx = np.where(~np.isnan(yw))[0]
        yw = yw[idx]
        w = w[idx]
        Xw = np.dot(X[:,idx], np.diag(w))
        Xw = Xw[np.where(b!=0)[0], :]
        lm = LinearRegression(fit_intercept = False)
        lm.fit(Xw.transpose(), yw)
        yhat = lm.predict(Xw.transpose())
	r2 = np.sum(np.power(yhat,2)) / np.sum(np.power(yw,2))
        return r2



def compute_R2(Yset, Wset, X, Bset):
    [T, K] = X.shape
    Yset = np.array(Yset)
    Wset = np.array(Wset)
    Bset = np.array(Bset)
    R2, nonzeroL, nonzeroT, minYhat, nanY = [], [], [], [],[]

    for k in range(len(Bset)):
        y = Yset[k]
        w = Wset[k]
        b = Bset[k]
        nonzeroL.append(np.sum(b!=0))
        nonzeroT.append(T - np.sum(y[~np.isnan(y)]!=0))
        R2.append(compute_R2_one_pair(y,w,X,b))
        
    return [np.array(R2), np.array(nonzeroL), np.array(nonzeroT)]




if __name__ == '__main__':

    [Y, W, B_correct] = readin_YW_BP(prefix, Nrows_to_read=None)

    np.random.seed(0)
    random_idx = np.random.choice(range(len(Y)), 10000, replace =False)
    pairs = Y.index[random_idx]

    [R2, nonzeroL, nonzeroT] = compute_R2(Y.loc[pairs], W.loc[pairs], np.array(X), B_correct.loc[pairs])
    a = np.array(R2)
    [thr0, thr02, thr06] = [np.sum(a>0) / float(len(a)), np.sum(a>0.2) / float(len(a)), np.sum(a>0.6) / float(len(a))]
    print ['thr0', 'thr02', 'thr06']
    print [thr0, thr02, thr06]

