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

    B = pd.read_csv('%s/%s.txt' % (ll_dir, LMfn.replace('_corrected','')), sep='\t', nrows = Nrows_to_read)
    Pairs = np.array(B.index)
    genes = [p.split(' ')[0] for p in Pairs]
    SNPs  = [p.split(' ')[1] for p in Pairs]
    B.index = pd.MultiIndex.from_arrays([genes, SNPs])

    Y = Y.loc[B.index]    
    W = W.loc[B.index]

    print 'Number of pairs with >= 1 non-zero loading: %.5f' % (np.sum(np.sum(B!=0, axis=1) > 0) / float(len(B)))
    return [Y, W, B]




def catch_compensate_loadings_one_pair(y, w, b):
    
    idx = list(np.where(b!=0)[0])
    compensate_loadings = []
    if 0 in idx:
        idx.remove(0)
    N  = len(idx)

    for i in range(N):
        k = idx[i]
        b_each = np.zeros(len(b))
        b_each[k] = b[k]

	y_each = np.dot(b_each, X) * w
	y_each[np.isnan(y_each)] = 0
        
	compare_tissues = np.where(np.abs(y_each) > 0.01)[0]
	y_each = y_each[compare_tissues]
	yw       = (y*w)[compare_tissues]
        
	compensate = np.zeros(len(compare_tissues))
	for j in range(len(compare_tissues)):
		if (np.sign(y_each[j]) * np.sign(yw[j]) == -1):
			compensate[j] = 1
		elif (np.sign(y_each[j]) * np.sign(yw[j]) == 1) and (np.abs(yw[j]) < 3):
			compensate[j] = 1

        if np.sum(compensate==1) > len(compensate) * 0.8:
            compensate_loadings.append(k)

    if len(compensate_loadings) > 0:
    	b[np.array(compensate_loadings)]  = 0

    return b
                    
        
        


def correct_B():

    B_correct = B.copy()
 
    for i in range(len(Y)):
        y = np.array(Y.iloc[i])
        w = np.array(W.iloc[i])
        b = np.array(B.iloc[i])

	#genei = 'ENSG00000000457.13'
	#snpi  = 'chr1_169877742'
	#y = np.array(Y.loc[genei].loc[snpi])
	#w = np.array(W.loc[genei].loc[snpi])
	#b = np.array(B.loc[genei].loc[snpi])

        nonEmpty = np.where(b!=0)[0]
        if len(nonEmpty) == 0:
            continue

        ## re-fit using the non-zero coefficients
        nonnan     = np.where(~np.isnan(y*w))[0]
	Xw = np.dot(X[:,nonnan], np.diag(w[nonnan]))
        lm = LinearRegression(fit_intercept = False)
        lm.fit(Xw[nonEmpty, :].transpose(), (y*w)[nonnan])
        fitted_b = np.zeros(len(b))
        fitted_b[nonEmpty] = lm.coef_

        b_correct = catch_compensate_loadings_one_pair(y, w, fitted_b)
	B_correct.iloc[i] = b_correct.copy()

    x1 = np.sum(np.abs(B)!=0,axis=0) / len(B)
    x2 = np.sum(np.abs(B_correct)!=0,axis=0) / len(B_correct)
    print 'Proportion of non-zero = %.5f / %.5f' % (np.mean(x1), np.mean(x2))
    
    return B_correct



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

    [Y, W, B] = readin_YW_BP(prefix, Nrows_to_read=None)

    B_correct = correct_B()
    B_correct.to_csv('%s/%s.txt' % (ll_dir, LMfn), sep='\t')


    np.random.seed(0)
    random_idx = np.random.choice(range(len(Y)), 1000, replace =False)
    pairs = Y.index[random_idx]

    [R2, nonzeroL, nonzeroT] = compute_R2(Y.loc[pairs], W.loc[pairs], X, B_correct.loc[pairs])
    a = np.array(R2)
    [thr0, thr02, thr06] = [np.sum(a>0) / float(len(a)), np.sum(a>0.2) / float(len(a)), np.sum(a>0.6) / float(len(a))]
    print ['thr0', 'thr02', 'thr06']
    print [thr0, thr02, thr06]

