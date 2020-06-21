from GLOBAL_VAR import *
from METHOD_VAR import *

from scipy.stats import fisher_exact 
from scipy.stats import chi2_contingency
from statsmodels.stats import multitest



def one_factorActive(outlier_snps, in_tissues, sub_feature = None):
        outlier_snps_ts = np.array(outlier_snps[in_tissues])
        if ',' in str(outlier_snps_ts[0]):
                snp_total = np.array([xx != '.,' for xx in outlier_snps_ts]) * 1
        else:
                snp_total = np.array([xx != '.' for xx in outlier_snps_ts]) * 1

        Active_fc = np.array([sub_feature in xx for xx in outlier_snps_ts]) * 1
        return [snp_total, Active_fc]



def factorsActive(featureName, sub_feature = None):
	featureSNpset, featureSNpset_rd, SNPset_size, SNPset_size_rd = dict(), dict(), dict(), dict()
	mT, unmT = dict(), dict()

        for group in range(len(Comp_tissues)):
	#for group in range(3):
		print(group)
                try:
                        outlier_fc = pd.read_csv('%s/%s_outlierSNPs_%s_group%d.txt' % (sigSNPfeaturedir, FMfn, featureName, group), sep='\t', index_col = 0)
			outlier_fc.columns = [x.replace('chr1_','') for x in outlier_fc.columns]
                except:
                        continue
		# deal with duplicated column names
		if group != 0:
			existing_tissues = [x for x in outlier_fc.columns if x.split('.')[0] in '_'.join(Comp_tissues[group])]
		else:
			existing_tissues = outlier_fc.columns
                if len(existing_tissues) == 0:
                        continue
                ### background: random SNPs matched for MAF
		outlier_rd = pd.read_csv('%s/%s_outlierSNPs_random_matched_%s_group%d_4folds.txt' % (sigSNPfeaturedir, FMfn, featureName, group), sep='\t', index_col = 0)
		outlier_rd.columns = [x.replace('chr1_','') for x in outlier_rd.columns]

		## intersect with open regions
		if intersectWithOpen:
			try:
				outlier_fc_open = pd.read_csv('%s/%s_outlierSNPs_%s_group%d.txt' % (sigSNPfeaturedir, FMfn, 'DNase', group), sep='\t', index_col = 0)
                        	outlier_fc_open.columns = [x.replace('chr1_','') for x in outlier_fc_open.columns]
			except:
				continue
			outlier_rd_open = pd.read_csv('%s/%s_outlierSNPs_random_matched_%s_group%d_4folds.txt' % (sigSNPfeaturedir, FMfn, 'DNase', group), sep='\t', index_col = 0)
			outlier_rd_open.columns = [x.replace('chr1_','') for x in outlier_rd_open.columns]
              
			assert np.sum(outlier_rd.index != outlier_rd_open.index) == 0
			assert np.sum(outlier_fc.index != outlier_fc_open.index) == 0

			tissue_with_open = np.intersect1d( outlier_fc_open.columns, outlier_fc.columns)
			outlier_fc_open = outlier_fc_open[tissue_with_open]
			outlier_rd_open = outlier_rd_open[tissue_with_open]
			outlier_fc      = outlier_fc[tissue_with_open]
			outlier_rd      = outlier_rd[tissue_with_open]

			if group != 0:
				existing_tissues = np.intersect1d(Comp_tissues[group], outlier_fc.columns)
			else:
				existing_tissues = outlier_fc.columns
	                if len(existing_tissues) == 0:
				continue

			fc_array = np.array(outlier_fc)
			fc_array[np.where(np.array(outlier_fc_open) == '0,')] = '0'
			fc_df = pd.DataFrame(fc_array)
			fc_df.columns = outlier_fc.columns
			fc_df.index = outlier_fc.index

                        rd_array = np.array(outlier_rd)
                        rd_array[np.where(np.array(outlier_rd_open) ==  '0,')] = '0'
                        rd_df = pd.DataFrame(rd_array)
                        rd_df.columns = outlier_rd.columns
                        rd_df.index = outlier_rd.index

			outlier_fc = fc_df
			outlier_rd = rd_df

		featureSNpset[group]            = [[], []]
		featureSNpset_rd[group]         = [[], []]
		SNPset_size[group]		= [[], []]
		SNPset_size_rd[group]           = [[], []]

		mT[group] = []
		unmT[group] = []

		## matched tissues
		for matched_tissue in existing_tissues:
                	[N, Active_fc] = one_factorActive(outlier_fc, matched_tissue, sub_feature)
                	featureSNpset[group][0].append(np.sum(Active_fc))
			SNPset_size[group][0].append(np.sum(N))

                	[N, Active_rd] = one_factorActive(outlier_rd, matched_tissue, sub_feature)
                	featureSNpset_rd[group][0].append(np.sum(Active_rd))
			SNPset_size_rd[group][0].append(np.sum(N))

			mT[group].append(matched_tissue)

		## unmatched tissues
		non_in_tissues = np.intersect1d(list(set(tissues) - set(existing_tissues)), outlier_fc.columns)
		if len(non_in_tissues) == 0:
			continue
		for unmatched_tissue in non_in_tissues:
                        [N, Active_fc] = one_factorActive(outlier_fc, unmatched_tissue, sub_feature)
                        featureSNpset[group][1].append(np.sum(Active_fc))
			SNPset_size[group][1].append(np.sum(N))

                        [N, Active_rd] = one_factorActive(outlier_rd, unmatched_tissue, sub_feature)
                        featureSNpset_rd[group][1].append(np.sum(Active_rd))
			SNPset_size_rd[group][1].append(np.sum(N))

			unmT[group].append(unmatched_tissue)
	return [featureSNpset, featureSNpset_rd, SNPset_size, SNPset_size_rd, mT, unmT]




def contig_test(a, b, c, d):
        #print [[a,b],[c,d]]
	a = float(a)
	b = float(b)
	c = float(c)
	d = float(d)
	if c == 0:
		a = a+1
		b = b+1
		c = c+1
		d = d+1
        return [(a/b)/(c/d), chi2_contingency([[a,b],[c,d]])[1]]



def active_enrichment_compared_to_random(active, rd_active, SNPset_size, rd_SNPset_size, mT, unmT):
        inoutN, oddsRatio, pvalue, groupName, matchedGroup, match_or_not = [], [], [], [], [], []
        for group in list(active.keys()):
		# matched tissues
		for t in range(len(active[group][0])):
			A  = SNPset_size[group][0][t]
			B  = rd_SNPset_size[group][0][t]
			pN = active[group][0][t]
			qN = rd_active[group][0][t]
			inoutN.append([group, 1, pN, A-pN, qN, B-qN])
			enrichTest1 = contig_test(pN, A-pN, qN, B-qN)
                	oddsRatio.append(enrichTest1[0])
                	pvalue.append(enrichTest1[1])
                	groupName.append(group)
			matchedGroup.append(mT[group][t])
			match_or_not.append(1)

                # unmatched tissues
                for t in range(len(active[group][1])):
                        A  = SNPset_size[group][1][t]
                        B  = rd_SNPset_size[group][1][t]
                        pN = active[group][1][t]
                        qN = rd_active[group][1][t]
			inoutN.append([group, -1, pN, A-pN, qN, B-qN])
			enrichTest1 = contig_test(pN, A-pN, qN, B-qN)
                        oddsRatio.append(enrichTest1[0])
                        pvalue.append(enrichTest1[1])
                        groupName.append(group)
                        matchedGroup.append(unmT[group][t])
			match_or_not.append(0)

        ### format
        enrichment = pd.DataFrame([oddsRatio, pvalue, groupName, matchedGroup, match_or_not]).transpose()
        enrichment.columns = ['Oddsratio', 'pvalue', 'group', 'inTissue', 'matched']
        enrichment = enrichment.sort_values('pvalue')
	tested_N   = pd.DataFrame(np.array(inoutN))
	tested_N.columns = ['group', 'matched', 'inTest', 'notinTest', 'bg_inTest', 'bg_notinTest']
        return [enrichment, tested_N] 




def wrap_up(featureName, sub):
        [ts_pairs, rd_pairs, SNPset_size, rd_SNPset_size, mT, unmT] = factorsActive(featureName, sub_feature = sub)
	if sub == 'None':
		print SNPset_size
        print featureName, sub

        [enrichment_results, tested_N] = active_enrichment_compared_to_random(ts_pairs, rd_pairs, SNPset_size, rd_SNPset_size, mT, unmT)
	print(enrichment_results)
	background = 'random'
        if intersectWithOpen:
                enrichment_results.to_csv('%s/%s_Active_enrichment_%s_%s_to%s_WithInOpen.txt' % (outdir, FMfn, featureName, sub.replace('/','_'), background), index=False)
                tested_N.to_csv('%s/%s_Active_enrichment_N_%s_%s_to%s_WithInOpen.txt' % (outdir, FMfn, featureName, sub.replace('/','_'), background), index=False)
        else:
                enrichment_results.to_csv('%s/%s_Active_enrichment_%s_%s_to%s.txt' % (outdir, FMfn, featureName, sub.replace('/','_'), background), index=False)
                tested_N.to_csv('%s/%s_Active_enrichment_N_%s_%s_to%s.txt' % (outdir, FMfn, featureName, sub.replace('/','_'), background), index=False)


if __name__ == '__main__':
        outdir       = '%s/active_regions' % enrichment_dir
	intersectWithOpen = int(sys.argv[1])

	for sub in ['1_TssA', '7_Enh']:
		wrap_up('ROADMAP', sub)
   

