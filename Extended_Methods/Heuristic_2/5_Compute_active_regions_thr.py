from GLOBAL_VAR import *


### Produce:
### for each loading: all active SNPs (including zero & nonzero) in matched tissues of this loading



def readInGenomicAnnotation_factors(active_regions, featureName, sub):
	print(featureName)
	active_regions.columns = [x.replace('chr1_','') for x in active_regions.columns]
	#for group in range(len(Comp_tissues)):
	for group in list(range(len(Comp_tissues))) + [-1]:
		print(group)
		if group != -1:
                	existing_tissues = np.intersect1d(Comp_tissues[group], active_regions.columns)
		else:
			existing_tissues = active_regions.columns

                if len(existing_tissues) == 0:
                        continue
                peak_annotated = active_regions[existing_tissues]
                peak_annotated = np.array(peak_annotated)

                if 'DNase' in featureName:
			Active = [float(xx.split(',')[0]) > 0 for x in peak_annotated for xx in x]
                elif 'ROADMAP' in featureName:
			Active = [sub in xx for x in peak_annotated for xx in x]
		Active = np.reshape(Active, peak_annotated.shape)
		Active = np.sum(Active, axis=1) > (len(existing_tissues) * tissue_prop)
		active_snps = np.array(active_regions.index[np.where(Active)[0]])

		outfile = open('%s/%s_Active_SNPset_%s_%s_group%d.txt' % (activeSNPdir, FMfn, featureName, sub.replace('/','_'), group), 'w')
		for snp in active_snps:
			outfile.write('%s\n' % snp)
		outfile.close()




if __name__ == '__main__':
	alldatasetName = 'SNP_loc'
	tissue_prop = 0.49

	featureName = 'DNase'
	active_regions = pd.read_csv('%s/tissue_%s/%s_%s_onlyhits.bed' % (allSNPfeaturedir, featureName, alldatasetName, featureName), sep='\t', index_col = 0)
	readInGenomicAnnotation_factors(active_regions, featureName, 'NONE')

	featureName = 'ROADMAP'
	active_regions = pd.read_csv('%s/tissue_%s/%s_%s_onlyhits.bed' % (allSNPfeaturedir, featureName, alldatasetName, featureName), sep='\t', index_col = 0)
	for sub in ['1_TssA', '6_EnhG', '7_Enh', '4_Tx', '5_TxWk']:
		readInGenomicAnnotation_factors(active_regions, featureName, sub)


