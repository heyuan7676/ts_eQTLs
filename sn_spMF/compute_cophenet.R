##########################################################################
## Compute cophenetic coefficient for factors from multiple runs
## See Brunet et al. (2004)
##########################################################################

library(plyr)

compute_cophenet <- function(M){
	#compute the consensus matrix C, allowing samples to be assigned to one factor
	C = list();
	D = dim(M[[1]])[1]

	for(m in 1:length(M)){
		assignment = apply(M[[m]], 1, which.max)
		ag_rep_rows = matrix(rep(assignment, D), ncol=D)
		C[[m]] = (ag_rep_rows == t(ag_rep_rows)) * 1;
	}

	#take the average to get the fraction of consistent assignment
	Cbar = aaply(laply(C, as.matrix), c(2, 3), mean)

	#compute the correlation between original distances and the cophenetic distances from a hierachical clustering based on average linkage
	Ds = Cbar - diag(rep(1,D));
	d1 = dist(Ds)
	hc = hclust(d1, "ave")
	d2 = cophenetic(hc)
	coph = cor(d1, d2)

	return(coph)
}
