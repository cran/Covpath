thresholdCondCorrMatrix <- function(C,thrsQuantile){
	
	pnnz = thrsQuantile/100;
	thrc=as.vector(quantile(c(abs(C)),1-pnnz))
	coefs = C*(abs(C) > thrc)
	m2 = abs(coefs)
	m2[lower.tri(m2,diag =TRUE)] = 0
	temp = colSums(m2) + colSums(t(m2))
	ikeep = which(temp>0)
	X3 = coefs[ikeep,ikeep]
	soln = list(ikeep = ikeep , prunedC = X3)
	return(soln)
	
	}