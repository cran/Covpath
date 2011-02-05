
makeCondMat <- function(Umat){

n = ncol(Umat)
condmat = matrix(0, nrow = n, ncol =n)
#condmat = rep(0,n*n)
for ( i in 1:n){
	for (j in (i+1):n){
		if ( i < n ){
				if (abs(Umat[i,j]) > 0){
				#	if (abs(Umat[i+ (j-1)*n]) > 0){
			#bufm = solve(Umat[c(i,j),c(i,j)])
			a = Umat[i,i]
			b = Umat[i,j]
			c = Umat[j,i]
			d = Umat[j,j]	
			#a = Umat[i +  (i-1)*n]
			#b = Umat[i+ (j-1)*n]
			#c = Umat[j + (i-1)*n]
			#d = Umat[j+(j-1)*n]
			e = a*d - b*c
			b11 = d/e
			b22 = a/e
			b12 = -b/e
			cc = b12/(sqrt(b11*b22))
			#condmat[i,j] = bufm[1,2]/(sqrt(bufm[1,1]*bufm[2,2]))
			#condmat[j,i] = bufm[1,2]/(sqrt(bufm[1,1]*bufm[2,2]))
			condmat[i,j] = cc
			condmat[j,i] = cc
			#condmat[j + (i-1)*n] = cc
			#condmat[i + (j-1)*n] = cc
			
}}
}
}
		
		return(condmat)
		}