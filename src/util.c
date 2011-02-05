
/*
 *  util.c
 *  Covpath R Package
 *  Created by Vijay Krishnamurthy on 3/25/11.
 */




#include "bcd.h"


double doubsum(double *xmat, int n)
{
	int i;
	double res=0.0;
	for (i=0;i<n;i++){res+=xmat[i];};
	return res;
}

double doubdot(double *xvec, double *yvec, int n)
{
	int i;
	double res=0.0;
	for (i=0;i<n;i++){res+=xvec[i]*yvec[i];};
	return res;
}


double doubasum(double *xmat, int n)
{
	int i;
	double res=0.0;
	for (i=0;i<n;i++){res+=dabsf(xmat[i]);};
	return res;
}
 

double cblasR_dasum(int n,  double *x, int incx1)
{
	double	res;
	res  = doubasum(x,n);
	return res;

}
double cblasR_ddot( int n, double *xvec,int incx1, double *yvec, int incx2)
{
	double dp;
	dp = doubdot(xvec, yvec, n);
	return dp;
}
void cblasR_dscal( int N, double alpha, double *X, int incX)
{
	dscal_(&N,&alpha,X,&incX);
}

void cblasR_dcopy(int N,double *X,int incX,double *Y,int incY)
{
	dcopy_(&N,X,&incX,Y,&incY);
}

void cblasR_dgemm(enum CBLAS_ORDER Order,enum CBLAS_TRANSPOSE transA, enum CBLAS_TRANSPOSE transB,
				  int M, int N, int K, double alpha, double *A, int lda,
				  double *B, int ldb, double beta, double *C, int ldc)
{
	char ta[1],tb[1];
	if (transA==111)
	{
		*ta='N';
	}
	else
	{
		*ta='T';
	};
	if (transB==111) 
	{
		*tb='N';
	}
	else
	{
		*tb='T';
	};
	dgemm_(ta,tb,&M,&N,&K,&alpha,A,&lda,B,&ldb,&beta,C,&ldc);
}

void cblasR_dgemv(enum CBLAS_ORDER Order,enum CBLAS_TRANSPOSE transA,
				  int M, int N, double alpha, double *A, int lda,
				  double *B, int incB, double beta, double *C, int incC)
{
	char ta[1];
	if (transA==111)
	{
		*ta='N';
	}
	else
	{
		*ta='T';
	};
	dgemv_(ta,&M,&N,&alpha,A,&lda,B,&incB,&beta,C,&incC);
}

void cblasR_daxpy(int N,double alpha,double *X,int incX,double *Y,int incY)
{
	daxpy_(&N,&alpha,X,&incX,Y,&incY);
}



void real_lineareq(double *pra, double *prb, int n, double *prx)
{
	
	
	// To compute solution to linear equation 
	// Ax = b ; 
	// on input pra is pointer to array of size n*n
	//			prb is pointer to array of size n
	//			prx is pointer to empty array of size n
	// on output pra, prb is unchanged 
	//			 prx is pointer to the solution
	int *ipivot,  incx = 1, lapackinfo, colb; 
	
	
	ipivot = (int *)malloc(n*sizeof(int));
	colb = 1;
	cblasR_dcopy(n,prb,incx,prx,incx);
	dgesv_(&n, &colb, pra, &n, ipivot, prx, &n, &lapackinfo);  
	
	
	
}



void real_inverse(double *prx,int nrows,double *pry)
{
	
	
	// To compute Matrix Inverse. 
	
	// y = inv(x)
	// Make sure parameters are of right size
	// on input pry is pointer to a empty array of size nrows*nrows 
	//			prx is pointer to array of size nrows*nrows	
	// on output prx is unchanged 
	//			 pry is pointer to inverse of prx. 
	
	
	
	double  *temp, *work;
	int *ipivot, iwork,	k, info;
	
	iwork = nrows;
	work = (double *)malloc(nrows*sizeof(double));
	ipivot = (int *)malloc(nrows*sizeof(int));
	
	temp = prx; 
	for (k = 0; k < nrows*nrows; k++) { 
		pry[k] = temp[k]; 
	}
	
	dgetrf_(&nrows, &nrows, pry, &nrows, ipivot, &info);
	dgetri_(&nrows, pry,  &nrows, ipivot,work, &iwork, &info);
	
	
	free(ipivot);
	free(work);
	
}




// Sorting function, returns indices and sorted vector
void modsort(double *numbers, int *index, int array_size)
{
	int incx=1;
	extern double *cmpvec;
	
	cblasR_dcopy(array_size,numbers,incx,cmpvec,incx);
	qsort(index,(size_t) array_size, sizeof(int),cmp);
}

int cmp(const void *keyvali, const void *datumi)
{
	extern double *cmpvec;
	int* keyval = (int*) keyvali;
	int* datum = (int*) datumi;
	
	if (cmpvec[*keyval]==cmpvec[*datum])
	{
		return 0;
	}
	else
	{
		if (cmpvec[*keyval]<cmpvec[*datum]){return 1;}
		else {return -1;}
	}
	
}




// mlog 
double mlog(double x)
{
	
	if ( x <= 0)
	{
		return(-INFINITY);
	}
	else
		return(log(x));
}



void planerot( double *G, double *y, double *x)
{
	//       [G,Y] = PLANEROT(X), where X is a 2-component column vector,
	//       returns a 2-by-2 orthogonal matrix G so that Y = G*X has Y(2) = 0.
	//       See also QRINSERT, QRDELETE.
	
	double r; 
	if (x[1] != 0)
    {
		r  = doubnorm2(x, 2);
		G[0] = x[0]/r;
		G[1] = -x[1]/r;
		G[2] = x[1]/r;
		G[3] = x[0]/r;
		y[0] = r;
		y[1] = 0;
    }
	else
    {
		G[0] = 1.0;
		G[1] = 0.0;
		G[2] = 0.0;
		G[3] = 1.0;
		y[0] = x[0];
		y[1] = x[1];
    }
	
}




int isignf(double x)
{
	if (x>=0)
		return 1;
	else
		return -1;
}




// Some useful functions ...

int idxmax(double *xmat, int n)
{
	int i;
	int res=0;
	for (i=0;i<n;i++)
	{
		if (xmat[i]>xmat[res]) {res=i;}
	}
	return res;
}



int idxabsmax(double *xmat, int n)
{
	int i;
	int res=0;
	for (i=0;i<n;i++)
	{
		if ( dabsf(xmat[i]) > dabsf(xmat[res]) ) {res=i;}
	}
	return res;
}



double doubnorm2(double *xmat, int n)
{
	int i;
	double res=0.0;
	for (i=0;i<n;i++){res+=xmat[i]*xmat[i];};
	return sqrt(res);
}

double infnorm(double *xmat, int n)
{
	int i,j;
	double res=0.0,sum;
	
	for (j=0;j<n;j++){
		sum=0.0;
		for(i=0;i<n;i++)
			sum+=dabsf(xmat[j+i*n]);
		if(sum>=res) res=sum;
	}	
	return res;
}

double frobnorm(double *xmat, int n)
{
	int i,j;
	double res=0.0;
	
	for (i=0;i<n;i++)
		for (j=0;j<n;j++)
			res+=(xmat[i*n+j]*xmat[i*n+j]);
	
	return pow(res,.5);
}

double dsignf(double x)
{
	if (x>=0)
		return 1.0;
	else
		return -1.0;
}

double dminif(double x, double y)
{
	if (x>=y)
		return y;
	else
		return x;
}

int iminif(int x, int y)
{
	if (x>=y)
		return y;
	else
		return x;
}


double dmaxf(double x, double y)
{
	if (x>=y)
		return x;
	else
		return y;
}

int imaxf(int x, int y)
{
	if (x>=y)
		return x;
	else
		return y;
}

double dabsf(double x)
{
	if (x>=0)
		return x;
	else
		return -x;
}

int iroundf(double number)
{
    return (number >= 0) ? (int)(number + 0.5) : (int)(number - 0.5);
}

void dispmat(double *xmat, int n, int m)
{
	int i,j;
	
	for (i=0; i<n; i++)
	{
		for (j=0;j<m;j++)
		{
			printf("%+.4f ",xmat[j*n+i]);
		}
		printf("\n");
	}
	printf("\n");
}

//

// symmetrize a matrix xmat and return in matrix ymat
void symmetrize(double *xmat,double *ymat,int n)
{
	int i,j,incx=1,n2=n*n;
	double alpha=0.0;
	
	cblasR_dscal(n2,alpha,ymat,incx);	
	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			ymat[n*i+j]=(xmat[n*i+j]+xmat[n*j+i])/2;	
}

 
