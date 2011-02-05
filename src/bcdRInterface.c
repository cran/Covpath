

/*
 *  bcdRInterface.c
 *  Covpath R Package
 *  Created by Vijay Krishnamurthy on 3/25/11.
 */


#include <R.h>
#include <Rinternals.h>
#include <R.h> 
#include "bcd.h"

double *cmpvec;
SEXP bcdCorrector(SEXP Sigma, SEXP Uk, SEXP rho, SEXP t, SEXP tol, SEXP colFraction, 
		  SEXP Xk, SEXP origGap) 

{
  // Variables to accept R Inputs
  double *adUmat, *adSigmat, *adUk, *adXk, *dRelGap, dscalart,dRho, dTol, *adXmat, dColFraction,  dOrigGap; 
  int i,n,ibcdOuterIter,iIndexLength,  incx = 1, *dim;
  int *adIndexList,*adRowIndexList; // wrong name convention 
  double *prho, *pt, *ptol, *pColFr, *pOrigGap;
  SEXP Uknew, Xknew,relGapnew, list, list_names;
  char *names[3] = {"Uk", "Xk", "relGap"};
  
	
  /* Read R Inputs */
  adSigmat = REAL(Sigma);
  adUk = REAL(Uk);
  prho  = REAL(rho);dRho = *prho;
  pt = REAL(t);dscalart  = *pt;
  ptol=  REAL(tol);dTol = *ptol;
  pColFr =  REAL(colFraction);dColFraction = *pColFr;
  adXk = REAL(Xk);
  pOrigGap =  REAL(origGap);dOrigGap= *pOrigGap;
  ibcdOuterIter = 5;
  dim =  INTEGER(coerceVector(getAttrib(Sigma, R_DimSymbol), INTSXP));
  n = dim[0];
	
  // outputs 
  PROTECT(Uknew =allocMatrix(REALSXP, n ,n));
  PROTECT(Xknew = allocMatrix(REALSXP, n ,n));
  PROTECT(relGapnew = allocVector(REALSXP,1));
  PROTECT(list_names = allocVector(STRSXP,3));

  for(i = 0; i < 3; i++)   
    SET_STRING_ELT(list_names,i,mkChar(names[i])); 

  adUmat = REAL(Uknew);
  adXmat = REAL(Xknew);
  dRelGap = REAL(relGapnew) ;
	
  cblasR_dcopy((n*n),adUk,incx,adUmat,incx);
  cblasR_dcopy((n*n),adXk,incx,adXmat,incx);

	
  //Memory Allocations
  cmpvec=(double *) calloc(n+1,sizeof(double));
  adIndexList = (int*)calloc(n, sizeof(double));
  iIndexLength = (int)(floor(n*dColFraction));
  adRowIndexList = (int*)calloc(n, sizeof(double));
  //iRowIndexLength = (int)(floor((n-1)*dRowFraction));
    

  (*dRelGap) = coordDescent(adSigmat, adUmat, dRho,dscalart, dTol,adIndexList, iIndexLength,n, adXmat,dOrigGap,ibcdOuterIter);
  free(adIndexList);
  free(adRowIndexList);
  free(cmpvec);
  
  PROTECT(list = allocVector(VECSXP, 3)); 
  SET_VECTOR_ELT(list, 0, Uknew); 
  SET_VECTOR_ELT(list, 1, Xknew); 
  SET_VECTOR_ELT(list, 2, relGapnew); 

  // and attaching the vector names:
  setAttrib(list, R_NamesSymbol, list_names); 
  UNPROTECT(5);
  return list;


	
}	
	
	
                
  



