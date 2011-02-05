/*
 *  bcd.h
 *  Covpath R Package
 *  Created by Vijay Krishnamurthy on 3/25/11.
 */


#ifndef HEADERFORCOV
#define HEADERFORCOV
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <string.h>
#include <limits.h>
#include <complex.h>
#include <R_ext/Lapack.h>

enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102 };
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113,
	AtlasConj=114};
enum CBLAS_UPLO  {CblasUpper=121, CblasLower=122};
enum CBLAS_DIAG  {CblasNonUnit=131, CblasUnit=132};
enum CBLAS_SIDE  {CblasLeft=141, CblasRight=142};
#endif

double cblasR_dasum(int n,  double *x, int incx1);
double cblasR_ddot( int n, double *xvec,int incx1, double *yvec, int incx2);
void cblasR_dscal( int N, double alpha, double *X, int incX);
void cblasR_dcopy(int N,double *X,int incX,double *Y,int incY);
void cblasR_dgemm(enum CBLAS_ORDER Order,enum CBLAS_TRANSPOSE transA, enum CBLAS_TRANSPOSE transB,int M, int N, int K, double alpha, double *A, int lda,double *B, int ldb, double beta, double *C, int ldc);
void cblasR_dgemv(enum CBLAS_ORDER Order,enum CBLAS_TRANSPOSE transA, int M, int N, double alpha, double *A, int lda,double *B, int incB, double beta, double *C, int incC);
void cblasR_daxpy(int N,double alpha,double *X,int incX,double *Y,int incY);




// Local Functions
void modsort(double *numbers, int *index, int array_size);
int cmp(const void *keyvali, const void *datumi);

double coordDescent(double *adSigmat, double *adUmat, double dRho,double dscalart, double dTol, int *adIndexList, int iIndexLength, int n, double *adXmat, double dOrigGap, int ibcdOuterIter);

//void cdIndexUpdate(double *adVinmat, double *adbvec, double dscalarc, double t,  double dRho, double *aduo, double dscalarwo, double *aduvec, double *adw, int n, double dTolforIntitalGap, int *adRowIndices, int iIndexLength);
void cdIndexUpdate(double *adVinmat, double *adbvec, double dscalarc, double t,  double dRho, double *aduo, double dscalarwo, double *aduvec, double *adw, int n, double dTolforIntitalGap);

int rootsofCubic(double dscalara, double dscalarb, double dscalarc, double dscalard, double *adRoots);
double minwcdobj(double dTemp1, double dTemp2, double temp, double dscalarc, double t, double dRho);
double minucdobj(double *adRoots, int iNumRoots,double dalpha,double dbeta,double dgamma,double t, double b, double dRho);



void projectFeas(double *aduo, double dscalarwo, double dRho, double *adbvec, double dscalarc, double *adVinmat, double *aduvec, double *adw, double *buffveca, int n);
int newtonForCD(double *adVinmat, double *adbvec, double dscalarc, double t,  double dRho, double *aduo, double dscalarwo, double dtol, int iMaxIters, int iPeriod, double dalpha, double dbeta, double dTolforIntitalGap, double *aduvec, double *adw, int n);
double mlog(double x);
double bcdDual(double x2,double x3, double *y,double *z,double t, double rho, double *b, double c, int n );
double bcdPrimal(double x1,double x2,double x3, double *y,double *z,double t,int n);


void real_inverse(double *prx,int nrows,double *pry);
void real_lineareq(double *pra, double *prb, int n, double *prx);	
void planerot( double *G, double *y, double *x);
	
 
// Local functions 


int iminif(int x, int y);

double doubsum(double *xmat, int n);

double doubdot(double *xvec,  double *yvec, int n);

double doubasum(double *xmat, int n);

double doubnorm2(double *xmat, int n);

double infnorm(double *xmat, int n);

double frobnorm(double *xmat, int n);

int idxmax(double *xmat, int n);

int idxabsmax(double *xmat, int n);

double dsignf(double x);

double dminif(double x, double y);

double dmaxf(double x, double y);
int iminf(int x, int y);
int imaxf(int x, int y);

double dabsf(double x);

int isignf(double x);

void dispmat(double *xmat, int n, int m);

void symmetrize(double *xmat,double *ymat,int n);

double maxeig(double *xmat,double *bufveca,double *bufvecb,int n);


