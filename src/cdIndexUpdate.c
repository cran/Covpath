

/*
 *  cdIndexUpdatea.c
 *  Covpath R Package
 *  Created by Vijay Krishnamurthy on 3/25/11.
 */

#include "bcd.h"
void cdIndexUpdate(double *adVinmat, double *adbvec, double dscalarc, double t,  double dRho, double *aduo, double dscalarwo, double *aduvec, double *adw, int n, double dTolforIntitalGap)
{
	// Local Variables 
	double 	 dTemp, dTemp2, dTemp3, temp, dTemp1, p1, p2, p3, p4, adRoots[3];
	double  *buffveca, one = 1.0, zero = 0, dalpha, dbeta, dgamma;
	int j,i, incx =1, row1, row2, col1, col2,iNumRoots;
	double ady[n], adz[n], pobj, dobj, x1, x2, x3, dOrigGap;
	buffveca = calloc(n, sizeof(double));
	cblasR_dcopy(n,aduo,incx,aduvec,incx);
	
	*adw = dscalarwo;
	// Compute initial Duality Gap
	// x1 = w - u'T*V*u
	row1 = n; col1 = n; row2 = n; col2 = 1;
	cblasR_dgemv(CblasColMajor,CblasNoTrans, row1, col1, one, adVinmat, row1, aduvec, incx, zero, buffveca,incx);
	dTemp = -(cblasR_ddot(n, aduvec,incx, buffveca, incx)); 
	dTemp = dTemp + *adw;
	x1 = dTemp;x2 = dRho + dscalarc - *adw;x3 = dRho - dscalarc + *adw;
	for(i =0; i < n; i++)
	{
		ady[i] = (dRho - aduvec[i] + adbvec[i]);
		adz[i] = (dRho + aduvec[i] - adbvec[i]);
	}
	pobj = 	bcdPrimal(x1, x2, x3, ady, adz, t, n); dobj = bcdDual( x2,  x3, ady, adz, t, dRho, adbvec, dscalarc,n );
	dOrigGap = pobj - dobj;
	
	if ( dOrigGap > dTolforIntitalGap)
	{
		for(j = 0; j < n ; j++)
		{
			//dTemp = w - u'T*V*u
			row1 = n; col1 = n; row2 = n; col2 = 1;
			cblasR_dgemv(CblasColMajor,CblasNoTrans, row1, col1, one, adVinmat, row1, aduvec, incx, zero, buffveca,incx);
			dTemp = -(cblasR_ddot(n, aduvec,incx, buffveca, incx)); 
			dTemp = dTemp + *adw;
			dTemp1 = aduvec[j];
			dTemp2 = adbvec[j];
			//alpha = - Vin(j,j);
			dalpha = -adVinmat[j + n*j];
			//beta = -2*(Vin(j,:)*u - Vin(j,j)*u(j));
			cblasR_dcopy(n,(adVinmat + (n*j)),incx,buffveca,incx);
			dTemp3 = (cblasR_ddot(n, aduvec,incx, buffveca, incx));
			dbeta = -2*(dTemp3 + (dalpha*dTemp1));
			//gamma = w - u'*Vin*u  -alpha*u(j)*u(j) -beta*u(j);
			dgamma = dTemp - (dalpha*dTemp1*dTemp1) -(dbeta*dTemp1);
			p1 = 2*dalpha + 4*t*dalpha;
			p2 = -4*dalpha*dTemp2 + dbeta - 4*t*dalpha*dTemp2 + 4*t*dbeta; 
			p3 = -2*dRho*dRho*dalpha + 2*dalpha*dTemp2*dTemp2 -2*dTemp2*dbeta -4*t*dbeta*dTemp2 + 4*t*dgamma;
			p4 = -dbeta*dRho*dRho + dbeta*dTemp2*dTemp2 -4*t*dTemp2*dgamma; 
			iNumRoots =  rootsofCubic(p1,p2,p3,p4, adRoots);
			aduvec[j] = minucdobj(adRoots, iNumRoots,dalpha,dbeta,dgamma,t, dTemp2,dRho);
		}	
		
		// Compute diagonal entry separately
		//dTemp =  u'T*V*u, 
		row1 = n; col1 = n; row2 = n; col2 = 1;
		cblasR_dgemv(CblasColMajor,CblasNoTrans, row1, col1, one, adVinmat, row1, aduvec, incx, zero, buffveca,incx);
		dTemp = (cblasR_ddot(n, aduvec,incx, buffveca, incx)); 
		p1 = 1 + (2*t);
		p2 = -2*(dscalarc*t + dTemp*t + dscalarc);
		p3 = (dscalarc*dscalarc) + 2*t*dTemp*dscalarc - (dRho*dRho);
		temp = dTemp;
		dTemp = (p2*p2) - 4*p1*p3;
		if ( dTemp >= 0)
		{
			dTemp1  = (-p2 + sqrt(dTemp))/(2*p1);
			dTemp2  = (-p2 - sqrt(dTemp))/(2*p1);
			*adw = minwcdobj(dTemp1, dTemp2, temp, dscalarc, t, dRho);
		}
	}
	free(buffveca);
}

