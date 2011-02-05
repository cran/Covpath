
/*
 *  memberfunc.c
 *  Covpath R Package
 *  Created by Vijay Krishnamurthy on 3/25/11.
 */


#include "bcd.h"



int rootsofCubic(double dscalara, double dscalarb, double dscalarc, double dscalard, double *adRoots)
{
	// Computes the roots of a cubic equation 
	
	double p, q, Q, R,D, S, T, theta, dTemp1, c0, c1, c2, c3;
	int iNumDistinctRealRoots;
	c0 = dscalard/dscalara;
	c1 = dscalarc/dscalara; 
	c2 = dscalarb/dscalara;
	c3 = 1.0;
	p = (3*c1-(pow(c2,2)))/3;
	q = (9*c1*c2-27*c0-2*(pow(c2,3)))/27;
	Q=p/3;R=q/2;D=(pow(Q,3)) +(pow(R,2));
	if ( D > 0)
	{
		// Only one root real 
		iNumDistinctRealRoots = 1;
		S = cbrt(R+sqrt(D));
		T = cbrt(R -sqrt(D));
		adRoots[0] = -(c2/3)+(S+T);
		return iNumDistinctRealRoots; 
	}
	else if (D == 0 )
	{
		// All  Roots Real and atmost 2 distinct roots
		S = cbrt(R+sqrt(D));
		T = cbrt(R -sqrt(D));
		iNumDistinctRealRoots = 2;
		adRoots[0] = -(c2/3)+(S+T);
		adRoots[1] = -(c2/3)-((S+T)/2);
		return iNumDistinctRealRoots; 
	}
	else
	{
		// All real roots and distinct and unequal 
		iNumDistinctRealRoots = 3; 
		dTemp1 = R/(sqrt(-(pow(Q,3))));
		theta = acos(dTemp1);
		adRoots[0] = 2*sqrt(-Q)*(cos((theta/3))) - (c2/3);
		adRoots[1] = 2*sqrt(-Q)*(cos((theta+(2*M_PI))/3)) - (c2/3);
		adRoots[2] =  2*sqrt(-Q)*(cos((theta+(4*M_PI))/3)) - (c2/3);
		return iNumDistinctRealRoots; 
	}
}


double minucdobj(double *adRoots, int iNumRoots,double dalpha,double dbeta,double dgamma,double t, double b, double dRho)
{
	
	// Computes the update for a non diagonal entry by comparing objectives corresponding to the roots of the cubic equation
	int i;
	double dTemp[3];
	
	if ( iNumRoots == 1 )
		return adRoots[0];
	else if (iNumRoots ==2)
	{
		dTemp[0] = -mlog(dalpha*adRoots[0]*adRoots[0] + dbeta*adRoots[0] + dgamma) - 2*t*mlog(dRho+adRoots[0]-b) - 2*t*mlog(dRho-adRoots[0]+b);
		dTemp[1]= -mlog(dalpha*adRoots[1]*adRoots[1] + dbeta*adRoots[1] + dgamma) - 2*t*mlog(dRho+adRoots[1]-b) - 2*t*mlog(dRho-adRoots[1]+b);
		return(dTemp[0] <= dTemp[1] ?   adRoots[0]  :   adRoots[1]);
	}
	else
	{
		dTemp[0] = -mlog(dalpha*adRoots[0]*adRoots[0] + dbeta*adRoots[0] + dgamma) - 2*t*mlog(dRho+adRoots[0]-b) - 2*t*mlog(dRho-adRoots[0]+b);
		dTemp[1]= -mlog(dalpha*adRoots[1]*adRoots[1] + dbeta*adRoots[1] + dgamma) - 2*t*mlog(dRho+adRoots[1]-b) - 2*t*mlog(dRho-adRoots[1]+b);
	    i = (dTemp[0] <= dTemp[1]) ? 0 : 1;
	    dTemp[2]= -mlog(dalpha*adRoots[2]*adRoots[2] + dbeta*adRoots[2] + dgamma) - 2*t*mlog(dRho+adRoots[2]-b) - 2*t*mlog(dRho-adRoots[2]+b);
		return (dTemp[i] <= dTemp[2] ? adRoots[i]  : adRoots[2]) ;
	}
}




double minwcdobj(double dTemp1, double dTemp2, double temp, double dscalarc, double t, double dRho)
{
	// Computes the update for a diagonal entry by comparing objectives corresponding to the roots of the quadratic  equation

	double dTemp[2];
	dTemp[0] = -mlog(dTemp1-temp) - t*mlog(dRho+dTemp1-dscalarc) - t*mlog(dRho-dTemp1+dscalarc);
	dTemp[1] = -mlog(dTemp2-temp) - t*mlog(dRho+dTemp2-dscalarc) - t*mlog(dRho-dTemp2+dscalarc);
	return (dTemp[0] <= dTemp[1] ? dTemp1 : dTemp2 );
}



double bcdDual(double x2,double x3, double *y,double *z,double t, double rho, double *b, double c, int n )
{
	int nOrig, i;
	double alpha1, alpha2,alpha3, dobj;
	nOrig = n+1;
	alpha2 = t/x2;
	alpha3 = t/x3;
	alpha1 = alpha2 - alpha3; 
	dobj = 1+2*t*(2*nOrig -1) + mlog(alpha1) - alpha2*(rho+c) - alpha3*(rho -c) + t*mlog(alpha2/t) + t*mlog(alpha3/t);
	for(i = 0 ; i < n; i++)
	{
		dobj = dobj - 2*t*(rho+b[i])/y[i] - 2*t*(rho-b[i])/z[i] -2*t*log(y[i]) -2*t*log(z[i]);
	}
	return dobj;
}

double bcdPrimal(double x1,double x2,double x3, double *y,double *z,double t,int n)
{
	int i;
	double pobj = 0.0;
	pobj = -mlog(x1) - t*mlog(x2) - t*mlog(x3);
	for(i = 0; i < n ; i++)
	{
		pobj -= 2*t*(mlog(y[i]) + mlog(z[i]));
	}
	return pobj;
}



