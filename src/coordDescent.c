
/*
 *  coordDescent.c
 *  Covpath R Package
 *  Created by Vijay Krishnamurthy on 3/25/11.
 */


#include "bcd.h"



double coordDescent(double *adSigmat, double *adUmat, double dRho,double dscalart, double dTol, int *adIndexList, int iIndexLength, int n, double *adXmat, double dOrigGap, int ibcdOuterIter)
{
	int i,j, k,row1, col1, row2, col2, iTemp, row, iCol, iVal, incx = 1, iCounter;
 	double one = 1, zero = 0, *buffmat, *adevec, *adhvec, sNegOne = -1, dTemp, dTemp2, *adbvec, dscalarc, *aduo, *adwo, *aduvec, *adw, *adVinmat, *buffveca, *adfvec, *buffvecb, *adTempVec,dCurrGap, dRelGap; 

	adTempVec = calloc(n, sizeof(double));
	adevec = calloc(n, sizeof(double));
	adhvec = calloc(n, sizeof(double));
	adfvec = calloc(n, sizeof(double));
	adVinmat = calloc( (n-1)*(n-1) , sizeof(double));
	aduo = calloc((n-1), sizeof(double));
	adwo = calloc(1, sizeof(double));
	aduvec = calloc((n-1), sizeof(double));
	adw = calloc(1, sizeof(double));
	adbvec = calloc((n-1), sizeof(double));
	buffveca = calloc(n, sizeof(double));
	buffvecb = calloc(n, sizeof(double));
	buffmat = calloc (n*n, sizeof(double));
	k =0;
	/* DualityGap =  -n + trace(Sigma*Xk) + rho*(sum(sum(abs(Xk))));*/
		
	dCurrGap = - n + cblasR_ddot((n*n), adSigmat,incx, adXmat, incx) + dRho*(cblasR_dasum((n*n), adXmat,incx));
	dRelGap = dCurrGap/dOrigGap;
	
	/* Main loop: Iterate until converge or till a maximum bcdOuterIter iterations */
	iCounter = 0; 
	while ( (dRelGap > dTol) && (iCounter < ibcdOuterIter )){
		iCounter++;
		/* 
		 Choose Important Columns to Update.
		 ind = indfind(Uk,Xk,Sigma, rhonew,colFraction );
		 */
		for (i=0;i<n;i++){	
			adTempVec[i]= (dRho + adSigmat[i + (n*i)] - adUmat[i + (n*i)])*adXmat[i+(n*i)]; 
			adIndexList[i]=i;
		}  
		modsort(adTempVec,adIndexList,n);
		/* Update choosen columns */  
		i = 0; 
		while  (i < iIndexLength){	
			iCol = adIndexList[i]; 
			/* 
			 1. Form Vin
			 First Rank 1 update ; 
			 e(i) = 1 ; e(j) = 0 (j != i) 
			 f = -U(:,i); f(i) = 1 - U(i,i)
			 Vin=Xk-((Xk*f*e'*Xk)/(1+e'*Xk*'f);
			 */
			adevec[iCol] = 1;
			cblasR_daxpy( n , sNegOne, (adUmat + (n*iCol)) ,incx,adfvec,incx);
			cblasR_daxpy( n , sNegOne, (adUmat + (n*iCol)) ,incx,adhvec,incx);
			adfvec[iCol] = 1 - adUmat[iCol + n*iCol];
			adhvec[iCol] = 0.0;
			for(iVal = 0 ; iVal < n; iVal++)
			{
				buffveca[iVal] = adXmat[iCol + n*iVal];
			}
			dTemp = 1 +  cblasR_ddot(n, adfvec,incx, buffveca, incx);
			cblasR_dscal(n,(1/dTemp),adfvec,incx);
			row1 = n; col1 = n; row2= n; col2=1;
			cblasR_dgemv( CblasColMajor,CblasNoTrans, n, n, one, adXmat, n, adfvec ,incx, zero, buffveca, incx); // buffveca contains Xk*f
			for(j = 0 ; j < n ; j++){
				if ( j != iCol){
					if ( j < iCol)
						k = j;
					else if ( j > iCol)
						k = j-1;
					cblasR_dcopy( iCol , (adXmat + (n*j)) ,incx,(adVinmat + ((n-1)*(k))),incx);
					cblasR_dcopy( n-iCol-1 , (adXmat + (n*j + iCol+1)) ,incx,(adVinmat + (n-1)*k + iCol),incx);
					dTemp = (adXmat[iCol + n*j]);
					cblasR_daxpy( iCol , -dTemp, buffveca   ,incx,(adVinmat + ((n-1)*(k))),incx);
					cblasR_daxpy( n-iCol-1 , -dTemp, (buffveca + (iCol+1)) ,incx,(adVinmat + (n-1)*k + iCol),incx);
				}
				buffvecb[j] = adXmat[iCol + n*j] *(1 - buffveca[iCol]);
			}
			/*
			 Second Rank1 update 
			 ;h = -U(:,i); h(i) = 0;
			 Vin=Vin-((Vin*e*h'*Vin)/(1+e'*Vin*'h));Vin(i,:)=[];Vin(:,i)=[];
			 */
			dTemp2 = 1 +  cblasR_ddot(n, adhvec,incx, buffvecb, incx);
			cblasR_dscal(n,(1/dTemp2),adhvec,incx);
			cblasR_dscal( n ,-adXmat[iCol + n*iCol] ,buffveca,incx);
			cblasR_daxpy( n , one, (adXmat + (n*iCol)) ,incx,buffveca,incx);
			for ( j = 0 ; j < n ; j++){
				if ( j != iCol){
					if ( j < iCol)
						k = j;
					else if ( j > iCol)
						k = j-1;
					dTemp = cblasR_ddot(iCol, adhvec,incx,(adVinmat + ((n-1)*(k))),incx);
					dTemp = dTemp + cblasR_ddot(n-iCol-1, (adhvec + iCol +1),incx,(adVinmat + (n-1)*k + iCol),incx);
					cblasR_daxpy( iCol , -dTemp, buffveca   ,incx,(adVinmat + ((n-1)*(k))),incx);
					cblasR_daxpy( n-iCol-1 , -dTemp, (buffveca + (iCol+1)) ,incx,(adVinmat + (n-1)*k + iCol),incx);
				}
			}
			/*
			 2. Form b, c, uo, wo - inputs to problem of updating choosen column, 
			 then call cdIndexUpdate
			 b = Sigma(:,i); c = b(i); b(i) = [];
			 uo = Uk(indi,i);wo = Uk(i,i);
			 [u,w] = indexCD_mex(Vin,b,c,t,rho,uo,wo,tolForInitGap);
			 FOR LATER IF SPEED IS AN ISSUE- Code heuristic to choose rows within the chosen col -. Also change 4 accordingly 
			 */
			//AA: can elimate dcopy calls by using pointer arthmetic 
			// KV - I definitely need copies. 
			cblasR_dcopy( iCol , (adUmat + (n*iCol))   ,incx, aduo,incx);
			cblasR_dcopy( n-iCol-1 , (adUmat + (n*iCol) + (iCol+1)) ,incx,(aduo + iCol),incx);
			cblasR_dcopy( iCol , (adSigmat + (n*iCol))   ,incx, adbvec,incx);
			cblasR_dcopy( n-iCol-1 , (adSigmat + (n*iCol) + (iCol+1)) ,incx,(adbvec + iCol),incx);
			dscalarc = *((adSigmat + (n*iCol) + iCol));
			*adwo =  *((adUmat + (n*iCol) + iCol));
			cdIndexUpdate(adVinmat, adbvec,dscalarc, dscalart, dRho, aduo, (*adwo), aduvec, adw,(n-1),dTol);
			/*
			 3. Update X = inv(U) by application of SMW - Doing it before update U 
			 First Rank one update 
			 indi = 1:n;indi(indi==i)=[];
			 f= zeros(n,1); h = zeros(n,1);f(indi)=u-U(indi,i);f(i)= (w - U(i,i); h(indi)=u-U(indi,i);h(i)= 0;	
			 Xk = Xk - ((Xk*f*e'*Xk)/(e'*X*f));
			 */
			for ( j = 0 ; j <iCol ; j++){
				adhvec[j] = aduvec[j] - adUmat[(n*iCol) + j];
				adfvec[j] = aduvec[j] - adUmat[(n*iCol) + j];
			}
			for(j = iCol+1 ; j < n ; j++){
				adhvec[j] = aduvec[(j-1)] -  adUmat[(n*iCol) + j];
				adfvec[j] = aduvec[(j-1)] -  adUmat[(n*iCol) + j];
			}
			adfvec[iCol]  = *adw - adUmat[(n*iCol) + iCol];
			adhvec[iCol] = 0.0;
			for(iVal = 0 ; iVal < n; iVal++){
				buffveca[iVal] = adXmat[iCol + n*iVal];
			}
			dTemp = 1 +  cblasR_ddot(n, adfvec,incx, buffveca, incx);
			cblasR_dscal(n,(1/dTemp),adfvec,incx);
			row1 = n; col1 = n; row2= n; col2=1;
			cblasR_dgemv( CblasColMajor,CblasNoTrans, n, n, one, adXmat, n, adfvec ,incx, zero, buffveca, incx); // buffveca contains Xk*f
			
			for(j = 0 ; j < n ; j++){
				dTemp = adXmat[iCol + n*j];
				cblasR_daxpy( n , -dTemp, buffveca   ,incx,(adXmat + n*j) ,incx);
			}
			
			/*
			 Second Rank 1 Update 
			 Xk=Xk-((Xk*e*h'*Xk)/(1+h'*Xk*e);
			 */
			cblasR_dcopy( n , (adXmat + (n*iCol)) ,incx,buffveca,incx);
			dTemp = 1 +  cblasR_ddot(n, adhvec,incx, buffveca, incx);
			cblasR_dscal(n,(1/dTemp),adhvec,incx);
			cblasR_dcopy( n , (adXmat + n*iCol),incx,buffveca,incx);
			for(j = 0 ; j <n ; j++){
				dTemp =  cblasR_ddot(n, adhvec,incx, (adXmat + (n*j)), incx);
				cblasR_daxpy( n , -dTemp, buffveca  ,incx,(adXmat + n*j) ,incx);
			}
			symmetrize(adXmat,buffmat, n);
			cblasR_dcopy( n*n , buffmat,incx,adXmat,incx);
			
			/*
			 3. Update U 
			 Uk(i,i) = w;Uk(indi,i) = u;Uk(i,indi) = u';
			 */
			iTemp = 0;
			for(row = 0;row < n;row++){ 
				if ( row != iCol )
				{
					adUmat[row + n*iCol ] = aduvec[iTemp];
					adUmat[iCol + n*row ] = aduvec[iTemp];
					iTemp++;
				}
				else 
					adUmat[row + n*iCol ] = *adw;
			}
			adevec[iCol] = 0;
			i++ ; 
		}
		
		/*
		 4. Recompute Gap
		 DualityGap =  -n + trace(Sigma*Xk) + rho*(sum(sum(abs(Xk))));
		 */
		dCurrGap = - n + cblasR_ddot((n*n), adSigmat,incx, adXmat, incx) + dRho*(cblasR_dasum((n*n), adXmat,incx));
		dRelGap = dCurrGap/dOrigGap;
		
	}
	free(adTempVec);
	free(adevec); 
	free(adhvec);
	free(adfvec);
	free(adVinmat);
	free(aduo); 
	free(adwo); 
	free(aduvec);
	free(adw); 
	free(adbvec);
	free(buffveca); 
	free(buffvecb);
	free(buffmat); 	
	
	return dRelGap;
}

