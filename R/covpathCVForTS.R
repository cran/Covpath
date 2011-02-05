
covpathCVForTS <- function(data,colFraction,t, rholist,numFolds,dates,targetCard,scalingChoice){
	
	l = length(dates)
	ns = nrow(data)
	p = ncol(data)
	XSolTS = c()
	condcorrTS = c()
	networkConfTS =c()
	
	if (l > 1){
		
		for (i in 1:(l-1)) {
			startdt = dates[i]
			enddt = dates[i+1] -1
			datap = data[startdt:enddt,1:p] 
			cvsol = covpathCV(datap,rholist,numFolds,targetCard,scalingChoice,t,colFraction)
			XSolTS = c(XSolTS, cvsol$Xmat)
			condcorrTS = c(condcorrTS, cvsol$Cmat)
			networkConfTS =c(networkConfTS,cvsol$Confmat)
					}
		
		}
			cvsol = covpathCV(data,rholist,numFolds,targetCard,scalingChoice,t,colFraction)
			XSolTS = c(XSolTS, cvsol$Xmat)
			condcorrTS = c(condcorrTS, cvsol$Cmat)
			networkConfTS =c(networkConfTS,cvsol$Confmat)


		TSSoln = list(XSolTS = XSolTS, CmatTS = condcorrTS, ConfmatTS = networkConfTS)
 	 return(TSSoln)
		
		}
		
	
	
	 

	