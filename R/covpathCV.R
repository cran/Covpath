
covpathCV <- function(data,rholist,numFolds,targetCard,scalingChoice,t,colFraction){

	
  thrsQuantile = 100 - targetCard
                                        #Covpath Parameters
 
  rhomax = rholist[1]
  inputRhoList = rholist[-1]
  outputRhoList =  rholist[-1]
  l = length(outputRhoList);
  freqMatrix = c()
  networkConfMartix = c()
  sumLinkConf = c()
  ns = nrow(data)
  p = ncol(data)
  K = numFolds;
  kappa =  floor(ns/K)
  rp = sample.int(ns,ns)

  if (scalingChoice == 1) {
    fullcovm = cov(data)
    fullcovm = fullcovm/(norm(fullcovm,"F"))
  }
  else fullcovm = cor(data)
  n = ncol(fullcovm)
  freqMatrix = rep(0,l*n*n)
  networkConfMatrix = rep(0,l*n*n) 
  for (k in 1:K) {
  	if (k != K)
        trainidx = rp[-(((k-1)*kappa + 1):(k*kappa))]
    else
        trainidx = rp[-(((k-1)*kappa + 1):ns)]
    trainData = data[trainidx,1:p]
    if (scalingChoice == 1) {
      covm = cov(trainData)
      covm = covm/(norm(covm,"F"))
    }
    else covm = cor(trainData)
    sol = covpath(covm,rholist,t,colFraction,1)
    XSol = sol$invCovmat
    for (i in 1:l){
      st = (i-1)*n*n +1
      en = st + (n*n) -1 
      Umat = matrix(XSol[st:en], nrow =n)
     X = makeCondMat(Umat) 
      pnnz = (targetCard)/100
	  thrc=as.vector(quantile(c(abs(X)),1-pnnz))
      strind = (i-1)*n*n + 1
      endind = strind + n*n -1
      freqMatrix[strind:endind] = freqMatrix[strind:endind] +  c(abs(X) > thrc )

    }
  }
  freqMatrix = freqMatrix/K

                                        #Entire Data
  print(" COVPATH Performance Dimension of input matrix")
  print(n)
  fullsol = covpath(fullcovm,rholist,t,colFraction,rhomax)
  XSol = fullsol$invCovmat
  RelGap =fullsol$relativeGaps
  cards = c()
  sumLinkConf = rep(0,l)
  
  for (i in 1:l){
    
    st = (i-1)*n*n +1
    en = st + (n*n) -1 
    Xmat = matrix(XSol[st:en], nrow =n)
        
    #F = matrix(freqMatrix[st:en], nrow =n)
     F = freqMatrix[st:en]
	xx = hist(c(abs(Xmat)),breaks = 98 , plot = 0)

    thrshld = xx$mids[2]
    cardk = 100*(sum(abs(Xmat)>thrshld))/(n*n)
    cards = c(cards,cardk)
    X = makeCondMat(Xmat)
        
    pnnz = (targetCard)/100
    thrc=as.vector(quantile(c(abs(X)),1-pnnz))
	binmat=2*(c(abs(X)>thrc))*F
    binmat = binmat + (c(abs(X)<=thrc)-F)
    networkConfMartix[st:en] = c(binmat)
    sumLinkConf[i] = sum((binmat))/(n*n)
    
  }
  
  S.No = 1:l
rho = outputRhoList
RelativeImp = RelGap
CardPerc = cards
resc = cbind(S.No,rho,RelativeImp,CardPerc)
print(resc)



    scard = sort(cards,index.return =1)
    i = which(scard$x>targetCard)
    i = i[1]
    i = scard$ix[i]
  
  st = (i-1)*n*n +1
  en = st + (n*n) -1 
  X = XSol[st:en]
  Xcard = cards[i];
  maxConf = sumLinkConf[i]
  netConfMatrix = networkConfMartix[st:en]
  CondMat= c(makeCondMat(matrix(X,n,n)))
 pathSolnCV = list(Xmat = X,Cmat = CondMat,Confmat = netConfMatrix)
  return(pathSolnCV)

  
}




 










