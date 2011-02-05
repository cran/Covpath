

covpath <- function(Sigma,rholist,t,colFraction,rhomax){

n = nrow(Sigma)

inputRhoList = rholist[2:length(rholist)]
outputRhoList = rholist[2:length(rholist)]
m = length(outputRhoList);
tol = 10^-5
Xorig = solve(Sigma);
mu = sum(sum(abs(Xorig)))
OrigGaps = mu*inputRhoList
rhovec = c(rhomax,inputRhoList)
l = length(rhovec)

XSol = c()
USol = c()
relGaps =c()
U0 = diag(diag(Sigma)) + 0.999*rhomax*diag(n)
Uk = U0
Xk = solve(Uk)

for (i in 2 : l)
  {
    
    rhocurr = rhovec[i-1]
    rhonew = rhovec[i]
    h = rhonew - rhocurr
    scale = (rhonew/rhocurr)
    Uk = (1 - scale)*Sigma + scale*Uk   
    Xk = solve(Uk)
    dOrigGap = mu*rhonew;
    update = .Call("bcdCorrector",Sigma,Uk,rhonew,t,tol,colFraction,Xk,dOrigGap, PACKAGE = "Covpath")
    Uk = update$Uk
    Xk = update$Xk
    XSol = c(XSol,Xk)
    USol = c(USol,Uk)
    relGaps = c(relGaps, update$relGap)
     }

  pathSoln = list(Covmat = USol, invCovmat = XSol, relativeGaps = relGaps)
  return(pathSoln)
  
}
