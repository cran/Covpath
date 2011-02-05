
#Load Rates Data
rts = read.csv("../data/Forwards.csv", header = FALSE)
Sigma = cov(rts)
Sigma = Sigma/(norm(Sigma,"F"))
Sigma = as.matrix(Sigma)

#Covpath Paraameters
colFraction =1
t = 10^-7

#Path Algorithm 
nrhos = 10
rhomax = max(diag(Sigma))
logrhomax = log10(rhomax)
logrhomin=logrhomax-1.5
rholist = 10^(seq(from = logrhomax, to = logrhomin, length.out= (nrhos+1)))
covpathsol = covpath(Sigma,rholist,t,colFraction,rhomax)
XSol = covpathsol$invCovmat
RelGap = covpathsol$relativeGaps
cards = thresholdAndPrintResults(Sigma,rholist,XSol,RelGap)

#Pick a index/solution correponding to a penalty and Generate Conditional Correlation Matrix 

ind = 10 
n = nrow(Sigma)
st = (ind-1)*n*n +1
en = st + (n*n) -1 
X = matrix(XSol[st:en], nrow =n)
C = makeCondMat(X)

#Plot Using Protovis
targetCard = 10
ColNames =paste(1:n)
GroupLabels = rep(1,n)
jsFileName = "rates"
htmlFileName = "FwdRates"
numEdges = plotDependencyGraph(C,ColNames,GroupLabels,targetCard,jsFileName,htmlFileName)
