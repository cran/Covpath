

#Load Rates Data
rts = read.csv("../data/Forwards.csv", header = FALSE)
data = as.matrix(rts)
nrts = ncol(data)
fulldata = data

 
# scalingChoice 1  covm = cov(fulldata);covm=covm/norm(covm);
# scalingChoice 2  covm=corrcoef(fulldata);
scalingChoice =2;
 
#Path Algorithm 
nrhos = 50
rhomax = 1
logrhomax = log10(rhomax)
if (scalingChoice == 1) logrhomin=logrhomax-1.3 else logrhomin=logrhomax-0.1
rholist = 10^(seq(from = logrhomax, to = logrhomin, length.out= (nrhos+1)))

 
# Covpath Parameters
colFraction =1
t = 10^-7
 
# CrossValidation Parameter
targetCard = 10
plotmaxConf = 0
numFolds = 10

dates = c(1,200, 300)
# dates = numeric(0)

cvTSsol = covpathCVForTS(data,colFraction,t,rholist,numFolds,dates,targetCard,scalingChoice)

condcorrTS = cvTSsol$CmatTS
numsamples = nrow(data)
p = ncol(data)


#Plot Using Protovis
ColNames =paste(1:nrts)
GroupLabels = rep(1,nrts)
jsFileName = "ratesTS"
htmlFileName = "FwdRatesTS"

plotDependencyTimeSeries(condcorrTS,targetCard,dates,numsamples,p,ColNames,GroupLabels,jsFileName,htmlFileName)
	
