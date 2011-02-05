


#Load FX Data
data = read.csv("../data/FxData.csv", header =TRUE)
fulldata = data

 
# scalingChoice 2  covm=corrcoef(fulldata);
scalingChoice =2;
 
#Path Algorithm 
nrhos = 100
rhomax = 1
logrhomax = log10(rhomax)
if (scalingChoice == 1) logrhomin=logrhomax-5 else logrhomin=logrhomax-0.2
rholist = 10^(seq(from = logrhomax, to = logrhomin, length.out= (nrhos+1)))

 
# Covpath Parameters
colFraction =1
t = 10^-7
 
# CrossValidation Parameter
targetCard = 10
plotmaxConf = 0
numFolds = 10

dates = c(400,500,600)
#dates = numeric(0)
tic = proc.time()
cvTSsol = covpathCVForTS(data,colFraction,t,rholist,numFolds,dates,targetCard,scalingChoice)
toc = (proc.time() -tic)[3]
print(toc)

condcorrTS = cvTSsol$CmatTS
numsamples = nrow(data)
p = ncol(data)


#Plot Using Protovis

ColNames =   c("EUR" ,"JAP","GBP","CHF","AUD","NZD","CAD","NOK","SEK","PLN","CZK","HUF","TRY","ZAR","MXN","BRL","KRW","INR","IDR","SGD")
GroupLabels = c(1,2,1,1,2,2,4,1,1,1,1,1,3,3,3,4,4,2,2,2) 
jsFileName = "fxTS"
htmlFileName = "fxdataTS"

plotDependencyTimeSeries(condcorrTS,targetCard,dates,numsamples,p,ColNames,GroupLabels,jsFileName,htmlFileName)

