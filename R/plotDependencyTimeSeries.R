
plotDependencyTimeSeries <-function(CmatTS,targetCard,dates,numsamples,p,colNames,groupLabels,jsFileName,htmlFileName){

condcorrTS = CmatTS
allcolNames = colNames
allgroupLabels = groupLabels
thrsQuantile = targetCard
javaScriptFileName = jsFileName


jsfile = paste(javaScriptFileName,".js",sep ="")
htmlfile = paste(htmlFileName,".html",sep ="")
cat(file = jsfile, append = "FALSE", "\n")

l = length(dates)
if ( l < 2) numPlots = 1 else numPlots = l 

for (i in 1:numPlots) {
	
	st = (i-1)*p*p+1
  en = st + (p*p) -1 
  C = matrix(condcorrTS[st:en],p,p)
  
  sol  = thresholdCondCorrMatrix(C,targetCard)
X3 = sol$prunedC
ikeep = sol$ikeep
numEdges = 0

  colNames = allcolNames[sol$ikeep]
groupLabels = allgroupLabels[sol$ikeep]
tempstr = paste("data",i, sep = "")
	cat(file = jsfile, append = "TRUE", "var",tempstr, "= {","\n",sep = " ")
	cat(file = jsfile, append = "TRUE", "nodes:[", "\n")
	
n = nrow(X3)
matConf = 3*matrix(1,n,n)	
	for ( j in 1:n){
		str = colNames[j]
		label = groupLabels[j]
		cat(file = jsfile, append = "TRUE", "{nodeName:","\"",str,"\"",", group:",label,"},", "\n",sep = "")
			}
				cat(file = jsfile, append = "TRUE", "],", "\n")
			cat(file = jsfile, append = "TRUE", "links:[", "\n")
	for (src in 1:n){
		for (tar in (1:(src-1))){
			if ( src != 1 ) {{}
		if  ( (X3[src,tar]) != 0 ) {
			numEdges = numEdges +  1
			cat(file = jsfile, append = "TRUE", "{source:",(src-1),", target:",(tar-1), ", value:", 1,", linkconf:",matConf[src,tar],", condcorr:",X3[src,tar],"},",  "\n",sep = "")
			}
			}
			}
			}
						cat(file = jsfile, append = "TRUE", "]", "\n")
					cat(file = jsfile, append = "TRUE", "};", "\n")
}

numNodes =n
graphWidth = numNodes*50
graphHeight = numNodes*50

createProtovishtmlforCV(numPlots, numsamples,dates,htmlfile,jsfile,graphWidth,graphHeight)
system(paste("open",htmlfile))
return(numEdges)

}

  

 
	
	
	
	
	
	
	
	


