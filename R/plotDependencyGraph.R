
plotDependencyGraph <-function(Cmat,colNames,groupLables,targetCard,jsFileName,htmlFileName){

C =Cmat
allcolNames = colNames
allgroupLabels = groupLables
thrsQuantile = targetCard
javaScriptFileName = jsFileName
jsfile = paste(javaScriptFileName,".js",sep ="")
htmlfile = paste(htmlFileName,".html",sep ="")
cat(file = jsfile, append = "FALSE", "\n")

sol  = thresholdCondCorrMatrix(C,thrsQuantile)
X3 = sol$prunedC
ikeep = sol$ikeep

numEdges = 0
colNames = allcolNames[sol$ikeep]
groupLabels = allgroupLabels[sol$ikeep]
i =1
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

numNodes =n
graphWidth = numNodes*50
graphHeight = numNodes*50
	
createProtovishtml(htmlfile,jsfile,graphWidth,graphHeight)
system(paste("open",htmlfile))

return(numEdges)



}