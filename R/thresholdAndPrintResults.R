
thresholdAndPrintResults <- function(Sigma,rholist,XSol,RelGap){
	
n = ncol(Sigma)
l = length(rholist)
outputRhoList = rholist[2:l]

cards=c()
for (i in 1:(l-1)){

  st = (i-1)*n*n +1
  en = st + (n*n) -1 
  Xmat = matrix(XSol[st:en], nrow =n)
  xx = hist(c(abs(Xmat)),breaks = 98 , plot = 0)
  thrshld = xx$mids[2]
  cardk = 100*(sum(abs(Xmat)>thrshld))/(n*n)
  cards = c(cards,cardk)
}
print("Covpath Results ")
S.No = 1:(l-1)
rho = outputRhoList
CardPerc = cards
resc = cbind(S.No,rho,RelGap,CardPerc)
print(resc)
return(cards)
}