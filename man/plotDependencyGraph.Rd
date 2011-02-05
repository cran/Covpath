\name{plotDependencyGraph}
\alias{plotDependencyGraph}
\title{Visualizing Dependency Structure}
\description{
Generates a force directed layout of the conditional correlation graph using Protovis}
\usage{
plotDependencyGraph(Cmat, colNames, groupLables, targetCard, jsFileName, htmlFileName)
}
\arguments{
  \item{Cmat}{ Conditional Correlation Matrix for visualization} 
  \item{colNames}{ Character vector node names for Protovis iinput}
  \item{groupLables}{ Characted vector of group lables for Protovis Input}
  \item{targetCard}{ Percentage of edges in the resulting conditional correlation graph}
  \item{jsFileName}{ String specifying name of Java Scrpit file for Protovis Data}
  \item{htmlFileName}{ String specifying name of HTML file for Protovis Data}
}
\value{
\item{numEdges}{ Number of edges in the resulting dependency graph}
}
\references{
A Pathwise Algorithm for Covariance Selection.
Optimization for Machine Learning, S. Sra, S. Nowozin and S. Wright, editors, MIT Press, 2011. 
}
\author{ Vijay Krishnamurthy, Alexandre d'Aspremont}
\keyword{ multivariate }