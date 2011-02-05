\name{plotDependencyTimeSeries}
\alias{plotDependencyTimeSeries}
\title{Visualizing Dependency Structure across time}
\description{
Generates a panel of conditional correlation graph susing Protovis}
\usage{
plotDependencyTimeSeries(CmatTS, targetCard, dates, numsamples, p, colNames, groupLabels, jsFileName, htmlFileName)
}
\arguments{
  \item{CmatTS}{ List of Conditional Correlation Matrix for visualization}
  \item{targetCard}{ Percentage of edges in the resulting conditional correlation graph}
  \item{dates}{ list of row indices in the input data that define the time series}
  \item{numsamples}{ Number of samples in the data = number of rows of input data matrix}
  \item{p}{ Dimension of data = number of columns of input data matrix }
  \item{colNames}{ Character vector node names for Protovis iinput}
  \item{groupLabels}{ Characted vector of group lables for Protovis Input}
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