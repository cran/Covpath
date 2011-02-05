\name{covpathCVForTS}
\alias{covpathCVForTS}
\title{COVPATH for Time Series with Cross-validation}
\description{
Generates solutions to the covariance selection problem (time series data) along a path of penalty values with cross-validation.
}
\usage{
covpathCVForTS(data, colFraction, t, rholist, numFolds, dates, targetCard, scalingChoice)
}
\arguments{
  \item{data}{ Input Data}
  \item{colFraction}{ The fraction of columns that are updated in each block coordinate run.}
  \item{t}{ Parameter to associated regularized problem.}
  \item{rholist}{ List of penalty values ( in descending order) to the covariance selection problem.}
  \item{numFolds}{ Number of folds for the cross-vaalidation}
  \item{dates}{ list of row indices in the input data that define the time series}
  \item{targetCard}{ Percentage of non-zero entries in the resulting conditional correlation matrix}
  \item{scalingChoice}{ Scaling options, set = 1 to have a unit norm covariance matrix and = 2 for correlation matrix}
}
\value{
a list containing:

\item{XSolTS}{ List of Estimated Inverse covariance Matrices}
\item{CmatTS}{ List of Estimated Conditional Correlation Matrices}
\item{ConfmatTS}{ List of link confidence matrices to the edges in the dependency graph}
}
\references{
A Pathwise Algorithm for Covariance Selection.
Optimization for Machine Learning, S. Sra, S. Nowozin and S. Wright, editors, MIT Press, 2011. 
}
\author{ Vijay Krishnamurthy, Alexandre d'Aspremont}
\keyword{ multivariate }