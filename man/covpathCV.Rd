\name{covpathCV}
\alias{covpathCV}
\title{COVPATH with Cross-validation}
\description{
Generates solutions to the covariance selection problem along a path of penalty values with cross-validation.
}
\usage{
covpathCV(data, rholist, numFolds, targetCard, scalingChoice, t, colFraction)
}
\arguments{
  \item{data}{ Input data}
  \item{rholist}{ List of penalty values ( in descending order) to the covariance selection problem.}
  \item{numFolds}{ Number of folds for the cross-vaalidation}
  \item{targetCard}{ Percentage of non-zero entries in the resulting conditional correlation matrix}
  \item{scalingChoice}{ Scaling options, set = 1 to have a unit norm covariance matrix and = 2 for correlation matrix}
  \item{t}{ Parameter to associated regularized problem.}
  \item{colFraction}{ The fraction of columns that are updated in each block coordinate run.}
}
\value{
a list containing:

\item{Xmat}{ Estimated Inverse covariance Matrix}
\item{Cmat}{ Estimated Conditional Correlation Matrix}
\item{Confmat}{ Matrix of link confidences to the edges in the dependency graph}
}
\references{
A Pathwise Algorithm for Covariance Selection.
Optimization for Machine Learning, S. Sra, S. Nowozin and S. Wright, editors, MIT Press, 2011. 
}
\author{ Vijay Krishnamurthy, Alexandre d'Aspremont}
\keyword{ multivariate }
