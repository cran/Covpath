\name{covpath}
\alias{covpath}
\title{COVPATH}
\description{
Generates solutions to the covariance selection problem along a path of penalty values.
}
\usage{
covpath(Sigma, rholist, t, colFraction, rhomax)
}
\arguments{
  \item{Sigma}{ Sample symmetric covariance matrix.}
  \item{rholist}{ List of penalty values ( in descending order) to the covariance selection problem.}
  \item{t}{ Parameter to associated regularized problem.}
  \item{colFraction}{ The fraction of columns that are updated in each block coordinate run.}
  \item{rhomax}{ Maximum penalty value, beyond which the solution has an explicit diagonal solution}
}
\value{
a list containing:

\item{Covmat}{ Array of estimated covariance}
\item{invCovmat}{ Array of estimated inverse covariance matrices}
\item{relativeGaps}{ Ratio of duality gap of solution to the duality gap of the solution to the unpenalized problem}
}
\references{
A Pathwise Algorithm for Covariance Selection.
Optimization for Machine Learning, S. Sra, S. Nowozin and S. Wright, editors, MIT Press, 2011. 
}
\author{ Vijay Krishnamurthy, Alexandre d'Aspremont}
\keyword{ multivariate }
