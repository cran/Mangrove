\name{Mangrove-package}
\alias{Mangrove-package}
\alias{Mangrove}
\docType{package}
\title{
Mangrove: Risk prediction on trees
}
\description{
\code{Mangrove} is an R package for performing genetic risk
prediction from genotype data.  You can use it to perform risk prediction for
individuals, or for families with missing data.
}
\details{
\tabular{ll}{
Package: \tab Mangrove\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-10-24\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}

For a detailed example of how to use \code{Mangrove}, see \code{vignette("Mangrove")}.
}
\author{
Luke Jostins <lj4@sanger.ac.uk>
}
\keyword{ package }
\examples{
data(exampleORs)
data(famped)
tree <- initialiseTree()
tree$addPed(famped,exampleORs)
sam <- tree$getPrevs(exampleORs,K=0.02)
plotNaivePrev(famped,0.02)
plot(sam)

}
