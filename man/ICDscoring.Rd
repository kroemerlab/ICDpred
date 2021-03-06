% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ICDscoring.R
\name{ICDscoring}
\alias{ICDscoring}
\title{ICD inducing prediction}
\usage{
ICDscoring(CID, SDF)
}
\arguments{
\item{CID}{Pubchem CID#'}

\item{SDF}{Molecule 2D structure (must be sdf format)}
}
\value{
A list with computed scores and other useful information
}
\description{
Compute an ICD score
}
\examples{
test=ICDscoring(CID=4212) #Compute score for MTX

}
\author{
Allan Sauvat, \email{allan.sauvat@gustaveroussy.fr}
}
