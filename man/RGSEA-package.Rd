\name{RGSEA-package}
\alias{RGSEA-package}
\alias{RGSEA}
\docType{package}
\title{
Random Gene Set Enrichment Analysis (RGSEA)
}
\description{
This is the package for similarity identifucation and classification of transcriptome data
}
\details{
\tabular{ll}{
Package: \tab RGSEA\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-04-22\cr
License: \tab GPL(>=3)\cr
}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Chengcheng Ma
Maintainer: tcmcc1986@gmail.com
}
\references{
~~ Literature or other references for background information ~~
}

\keyword{ package }
\seealso{
Song L, Langfelder P, Horvath S. Random generalized linear model: a highly 
accurate and interpretable ensemble predictor[J]. BMC bioinformatics, 2013, 14(1): 5.
Subramanian A, Tamayo P, Mootha V K, et al. Gene set enrichment analysis: 
a knowledge-based approach for interpreting genome-wide expression profiles[J].
 Proceedings of the National Academy of Sciences of the United States of America,
  2005, 102(43): 15545-15550.
}
\examples{
if(interactive()) {
    data(e1)
    data(e2)
    RGSEAfix(e1,e2, queryclasses=colnames(e1), refclasses=colnames(e2),      
random=20000, featurenum=1000, iteration=100)->test
}
}
