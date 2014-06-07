\name{RGSEAsd}
\alias{RGSEAsd}
\title{
Random Gene Set Enrichment Analysis features selected based on standard 
deviation from the mean value
}
\description{
This is the function for classification with features from top and bottom 
deviating from the mean value of the whole transcriptome for a certain standard
 deviations of the subtset features.  
}
\usage{
RGSEAsd(query, reference, queryclasses, refclasses, random = 5000, sd = 2, iteration = 100)
}
\arguments{
  \item{query}{
A matrix, The query data. This is the data which the research wants to know 
the class.
}
  \item{reference}{
A matrix. The reference data. Based of the reference data, the research infer 
the class of query data. 
}
  \item{queryclasses}{
A character vector. It contains the classes of query data. If you don't know
 the classes of query data, just give it a character vector equal to the number of query data.
}
  \item{refclasses}{
A character vector. It contains the classes of reference data. You must know it.
}
  \item{random}{
A numeric variable. The number of features in the subset randomly sampled from 
the whole features each time.
}
  \item{sd}{
number of standard deviations the features selected from the subset deviate from 
the mean value of the subset.
}
  \item{iteration}{
A numeric varialbe. The times of random sampling.
}
}
\value{
[1] The times of each sample in the reference dataset is the most similar to
 the query data.
[2] The frequencey of features selected from the top and bottom of the subsets 
from the query data, if the query data is correcly classified.
}
\author{
Chengcheng Ma
}
\examples{
if(interactive()) {
    data(cmap)
    test <- RGSEAsd(cmap[,1],cmap[,2:6], queryclasses=colnames(cmap)[1], 
      refclasses=colnames(cmap)[2:6], random=5000, sd=2, iteration=100)
}
}
