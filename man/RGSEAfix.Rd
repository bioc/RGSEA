\name{RGSEAfix}
\alias{RGSEAfix}
\title{
Random Gene Set Enrichment Analysis with fixed number of features
}
\description{
This is the function for classification and feature selection with fixed number 
of features from top and bottom of the subtset features.  
}
\usage{
RGSEAfix(query, reference, queryclasses, refclasses, random = 5000, featurenum
 = 500, iteration = 100)
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
the classes of query data, just give it a character vector equal to the number
 of query data.
}
  \item{refclasses}{
A character vector. It contains the classes of reference data. You must know it.
}
  \item{random}{
A numeric variable. The number of features in the subset randomly sampled from 
the whole features each time.
}
  \item{featurenum}{
A numeric varialbe. The number of features selected from top and bottom of the 
subset respectivelly.
}
  \item{iteration}{
A numeric varialbe. The times of random sampling.
}
}
\value{
[1] The times of each sample in the reference dataset is the most similar to 
the query data.
[2] The frequencey of features selected from the top and bottom of the
 subsets from the query data, if the query data is correcly classified.
}
\author{
Chengcheng Ma
}
\examples{
if(interactive()) {
    data(e1)
    data(e2)
    RGSEAfix(e1,e2, queryclasses=colnames(e1), refclasses=colnames(e2),      
random=20000, featurenum=1000, iteration=100)->test
}
}