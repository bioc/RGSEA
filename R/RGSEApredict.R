RGSEApredict <-
    function(RGSEAresult, refclasses)
{
    RGSEAclass <- matrix(,dim(RGSEAresult)[1], length(unique(refclasses)))
    colnames(RGSEAclass) <- levels(factor(refclasses))
    rownames(RGSEAclass) <- rownames(RGSEAresult)
    for(i in 1:dim(RGSEAresult)[1])
    {
	tmp <- tapply(RGSEAresult[i, ], refclasses, mean)
	RGSEAclass[i,] <- tmp/sum(tmp)
    }
    RGSEAclass
}
