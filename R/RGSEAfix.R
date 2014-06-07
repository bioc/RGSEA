RGSEAfix <-
    function(query, reference, queryclasses, refclasses,random = 5000, 
    featurenum = 500, iteration = 100)
{
    query<-as.matrix(query)
    if(dim(query)[2] != length(queryclasses))
	stop("Query classes dose not equal to number of query data.")
    if(dim(reference)[2] != length(refclasses))
	stop("Reference classes dose not equal to number of reference data.")

    internames <- intersect(rownames(query), rownames(reference))
    e1 <- as.matrix(query[internames, ])
    e2 <- reference[internames, ]
    e1[which(is.na(e1))]<-0
    e2[which(is.na(e2))]<-0

    RGSEAresult<-matrix(0,dim(query)[2],dim(reference)[2])
    colnames(RGSEAresult)<-colnames(e2)
    rownames(RGSEAresult)<-colnames(e1)
    featuretimes<-numeric(dim(e1)[1])
    names(featuretimes)<-rownames(e1)

    if(random > dim(e1)[1])
	stop("Number of sample data exceed number of total features.")
    if(featurenum*2 > random)
	stop("Number of selected features exceed number of randomly sampled features.")
    if(!is.character(queryclasses))
	stop("Mode of queryclasses is not character.")
    if(!is.character(refclasses))
	stop("Mode of refclasses is not character.")

    for (a in 1:dim(e1)[2])
    {
	winnercollector <- numeric()
	cat(sprintf("Processing #%d...\n", a))
        for (iter in 1:iteration)
	{
	    escollector <- numeric()
	    sampleindex <- sample(1:dim(e1)[1], random)
	    s1 <- sort(e1[sampleindex,a], decreasing = TRUE)
	    up1 <- names(s1[1:featurenum])
	    down1 <- names(s1[(random-featurenum + 1):random])
            for (b in 1:dim(e2)[2])
	    {
		s2<-sort(e2[sampleindex,b], decreasing = TRUE)

		um1 <- sort(match(up1,names(s2)))
	        dm1 <- sort(match(down1,names(s2)))
		uscore1 <- (1:featurenum)/featurenum - (um1 - 1:featurenum)/
                               (random - featurenum)
		dscore1 <- (1:featurenum)/featurenum - (dm1 - 1:featurenum)/
                               (random - featurenum)
                ues1 <- max(uscore1)
                des1 <- min(dscore1)
                es1 <- 1 - (ues1 - des1)/2

                up2 <- names(s2[1:featurenum])
		down2 <- names(s2[(random - featurenum + 1):random])
		um2 <- sort(match(up2,names(s1)))
		dm2 <- sort(match(down2,names(s1)))
	        uscore2 <- (1:featurenum)/featurenum - (um2 - 1:featurenum)/(random - featurenum)
                dscore2 <- (1:featurenum)/featurenum - (dm2 - 1:featurenum)/(random - featurenum)
                ues2 <- max(uscore2)
                des2 <- min(dscore2)
                es2 <- 1 - (ues2 - des2)/2

                es <- (es1 + es2)/2
                escollector <- c(escollector, es)
            }
	    winner <- which.min(escollector)
	    winnercollector <- c(winnercollector, winner)
	    if(queryclasses[a] == refclasses[winner])
	    {
		tmp <- table(c(up1,down1))
		featuretimes[names(tmp)] <- featuretimes[names(tmp)] + tmp
	    }
	}
	winners <- table(winnercollector)
	RGSEAresult[a, as.numeric(names(winners))] <- table(winnercollector)
    }
    list(RGSEAresult, featuretimes)
}
