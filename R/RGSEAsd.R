RGSEAsd <-
    function(query, reference, queryclasses, refclasses,random = 5000, sd = 2,
    iteration = 100)
{
    query <- as.matrix(query)
    if(dim(query)[2] != length(queryclasses))
        stop("Query classes dose not equal to number of query data.")
    if(dim(reference)[2] != length(refclasses))
	stop("Reference classes dose not equal to number of reference data.")

    internames <- intersect(rownames(query), rownames(reference))
    e1 <- as.matrix(query[internames,])
    e2 <- as.matrix(reference[internames,])
    e1[which(is.na(e1))] <- 0
    e2[which(is.na(e2))] <- 0

    RGSEAresult <- matrix(0,dim(query)[2],dim(reference)[2])
    colnames(RGSEAresult) <- colnames(e2)
    rownames(RGSEAresult) <- colnames(e1)
    featuretimes <- numeric(dim(e1)[1])
    names(featuretimes) <- rownames(e1)

    if(random > dim(e1)[1])
	stop("Number of sample data exceed number of total features.")
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
	    sampleindex <- sample(1:dim(e1)[1],random)
	    s1 <- sort(e1[sampleindex,a], decreasing = TRUE)
	    m1 <- mean(s1)
	    sd1 <- sd(s1)
	    up1 <- names(s1[which(s1 > m1 + sd*sd1)])
	    down1 <- names(s1[which(s1 < m1 - sd*sd1)])
	    lup1 <- length(up1)
	    ldown1 <- length(down1)
            for (b in 1:dim(e2)[2])
	    {
	        s2 <- sort(e2[sampleindex,b], decreasing = TRUE)
	 	m2 <- mean(s2)
	        sd2 <- sd(s2)
		up2 <- names(s2[which(s2 > m2 + sd*sd2)])
		down2 <- names(s2[which(s2 < m2 - sd*sd2)])
		lup2 <- length(up2)
		ldown2 <- length(down2)

		um1 <- sort(match(up1,names(s2)))
		dm1 <- sort(match(down1,names(s2)))
		uscore1 <- (1:lup1)/lup1-(um1 - 1:lup1)/(random - lup1)
		dscore1 <- (1:ldown1)/ldown1-(dm1 - 1:ldown1)/(random - ldown1)
		ues1 <- min(uscore1)
		des1 <- min(dscore1)
		es1 <- 1- (ues1 - des1)/2

		um2 <- sort(match(up2,names(s1)))
		dm2 <- sort(match(down2,names(s1)))
		uscore2 <- (1:lup2)/lup2-(um2 - 1:lup2)/(random - lup2)
		dscore2 <- (1:ldown2)/ldown2-(dm2 - 1:ldown2)/(random - ldown2)
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
	RGSEAresult[a, as.numeric(names(winners))] <- winners
    }
    list(RGSEAresult, featuretimes)
}
