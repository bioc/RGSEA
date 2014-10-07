test_RGSEAsd<-function(){
    set.seed(1)
    
	data(cmap)
	checkTrue(which.max(RGSEAsd(cmap[,1],cmap[,2:6], 
         queryclasses=colnames(cmap)[1], 
         refclasses=colnames(cmap)[2:6], 
         random=5000, sd=2, iteration=100)[[1]])==3)
}

