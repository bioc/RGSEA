test_RGSEAfix<-function(){
	data(e1)
	data(e2)
	checkTrue(RGSEAfix(e1,e2, queryclasses=colnames(e1), refclasses=colnames(e2), random=20000, featurenum=1000, iteration=100)[[1]][1,1]>50)
}
