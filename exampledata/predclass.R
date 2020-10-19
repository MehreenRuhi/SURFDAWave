library("wavethresh")
library("glmnet")
library("matrixStats")
args = commandArgs(trailingOnly = TRUE)
infile <- args[1]
modelfile <- args[2]
nstat<-as.numeric(args[3])
cvfin<-readRDS(modelfile)
fam="DaubLeAsymm"
filt=as.numeric(8)
datatest=infile
#split into individual statistics
groups<-function(dat,nstat,min.scale){
	groupposition<-c(rep(1:nstat,(ncol(dat))/nstat))
	eachstat<-list()
	for(i in rep(1:nstat)){
		xid4<-which(groupposition==i)
		eachstat4<-dat[,c(xid4)]
		eachstat[[i]]<-wavetrandec(eachstat4,min.scale)	
	}
	groupstat<-do.call("cbind",eachstat)
	return(groupstat)}


#conduct discrete wavelet transform
#adapted from (Zhao Y, Ogden RT, Reiss PT. Wavelet-based LASSO in functional linear regression. J Comput Graph Stat. 2012;21(3):600–617. doi:10.1080/10618600.2012.679241))
wavetrandec<-function(eachstat=eachstat[[i]],min.scale=min.scale){ywdS<-apply(eachstat,1,wd,filter.number=8,family="DaubLeAsymm")
	wdt<-matrix(0, nrow(eachstat),ncol(eachstat))
	
		for(i in 1:nrow(eachstat)){
		base<-0
		for(j in ((log2(ncol(eachstat))-1):min.scale)){
			temp<-accessD(ywdS[[i]], level=j, boundary=F)
			wdt[i,((base+1):(base+2^j))]<-temp
			base<-base+2^j
		}
		wdt[i,((ncol(eachstat)-2^min.scale+1):ncol(eachstat))]<-accessC(ywdS[[i]], level=min.scale, boundary=F)
	}
return(wdt)
}

temptest<-as.matrix(read.csv(datatest,header=FALSE,sep=','))
train_info<-strsplit(modelfile,'.',fixed=TRUE)[[1]]
traindat<-as.matrix(read.csv(train_info[[1]],header=FALSE,sep=','))
centest<-sweep(temptest,2,colMeans(traindat))
normtest<-sweep(centest,2,colSds(traindat),"/")
wavtot<-groups(normtest,nstat,as.numeric(train_info[[2]]))
fit<-as.numeric(predict(cvfin,newx=wavtot,s="lambda.1se",type='class'))
fit2<-predict(cvfin,newx=wavtot,s="lambda.1se",type='response')
fit3<-as.data.frame(fit2)
write.table(cbind(fit,fit3), file=paste(datatest,'.predclass',sep=''),sep=',',col.names=FALSE)

