###Prediction
library(glmnet)
library(wavethresh)
library(matrixStats)
args = commandArgs(trailingOnly = TRUE)
infile <- args[1]
data1=infile
nstat<-as.numeric(args[2])
#number of windows
nwin<-128
#number of features
nfeat<-nwin*nstat

#Read data
tempdat<-as.matrix(read.csv(data1,header=FALSE,sep=','))
set.seed(1)
#Split data into stats
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

normFunc <- function(trdat){(trdat-mean(trdat))/sd(trdat)}
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
#elastic net values for cross validation
alphavals<-seq(0.9,1,0.1)
#scales for cross validation
minvals<-seq(0,4)

resmat<-matrix(0,nrow=length(alphavals)*length(minvals),ncol=3)

#shuffle data
shuftempdat<-tempdat[sample(nrow(tempdat)),]
#10-fold cross validation
cvfolds<-cut(seq(1,nrow(shuftempdat)),breaks=10,labels=FALSE)
lind<-0
for(minscale in minvals){
	for(alph in alphavals){
	lind<-lind+1
	countlev<-0
	for(f in seq(1,10)){
		testinds<-which(cvfolds==f,arr.ind=TRUE)
		traindat<-shuftempdat[-testinds,1:nfeat]
		testdat<-shuftempdat[testinds,1:nfeat]
		
		normtrain<-apply(traindat,2,normFunc)
		centest<-sweep(testdat,2,colMeans(traindat))
		normtest<-sweep(centest,2,colSds(traindat),"/")
		grouptest<-groups(normtest,nstat,minscale)
		grouptrain<-groups(normtrain,nstat,minscale)
		cv<-cv.glmnet(x=grouptrain,y=shuftempdat[-testinds,(nfeat+1):(ncol(shuftempdat))],alpha=alph,standardize=FALSE,family="mgaussian")
		preds<-predict(cv, newx = grouptest, s = "lambda.1se")
		preddf<-as.data.frame(preds)
		act<-shuftempdat[testinds,(nfeat+1):(ncol(shuftempdat))]
		difs<-preddf-act
		countlev<-countlev+sum(abs(difs))
		
	}
	resmat[lind,1]<-alph
	resmat[lind,2]<-minscale
	resmat[lind,3]<-countlev

	}}
#best model
finmin<-which.min(resmat[,3])
#final training
traindatfin<-apply(shuftempdat[,1:nfeat],2,normFunc)
traindatfin2<-groups(traindatfin,nstat,resmat[finmin,2])
finalmodel<-cv.glmnet(x=traindatfin2,y=shuftempdat[,(nfeat+1):(ncol(shuftempdat))],alpha=resmat[finmin,1],standardize=FALSE,family="mgaussian")
saveRDS(finalmodel,file=(paste(infile,'.',resmat[finmin,2],'.rds',sep='')))
