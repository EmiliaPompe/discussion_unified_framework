setwd("~/discussion_unified_framework/couplingdeduplication/inst/")
##################DATA PREPARATION############################
library(RecordLinkage)
data(RLdata500)
myRLDATA=cbind(id=identity.RLdata500,RLdata500)


M1=matrix(unlist(strsplit(soundex(myRLDATA$fname_c1),split="")),ncol=4,byrow=TRUE)
M2=matrix(unlist(strsplit(soundex(myRLDATA$lname_c1),split="")),ncol=4,byrow=TRUE)
M3=matrix(unlist(strsplit(as.character(myRLDATA[,6]),split="")),ncol=4,byrow=TRUE)
M4=matrix( character(nrow(myRLDATA)*2),ncol=2)
for (i in 1:2) M4[,i]=as.character(myRLDATA[,6+i])

V=cbind(M1,M2,M3,M4)
V=as.data.frame(V)

d=ncol(V);


ALPHA=list()
for (l in 1:d) ALPHA[[l]]=table(V[,l])/nrow(V)
dimV=c()
for ( l in 1:d) {dimV[l]=length(levels(V[,l]))

V[,l]=as.numeric(V[,l])}
V=as.matrix(V)

MATCH=outer(myRLDATA$id, myRLDATA$id, FUN="==")
MATCH=MATCH[upper.tri(MATCH)]
###############################################################################


##########################CALLING THE c FUNCTION###############################
system("R CMD SHLIB rlambda11.c")
dyn.load("rlambda11.so")


rlambda=function(N,lambda=(myRLDATA$id-min(myRLDATA$id)),V,nMCMC,par.up=c(1,1,1,1),prior=c(0),up.lambda=c(0.15),g,sigma){
    tmp=.C("rlambda11",Nlambda=as.integer(N),
                     lambda=as.integer(lambda),
                     p=as.double(unlist(ALPHA)),
                     V=as.integer(c(V)-1 ),
                     dime=as.integer(dimV),
                     a=as.double(rep(0.01,nrow(V)*ncol(V))),
		     am=as.double(rep(log(0.01/0.99),ncol(V))),
                     Hdim=as.integer(ncol(V)),
                     nMCMC=as.integer(nMCMC),
                     sima=double(nMCMC*ncol(V)),
                     simlambda=integer(nMCMC*nrow(V)),
                     simnz=integer(nMCMC),
                     simNpop=integer(nMCMC),
                     n=as.integer(nrow(V)),
                     par.up=as.integer(par.up),
                     prior=as.integer(prior),
                     up.lambda=as.double(up.lambda),
                     gzeta=as.double(g),
		     sigmal=as.double(sigma),
		     simp=double(nMCMC*sum(dimV)))
    return(list(LAMBDA=matrix(tmp$simlambda,ncol=nrow(V),byrow=T),
                NZ1=c(tmp$simnz),
                Npop=c(tmp$simNpop),
                PAR=matrix(tmp$sima,ncol=ncol(V),byrow=T),
		theta=matrix(c(tmp$simp),ncol=sum(dimV),byrow=T)))}

#########################################################################################################################




##############################RUNNING THE MCMC#####################################################
nMCMC=10000
N1=2500
g=1.02
sigma=c(0.5, 0.1, 0.01)



out1=rlambda(N=N1, #initial value for N
            lambda=sample(min(N1,500),size=500,rep=TRUE), # initial partition
	    V=V,     #   observed data
	    nMCMC=nMCMC,  #number of MCMC iteration
	    par.up=c(1,1,1,1), # which elements to update. 1,1,1,1 means all the unknowns: lambda, alpha, theta, N
	    prior=c(0), #
	    up.lambda=c(0.1), #probability to update a single lambda[i,j]
	    g=g,
	    sigma=sigma)

###################################################################################################################
filename <- tempfile(pattern = "authorcode", tmpdir = "~/discussion_unified_framework", fileext = ".RData")
save(out1, file = filename)


# 
# ###################### FIGURE 2 PRODUCTION #######################################################################
# burnin=5000
# NZ=out1$NZ[-(1:burnin)]
# Npop=out1$Npop[-(1:burnin)]
# 
# 
# ####generazione a priori su k########
# rzeta=function(a){
# b=2^{a-1}
# cond=TRUE
# while(cond){
# u=runif(1)
# v=runif(1)
# x=floor(u^(-1/(a-1)))
# if (x==Inf) x=.Machine$double.xmax
# t=(1+1/x)^(a-1)
# cond=(v*x*(t-1)/(b-1)>t/b)}
# return(x)}
# 
# 
# 
# NN2=c()
# for(i in 1:100000) NN2[i]=rzeta(1.02)
# 
# k2=c()
# for(i in 1:100000) {
#      if (NN2[i]>.Machine$integer) k2[i]=length(unique(sample(.Machine$integer,size=500,rep=T)))
#      else k2[i]=length(unique(sample(NN2[i],size=500,rep=T)))}
# 
# 
# 
# 
# par(mfrow=c(2,2))
# #posteriori di k
# p.post=table(NZ)/length(NZ)
# x.post=as.numeric(names(table(NZ)))
# plot(x.post,p.post,xlim=c(410,495),xlab="k",type="b",lwd=1,cex.lab=2,main="",ylab="")
# p.prior=table(k2)/length(k2)
# x.prior=as.numeric(names(table(k2)))
# lines(x.prior,p.prior,col=2,lwd=2)
# 
# #posteriori di N
# p.post=table(Npop)/length(Npop)
# x.post=as.numeric(names(table(Npop)))
# hist(Npop,xlim=c(500,3500),xlab="N",cex.lab=2,main="",prob=TRUE,ylab="",col="orange")
# 
# 
# 
# #posteriori dei FNR E FDR
# post.rates=function(out,MATCH){
# 
# PAR=out$PAR
# LAMBDA=out$LAMBDA
# 
# 
# UM=matrix(ncol=2,nrow=nrow(LAMBDA))
# for ( j in 1:nrow(LAMBDA)){
# print(j)
# C=outer(LAMBDA[j,], LAMBDA[j,], FUN="==")
# C=C[upper.tri(C)]
# 
# 
# U=length(unique(LAMBDA[j,]))
# T=sum(C)
# 
# CL=sum(MATCH==1 & C==1) #correct link
# FN=sum(MATCH==1 & C==0) #false negative
# FP=sum(MATCH==0 & C==1) #false positive
# CNL=sum(MATCH==0 & C==0) #correct non link
# 
# FNR=FN/(CL+FN) # quanti dei veri link non vengono individuati?
# FDR=FP/(CL+FP) # quanti dei link individuati non sono corretti?
# 
# UM[j,]=c(FNR,FDR)}
# return(UM)}
# 
# UM=post.rates(out=out1,MATCH=MATCH)[-(1:burnin),]
# UM[,2]=round(UM[,2],2)
# 
# barplot(table(UM[,1])/sum(table(UM[,1])),xlab="FNR",cex.lab=2,cex.axis=2,cex.names=2,col="purple")
# barplot(table(UM[,2])/sum(table(UM[,2])),xlab="FDR",cex.lab=2,cex.axis=2,cex.names=2,col="purple")
# 
# plot(out1$Npop)
