## this script essentially runs the authors' code and saves the output
## requires 'rlambda11modified.c' to be in the working directory

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
system("R CMD SHLIB rlambda11modified.c")
dyn.load("rlambda11modified.so")


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
nMCMC=1e4
N1=2500
g=1.02
sigma=c(0.5, 0.1, 0.01)


library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
##
print("update of eta, theta and N")
## run with updates of eta, theta and N 
## 5 chains independently
authors_runs <- foreach(irep = 1:5) %dorng% {
  out1=rlambda(N=N1, #initial value for N
               lambda=sample(min(N1,500),size=500,rep=TRUE), # initial partition
               V=V,     #   observed data
               nMCMC=nMCMC,  #number of MCMC iteration
               par.up=c(1,0,1,1), # which elements to update. 1,1,1,1 means all the unknowns: lambda, alpha, theta, N
               prior=c(0), #
               up.lambda=c(0.1), #probability to update a single lambda[i,j]
               g=g,
               sigma=sigma)
}
save(nMCMC, authors_runs, file = "authorcodemodified_eta_theta_N.RData")

