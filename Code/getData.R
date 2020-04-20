#######
# run this file to get the data 

##################DATA PREPARATION############################
rm(list = ls())
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
# n <- 20
# V <- V[1:n,]

ALPHA=list()
for (l in 1:d) ALPHA[[l]]=table(V[,l])/nrow(V)
dimV=c()
for ( l in 1:d) {
  dimV[l]=length(levels(V[,l]))
  V[,l]=as.numeric(V[,l])
}
V=as.matrix(V)
H <- ncol(V)
# V is a matrix of n * 14
# entries in column l takes values in 1:Vl (in the code the authors use h for l)

MATCH=outer(myRLDATA$id, myRLDATA$id, FUN="==")
MATCH=MATCH[upper.tri(MATCH)]
