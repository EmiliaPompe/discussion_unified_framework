#'@export
get_RLdata500 <- function(){
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
  p=ncol(V);
  fieldfrequencies=list()
  for (l in 1:p) fieldfrequencies[[l]]=table(V[,l])/nrow(V)
  Mvec <- c()
  for (l in 1:p) {
    Mvec[l]=length(levels(V[,l]))
    V[,l]=as.numeric(V[,l])
  }
  V=as.matrix(V)
  # V is a matrix of n * p
  # entries are values in 1:M_l, where M_l is the number of possible categories in field l
  # and field l in in 1:p
  # Vector 'Mvec' has entries 'M_l' for l in 1:p
  return(list(V = V, Mvec = Mvec, fieldfrequencies = fieldfrequencies))
}