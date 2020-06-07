IntermediateON <-
function(plist, ONCorrMat){

Cor_NNforON<-function(pvec, ON.cor) {
Z=rnorm(100000,0,1)
Y=rnorm(100000,0,1)
YORD = ordinalize(pvec,Y)
c = cor(YORD[order(YORD)],Z[order(Z)])/cor(Y[order(Y)],Z[order(Z)]) 
r = ON.cor/c 
return(r)
}

if ( length(plist)==1 & is.matrix(ONCorrMat)==FALSE ) {
  cmat.corrected = Cor_NNforON(plist[[1]], ONCorrMat)
}

if ( length(plist)>1 & is.matrix(ONCorrMat)==FALSE ) {
  for (i in 1:length(plist)) {
    cmat.corrected[i] = Cor_NNforON(plist[[i]], ONCorrMat[i])
  }
}

if ( length(plist)>1 & is.matrix(ONCorrMat)==TRUE ) {
  cmat.corrected = matrix(NA,nrow(ONCorrMat), ncol(ONCorrMat) )
  for (j in 1:ncol(ONCorrMat)){
    for (i in 1:nrow(ONCorrMat)){
      cmat.corrected[i,j] = Cor_NNforON(plist[[j]], ONCorrMat[i,j] )
    }
  }
}

return(cmat.corrected)
}
