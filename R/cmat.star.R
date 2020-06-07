cmat.star <-
function(plist, CorrMat,  no.ord, no.norm){
  
if (no.norm==0 & no.ord>1) {
Sigma = IntermediateOO(plist, CorrMat)
}

if (no.norm>1 & no.ord==0) {
Sigma = CorrMat
}

if (no.norm==1 & no.ord==1) {
  if ( validate.target.cormat(plist, CorrMat,  no.ord, no.norm)) { ## target correlation and provabilities are validated here.
    
    ON = IntermediateON(plist, CorrMat[(no.ord+1):nrow(CorrMat), 1:no.ord] )
    Sigma = diag(2)
    Sigma[lower.tri((Sigma))] = ON
    Sigma = Sigma + t(Sigma)
    diag(Sigma) = 1}
}
  
if (no.norm>1 & no.ord==1) {
    if ( validate.target.cormat(plist, CorrMat,  no.ord, no.norm)) { ## target correlation and provabilities are validated here.
      
      ON = IntermediateON(plist, CorrMat[(no.ord+1):nrow(CorrMat), 1:no.ord] )
      NN = CorrMat[(no.ord+1):ncol(CorrMat), (no.ord+1):ncol(CorrMat) ]
      Sigma = cbind(c(1,ON), rbind(ON,NN) )
      
      if(!is.positive.definite(Sigma)){
        warning( "Intermediate correlation matrix is not positive definite. A nearPD function is applied.")
        Sigma=as.matrix(nearPD(Sigma, corr = TRUE, keepDiag = TRUE)$mat)
      }
      Sigma = ( Sigma+t(Sigma) )/2
    }
  }

if (no.norm==1 & no.ord>1) {
    if ( validate.target.cormat(plist, CorrMat,  no.ord, no.norm)) { ## target correlation and provabilities are validated here.
      
      OO = IntermediateOO(plist, CorrMat[1:no.ord,1:no.ord])
      ON = IntermediateON(plist, CorrMat[(no.ord+1):nrow(CorrMat), 1:no.ord] )
      Sigma = cbind(rbind(OO,ON), c(ON,1) )
                                        
      
      if(!is.positive.definite(Sigma)){
        warning( "Intermediate correlation matrix is not positive definite. A nearPD function is applied.")
        Sigma=as.matrix(nearPD(Sigma, corr = TRUE, keepDiag = TRUE)$mat)
      }
      Sigma = ( Sigma+t(Sigma) )/2
    }
  }
  
if (no.norm>1 & no.ord>1) {
if ( validate.target.cormat(plist, CorrMat,  no.ord, no.norm)) { ## target correlation and provabilities are validated here.

OO = IntermediateOO(plist, CorrMat[1:no.ord,1:no.ord])
ON = IntermediateON(plist, CorrMat[(no.ord+1):nrow(CorrMat), 1:no.ord] )
NN = CorrMat[(no.ord+1):ncol(CorrMat), (no.ord+1):ncol(CorrMat) ]
Sigma = cbind(rbind(OO,ON), rbind(t(ON),NN) )

if(!is.positive.definite(Sigma)){
warning( "Intermediate correlation matrix is not positive definite. A nearPD function is applied.")
Sigma=as.matrix(nearPD(Sigma, corr = TRUE, keepDiag = TRUE)$mat)
}
Sigma = ( Sigma+t(Sigma) )/2
}
}
rownames(Sigma)<-NULL
return(Sigma)

}
