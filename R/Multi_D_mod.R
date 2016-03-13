#' Information Gain or splitting cost
#' 
#' For a specific spliting in a node for a spcific feature vector, cost has been calculated. 
#'  
#' @param y Feature Vector
#' @param V_inv Covariance matrix of response matrix
#' @param Command 1 for RF and 2 for MRF
#' @return threshold_value 
#' @export
Multi_D_mod  <- function(y,V_inv,Command){
  NN=nrow(y)
  ybar5=colSums (y, na.rm = FALSE, dims = 1)/nrow(y)
  ybar2=matrix(ybar5,nrow=1,ncol=ncol(y))
  
  if (Command==2){ #using VMRF
    yhat=y-kronecker(matrix(1,nrow(y),1),ybar2)
    
    D = sum(diag(yhat %*% V_inv %*% t(yhat)))
    
  }else if(Command==1){  #using rf
    ybar=sum(y)/nrow(y)
    D=sum((y-ybar)^2)
  }
  return(D)
}