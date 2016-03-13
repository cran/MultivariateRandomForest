#' Prediction of testing sample
#' 
#' For a specific node, which way should the testing sample goes has been determined using spliting criteria stored in the model. 
#'  
#' @param Single_Model Model of a single particular tree
#' @param i number of which split needs to be done
#' @param X_test Testing samples of Q x N, Q is the number of testing samples and N is the number of features(same order and
#' size used as training) 
#' @param Variable_number How many drugs are given for computation
#' @return Predicted response of a Testing samples 
#' @export
predicting <- function(Single_Model,i,X_test,Variable_number){
  
  Result=NULL
  
  if(length(Single_Model[[i]])==5){
    feature_no=Single_Model[[i]][[3]]
    feature_value=X_test[feature_no]
    if(feature_value<Single_Model[[i]][[4]]){  #feature value less than threshold value
      #i=i*2+1
      Result=predicting(Single_Model,Single_Model[[i]][[5]][1],X_test,Variable_number)
    }else{                                    #feature value greater than threshold value
      #i=i*2+2
      Result=predicting(Single_Model,Single_Model[[i]][[5]][2],X_test,Variable_number)
    }
  }else{
    Result=matrix(  0*Variable_number,  ncol=Variable_number)
    if (Variable_number>1){
      for (jj in 1:Variable_number){
        Result[,jj]=mean(Single_Model[[i]][[1]][,jj])
      }
    }else {
      for (jj in 1:Variable_number){
        Result[,jj]=mean(unlist(Single_Model[[i]][[1]]))
      }
    }
    
  }
  return(Result)
}