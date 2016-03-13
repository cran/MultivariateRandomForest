#' single_tree_prediction
#' 
#' Using the model of the tree, for testing samples, prediction has been calculated
#'  
#' @param Single_Model Model of a single particular tree
#' @param X_test Testing samples of Q x N, Q is the number of testing samples and N is the number of features(same order and
#' size used as training) 
#' @param Variable_number Number of drugs which response needs to be calculated 
#' @return predicted response of Testing samples for a single tree
#' @details A model contrains for each split of the node what criteria has been used. For the testing samples, using the criteria
#' testing sample reaches a leaf node, which have some response value stored. Average of these values are used as prediction
#' for the testing samples
#' @export
single_tree_prediction <- function(Single_Model,X_test,Variable_number){
  
  
  Y_pred=matrix(  0*(1:nrow(X_test)*Variable_number)  ,nrow=nrow(X_test),  ncol=Variable_number)
  
  for (k in 1:nrow(X_test)){
    xt=X_test[k, ]
    i=1
    Y_pred[k,]=predicting(Single_Model,i,xt,Variable_number)
    
  }
  #Y_pred1=unlist(Y_pred, recursive = TRUE)
  #Y_pred1=matrix(Y_pred1,nrow=nrow(X_test))
  return(Y_pred)
}