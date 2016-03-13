#' Build single tree of a forest
#' 
#' Single tree using the training samples has been built and has been returned as a Model 
#'  
#' @param X Input Training matrix of M x N, M is the number of training samples and N is the number of features
#' @param Y Output Training response of M x T, M is the number of samples and T is number of ouput Features(Response)
#' @param mtree number of features will be used for each split
#' @param min_leaf minimum number of samples in the leaf node 
#' @param V_inv Covariance matrix of response matrix
#' @param Command 1 for RF and 2 for MRF
#' @return Model of single tree of the forest 
#' @export
#'
build_single_tree <- function(X, Y, mtree, min_leaf,V_inv,Command){
  model=rep( list(NULL), 10000 )
  i=1
  Index=1:nrow(X)
  
  model=split_node(X,Y,mtree,Index,i,model,min_leaf,V_inv,Command)
  return(model)
}