#' split_node
#' 
#' The Splitting criteria in a node of a tree has been calculated 
#'  
#' @param X Input Training matrix of M x N, M is the number of training samples and N is the number of features
#' @param Y Output Training response of M x T, M is the number of samples and T is number of ouput Features(Response)
#' @param mtree number of features will be used for each split
#' @param Index Index of training samples
#' @param i number of split
#' @param model The list of numbers where different spliting criteria has been stored
#' @param min_leaf minimum number of samples in the leaf node 
#' @param V_inv Covariance matrix of response matrix
#' @param Command 1 for RF and 2 for MRF
#' @return Model Updated Model after split has been done
#' @export
split_node <- function(X,Y,mtree,Index,i,model,min_leaf,V_inv,Command){
  ii=NULL
  Index_left=NULL
  Index_right=NULL
  if(length(Index)>min_leaf){ #create problem with 2
    Result = spliting(X,Y,mtree,Index,V_inv,Command)
    Index_left=Result[[1]]
    Index_right=Result[[2]]
    if(i==1){
      Result[[5]]=c(2,3)
    }else{
      j=1
      while (length(model[[j]])!=0){
        j=j+1
      }
      Result[[5]]=c(j,j+1)
    }
    
    model[[i]]=Result
    k=i
    i=1 #maybe unnecessary
    while (length(model[[i]])!=0){
      i=i+1
    } 
    model[[Result[[5]][1]]]=Result[[1]]
    model[[Result[[5]][2]]]=Result[[2]]
    
    model=split_node(X,Y,mtree,Index_left,model[[k]][[5]][1],model,min_leaf,V_inv,Command)#model[[i]][[5]][1]
    
    model=split_node(X,Y,mtree,Index_right,model[[k]][[5]][2],model,min_leaf,V_inv,Command)
    
    
  }else{
    ii[[1]]=matrix(Y[Index,],ncol=ncol(Y))
    model[[i]]=ii
  }
  
  
  return(model)
}