#' spliting
#' 
#' Spliting of training samples in a node using different feature vector has been calculated
#'  
#' @param X Input Training matrix of M x N, M is the number of training samples and N is the number of features
#' @param Y Output Training response of M x T, M is the number of samples and T is number of ouput Features(Response)
#' @param mtree number of features will be used for each split
#' @param Index Index of training samples
#' @param V_inv Covariance matrix of response matrix
#' @param Command 1 for RF and 2 for MRF
#' @return Spliting criteria for the splitting of the node
#' @details In a node of training samples, which spliting gives best cost reduction has been calculted. mtree number of feature
#' are considered and for a specific feature and specific split point, leaf node and right node split cost has been calculated.
#' Whichever feature and split has been given minimum cost has been selected and the criteria are returned
#' @export
spliting <- function(X,Y,mtree,Index,V_inv,Command){
  x=X[Index, ]
  if (Command==1){
    y=matrix(Y[Index, ],ncol=1)
  }else {
    y=Y[Index, ]
  }
  f = ncol(x) # number of features
  f_1 =sort(sample(f, mtree)) #randomly taken 10 features, for each splits vary
  N = nrow(x)
  min_score =NULL
  DL=matrix(data=NA, nrow=(N-1), ncol=length(f_1))
  DR=matrix(data=NA, nrow=(N-1), ncol=length(f_1))
  id =NULL
  for(j in 1:length(f_1)){
    xj = x[ ,f_1[j]]
    tmp=sort(xj,index.return=TRUE)   # sort the values
    s=tmp$x
    idx=tmp$ix 
    id2=NULL
    for (k1 in 1:(N-1)){
      k=k1
      yk_left=matrix(y[idx[1:k1],],ncol=ncol(y)) # will create problem with mrf
      yk_right=matrix(y[idx[(k1+1):length(idx)],],ncol=ncol(y))
      D_left= Multi_D_mod(yk_left,V_inv,Command)
      D_right= Multi_D_mod(yk_right,V_inv,Command)
      D=D_left+D_right
      DL[k,j]=D_left
      DR[k,j]=D_right
      if(j==1 && k==1){
        min_score = D
        which_feature = f_1[1]
        indexX=c(k,j)
        threshold_feature = (s[k] + s[k+1])/2
        ij_i = 1
        ij_j = 1
      }
      
      if(D< min_score){
        min_score = D
        which_feature = f_1[j]
        indexX=c(k,j)
        threshold_feature = (s[k] + s[k+1])/2
        ij_i = k
        ij_j = j
      }
      id2[[k1]]= rep(NA, times=(N-1))
      id2[[k1]] = idx[1:k1]
    }
    id[[j]] = id2
  }
  DD=rbind(DL,DR)
  index_left =id[[ij_j]][[ij_i]] # column, then row
  index_right = 1:N
  index_right=index_right[-index_left] 
  
  index_left=Index[sort(index_left)]
  index_right = Index[sort(index_right)]
  #IN=rbind(index_left, index_right)
  
  result=list(index_left, index_right, which_feature, threshold_feature)
  return(result)
}