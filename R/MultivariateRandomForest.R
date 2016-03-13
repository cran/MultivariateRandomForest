#' Build Random or Multivariate Random Forest and Predicttion for Testing Samples
#' 
#' Random or Multivariate Random Forest using Training samples is built and the model is used to do
#' the prediction of Testing Samples
#'  
#' @param trainX Input matrix of M x N, M is the number of training samples and N is the number of features 
#' @param trainY Output response of M x T, M is the number of samples and T is number of ouput Features(Response)
#' @param n_tree number of trees in the forest
#' @param mtree number of features will be used for each split
#' @param min_leaf minimum number of samples in the leaf node 
#' @param testX Testing samples of Q x N, Q is the number of testing samples and N is the number of features(same order and
#' size used as training)
#' @param Command 1 for RF and 2 for MRF
#' @return Predicted response of Testing samples. For RF, prediction of the specific feature and for MRF. prediction of all the output features.
#' @examples
#' trainX=matrix(runif(50*100),50,100)
#' trainY=matrix(runif(50*5),50,5)
#' n_tree=5
#' mtree=10
#' min_leaf=5
#' testX=matrix(runif(10*100),10,100)
#' Command=2#2 for MRF method
#' #Prediction size is 10 x 5, where 10 is the number 
#' #of testing samples and 5 is the number of output features
#' Prediction=MultivariateRandomForest(trainX, trainY, n_tree, mtree, min_leaf, testX,Command)
#' 
#' @references Breiman, Leo. "Random forests." Machine learning 45.1 (2001): 5-32.
#' @export
#' 
MultivariateRandomForest <- function(trainX, trainY, n_tree, mtree, min_leaf, testX,Command){
  theta <- function(trainX){trainX}
  results <- bootstrap::bootstrap(1:nrow(trainX),n_tree,theta) #no indics, gives number
  b=results$thetastar
  Variable_number=ncol(trainY)
  Y_HAT=matrix(  0*(1:Variable_number*nrow(testX)),  ncol=Variable_number,   nrow=nrow(testX)  )
  Y_pred=NULL
  
  for (i in 1:n_tree){
    Single_Model=NULL
    X=trainX[ b[ ,i],  ]
    Y=matrix(trainY[ b[ ,i],  ],ncol=Variable_number)
    V_inv = solve(stats::cov(Y)) # calculate the V inverse
    Single_Model=build_single_tree(X, Y, mtree, min_leaf,V_inv,Command)
    Y_pred=single_tree_prediction(Single_Model,testX,Variable_number)
    for (j in 1:Variable_number){
      Y_HAT[,j]=Y_HAT[,j]+Y_pred[,j]
    }
  }
  Y_HAT=Y_HAT/n_tree
  return(Y_HAT)
}