##' fit a classification model with "LDA", "KNN" or "SVM"
##'
##' @param data.train traing data that is used to build the classifier.
##' @param data.test testing data that need to be classified.
##' @param class.train Label of training data.
##' @param method "LDA", "KNN" or "SVM", default is "LDA".
##' @param K parameter for k-nearest neighbor method.
##' @return a list with pred.test, pred.score, pred.train and
##' error.train. when method is "LDA", will also return the scaling
##' vector "sca" that could be used to generate the Linear
##' Discriminants.
##' @export
classification <- function(data.train, data.test, class.train, 
                           method=c("LDA", "KNN", "SVM"), K=50){
  
  method <- match.arg(method)
  
  if(method == "LDA"){
    suppressMessages(require(MASS))|| stop("the required package 'MASS' is not installed. ")
    fit <- lda(data.train, grouping=class.train)
    pred <- predict(fit, data.test)
    pred.test <- pred$class
    names(pred.test) <- rownames(data.test)
    pred.score <- apply(pred$posterior, 1, max)
    pred.train <- predict(fit)$class
    error.train <- mean(pred.train != class.train)
    sca <- fit$scaling
    return(list(pred.test=pred.test, pred.score=pred.score, 
                pred.train=pred.train, error.train=error.train, sca=sca))
  }else if(method == "KNN"){
    suppressMessages(require(class))|| stop("the required package 'class' is not installed. ")
    knnfit <- knn(data.train, data.test, cl = class.train, k = K, prob=TRUE)
    pred.test <- knnfit
    names(pred.test) <- rownames(data.test)
    pred.score <- attr(knnfit, "prob")
    pred.train <- knn(data.train, data.train, cl = class.train, k = K, prob=TRUE)
    error.train <- mean(pred.train != class.train)
    return(list(pred.test=pred.test, pred.score=pred.score, 
                pred.train=pred.train, error.train=error.train))
  }else if(method == "SVM"){
    suppressMessages(require(e1071))|| stop("the required package 'e1071' is not installed. ")
    fit <- svm(x=data.train, y=as.factor(class.train), probability=TRUE)
    pred.test <- predict(fit, data.test, probability=TRUE)
    names(pred.test) <- rownames(data.test)
    pred.score <- apply(attr(pred.test, "probabilities"), 1, max)
    pred.train <- predict(fit)
    error.train <- mean(pred.train != class.train)
    return(list(pred.test=pred.test, pred.score=pred.score, 
                pred.train=pred.train, error.train=error.train))
  }
  
}  
