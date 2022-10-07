# Internal functions needed by miss-SNF



#' Normalization method for affinity matrices
#' 
#' @description Normalization method used to compute a "global" similarity
#' matrix in SNF algorithm, capturing the overall relationships between patients.
#' It corresponds to the formula:
#' 
#' \loadmathjax
#' \mjsdeqn{P(i,j) = 
#' \begin{cases}
#' \frac{W(i,j)}{2 \sum_{k \neq i} W(i,k)} & \text{, if $j \neq i$}\\
#' 
#' 1/2 & , \text{if $j = i$}
#' \end{cases}}
#' 
#' NOTE: This function is taken from SNFtool 
#' (https://cran.r-project.org/web/packages/SNFtool/index.html), but it was
#' defined inside the main function SNFtool::SNF. Here it is used as internal
#' function to make the code more readable.
#'
#' @param X matrix. Similarity matrix to normalize.
#'
#' @return Normalized matrix.
.normalize <- function(X){
    
    row.sum.mdiag <- rowSums(X) - diag(X)
    #If rowSumx(X) == diag(X), set row.sum.mdiag to 1 to avoid div by zero
    row.sum.mdiag[row.sum.mdiag == 0] <- 1
    X <- X/(2*(row.sum.mdiag))
    diag(X) <- 0.5
    
    return(X)
}
