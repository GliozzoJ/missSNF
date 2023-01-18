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
#' @keywords internal
#' @noRd
.normalize <- function(X){

    row.sum.mdiag <- rowSums(X) - diag(X)
    #If rowSumx(X) == diag(X), set row.sum.mdiag to 1 to avoid div by zero
    row.sum.mdiag[row.sum.mdiag == 0] <- 1
    X <- X/(2*(row.sum.mdiag))
    diag(X) <- 0.5

    return(X)
}


#' Compute local similarity matrix
#'
#' @description This function computes the local similarity matrix that
#' takes into account local affinity through K-Nearest Neighbor.
#' It corresponds to the formula:
#'
#' \loadmathjax
#' S(i,j) =
#' \begin{cases}
#'      \frac{W(i,j)}{\sum_{k \in N_i} W(i,k)} & \text{, if $j \in N_i$}\\
#'      0 & , \text{otherwise}\\
#' \end{cases}
#'
#' where $N_i = \{ x_k | x_k \in kNN(x_i) \cup \{ x_i \}\}$.
#' Note: this function corresponds to the internal function ".dominateset" in
#' SNFtool package (https://cran.r-project.org/web/packages/SNFtool/index.html)
#' but we modified the normalization of the function to avoid division by zero
#' when miss-SNF ZERO (i.e. mode="ignore") is used.
#'
#' @param xx matrix. Similarity matrix.
#' @param KK numeric. Number of neighbors in K-nearest neighbors.
#'
#' @return Local similarity matrix.
#' @keywords internal
#' @noRd
.localSimMat <- function(xx, KK=20) {

    zero <- function(x) {
        s = sort(x, index.return=TRUE)
        x[s$ix[1:(length(x)-KK)]] = 0
        return(x)
    }

    # Modified "normalize" function to avoid division by
    # zero for missing patients when mode="ignore"
    normalize <- function(X) {

        rSums <- rowSums(X)
        rSums[rSums == 0] <- 1
        X <- X / rSums

        return(X)
    }


    A = matrix(0,nrow(xx),ncol(xx));
    for(i in 1:nrow(xx)){
        A[i,] = zero(xx[i,]);
    }

    return(normalize(A))
}





