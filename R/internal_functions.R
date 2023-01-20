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


#' Local similarity matrix
#'
#' @description Computation of the local similarity matrix based on the first
#' KNN elements.
#' It is a modified version of the .dominateset function implemented in SNFtool
#' library (https://cran.r-project.org/web/packages/SNFtool/index.html) to avoid
#' division by zero when the "ignore" version of the miss-SNF algorithm is used.
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
#
#'
#' @param xx matrix. Square matrix of pairwise similarities.
#' @param KK integer. Number of nearest neighbours.
#'
#' @return Local similarity matrix.
#' @export
#'
#' @examples
#' # Create a matrix
#' set.seed(123);
#' M1 <- matrix(runif(100, min = 0, max = 1), nrow = 10); # 10 samples
#' rownames(M1) <- paste0("ID_", 1:nrow(M1));
#'
#' # Compute similarity matrix
#' sim <- scaled.exp.euclidean(M1, kk=3, sigma=0.5);
#' # Compute local similarity matrix
#' loc.sim <- missSNF::.local.similarity.matrix(sim, 5);
#'
#' @keywords internal
#' @noRd
.local.similarity.matrix <- function(xx,KK=20) {

    zero <- function(x) {
        s = sort(x, index.return=TRUE)
        x[s$ix[1:(length(x)-KK)]] = 0
        return(x)
    }

    normalize <- function(X) {
        row.sum <- rowSums(X);
        row.sum[row.sum == 0] <- 1;
        return(X / row.sum)
    }

    A = matrix(0,nrow(xx),ncol(xx));
    for(i in 1:nrow(xx)){
        A[i,] = zero(xx[i,]);
    }

    return(normalize(A))
}

