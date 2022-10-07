#' miss-SNF: integration of patients' networks considering missing cases
#'
#' @description Extension of the algorithm Similarity Network Fusion
#' (\url{https://doi.org/10.1038/nmeth.2810}) able to handle the absence (complete or
#' nearly complete) of a specific data source for a given patient by:
#' \loadmathjax
#' \itemize{
#' \item "reconstruct strategy" that can partially reconstruct missing
#'   data by using information from different sources
#' \item "ignore strategy" which simply ignores missing data during the integration
#'   process
#' }
#'
#' @param Mall list of named matrices/dataframes (samples x features).
#' @param sims vector of strings. It is a vector containing the names of the
#'             similarity measures to apply to the matrices in Mall.
#'             "scaled.exp.euclidean" is the scaled exponential euclidean distance;
#'             "scaled.exp.chi2" is the scaled exponential chi-square distance
#' @param mode string. If you want to partially reconstruct missing data use
#'             "reconstruct", otherwise ignore them during integration using
#'             "ignore".
#' @param perc.na percentage of NAs above which patient is considered missing.
#' @param miss.symbols vector of strings. If not NULL, the provided symbols
#'                     in matrices are converted to NA.
#' @param K Number of neighbors in K-nearest neighbors part of the algorithm.
#' @param t Number of iterations for the diffusion process.
#' @param impute string. Kind of imputation method to apply in case of samples
#'               with few missing values. Options are: NULL, "mean", "median".
#'               NOTE: if impute=NULL, then perc.na has to be 0.
#' @param ... additional arguments for similarity measures.
#'
#' @return A list with two elements:
#' \itemize{
#' \item W : integrated similarity matrix.
#' \item removed.pts : vector with names of removed patients (i.e. patients
#' present only in one matrix of Mall and having too much NAs).
#' }
#' @export
#'
#' @examples
#'
#' # Create list of imput matrices
#' set.seed(123);
#' M1 <- matrix(runif(50, min = 0, max = 1), nrow = 5); # 5 samples
#' rownames(M1) <- paste0("ID_", 1:nrow(M1));
#' M1[1, 1:4] <- c(NA, NA, NA, NA);
#' M1[4, ] <- rep(NA, ncol(M1));
#'
#' M2 <- matrix(runif(80, min = 0, max = 1), nrow = 8); # 8 samples
#' rownames(M2) <- paste0("ID_", 3:10);
#' M2[2, ] <- rep(NA, ncol(M2));
#' M2[4, ] <- rep(NA, ncol(M2));
#' M2[8, 1] <- NA;
#'
#' M3 <- matrix(runif(80, min = 0, max = 1), nrow = 8); # 8 samples
#' rownames(M3) <- c(paste0("ID_", 1:5), paste0("ID_", 11:13));
#' M3[4, 1] <- NA;
#' M3[7, 1:8] <- rep(NA, 8);
#'
#' Mall <- list("M1"=M1, "M2"=M2, "M3"=M3);
#'
#' # Call miss.snf using "reconstruct strategy"
#' W.r <- miss.snf(Mall, sims=rep("scaled.exp.euclidean", 3), mode="reconstruct",
#'                K=3, kk=2);
#'
#' # Call miss.snf using "ignore strategy"
#' W.i <- miss.snf(Mall, sims=rep("scaled.exp.euclidean", 3), mode="ignore",
#'                K=3, kk=2);
miss.snf <- function(Mall, sims, mode="reconstruct", perc.na=0.2,
                     miss.symbols=NULL, K=20, t=20, impute="median", ...) {


    # Check that a similarity measures is provided for each matrix
    if(length(Mall) != length(sims)){

        stop("miss.snf: number of matrices to integrate has to be equal to
             number of similarity measures provided.");
    }

    # Get list of missing patients for each matrix
    miss.pts <- get.miss.pts(Mall, perc.na=perc.na, miss.symbols=miss.symbols);

    # Convert miss.symbols to NA
    if(!is.null(miss.symbols)){
        for(i in 1:length(Mall)){
            for(j in 1:length(miss.symbols)){
                Mall[[i]][Mall[[i]] == miss.symbols[j]] <- NA;
            }
        }
    }

    # Convert Mall elements to matrix
    for (i in 1:length(Mall)){
        Mall[[i]] <- as.matrix(Mall[[i]]);
        class(Mall[[i]]) <- "numeric";
    }

    # Impute patients' features with NAs
    if(is.null(impute)){
        if(perc.na != 0){

            stop("miss.snf: set perc.na=0 to remove patients having NAs and avoid imputation.")
        }
    }else if(impute == "median"){
        for(i in 1:length(Mall)){
            Mall[[i]] <- apply(Mall[[i]], 2, function(x) ifelse(is.na(x), median(x, na.rm=T), x))
        }
    } else if(impute == "mean"){
        for(i in 1:length(Mall)){
            Mall[[i]] <- apply(Mall[[i]], 2, function(x) ifelse(is.na(x), mean(x, na.rm=T), x))
        }
    }

    # Compute similarity matrices W (without considering missing patients)
    Wall <- list();
    for(i in 1:length(Mall)){

        idx.miss.pts <- which(rownames(Mall[[i]]) %in% miss.pts[[i]]); #idx patients with NAs
        sim.fun <- get(sims[[i]]);

        # Handle presence of matrices without patients to remove for NAs
        if(length(idx.miss.pts) != 0){
            Wall[[i]] <- sim.fun(Mall[[i]][-idx.miss.pts, ], ...);
        } else{
            Wall[[i]] <- sim.fun(Mall[[i]], ...);
        }
    }

    # Align networks replacing missing patients with 0
    Wall_aligned <- NetInt::lalign.networks(fill=0, Wall);

    # Get names of patients removed from matrices because present in only
    # one matrix but having too much NAs
    removed.pts.tab <- table(unlist(miss.pts)) == length(Mall);
    removed.pts <- rownames(removed.pts.tab)[removed.pts.tab];

    # Set diagonal of missing elements to 1 if "reconstruct" mode is used
    # If "ignore" mode, missing patients are already set to vector of zeros
    idx.miss.aligned <- list();
    for(i in 1:length(Wall_aligned)){
        idx <- which(!(rownames(Wall_aligned[[i]]) %in% rownames(Wall[[i]])));
        idx.miss.aligned[[i]] <- idx;
        if(mode == "reconstruct"){

            diag(Wall_aligned[[i]])[idx] <- 1;

        }
    }

    # Use Wall_aligned in SNF code
    Wall <- Wall_aligned;

    # FROM NOW ON STARTS SNF CODE
    ############################################################################
    ############################################################################

    wall.names <- dimnames(Wall[[1]])

    LW <- length(Wall)

    #Normalize different networks to avoid scale problems.
    #### [J]: Wall after normalization contains "global" similarity matrices P
    newW <- vector("list", LW)
    nextW <- vector("list", LW)
    for(i in 1:LW){
        Wall[[i]] <- .normalize(Wall[[i]])
        Wall[[i]] <- (Wall[[i]] + t(Wall[[i]]))/2

        ### [J]: set diagonal to 1 if mode="reconctruct", 0 otherwise for missing
        ### patients. The other elements for missing patients should remain zero
        ### even after normalization.

        if(mode == "reconstruct"){

            diag(Wall[[i]])[idx.miss.aligned[[i]]] <- 1;
        } else {

            diag(Wall[[i]])[idx.miss.aligned[[i]]] <- 0;
        }
    }

    ### Calculate the local transition matrix. (KNN step?)
    #### [J]: newW contains the local transition matrices S
    #### [J]: since missing patients has no neighbours by design
    #### (except for itself), in case of recostruction strategy the weights
    #### are already set to zero and the diagonal to 1. However,
    #### if the strategy is "ignore", the computation leads to NaN
    #### (due to division by zero) for missing patients
    for(i in 1:LW){
        newW[[i]] <- (SNFtool:::.dominateset(Wall[[i]], K))

        if (mode == "ignore"){
            # Set to zero rows of missing patients
            newW[[i]][idx.miss.aligned[[i]], ] <- rep(0, ncol(newW[[i]]));
        }
    }

    #Perform the diffusion for t iterations
    for (i in 1:t) {
        for(j in 1:LW){
            sumWJ <- matrix(0,dim(Wall[[j]])[1], dim(Wall[[j]])[2])
            for(k in 1:LW){
                if(k != j) {
                    sumWJ <- sumWJ + Wall[[k]]
                }
            }
            nextW[[j]] <- newW[[j]] %*% (sumWJ/(LW-1)) %*% t(newW[[j]])
        }

        #Normalize each new obtained networks.
        for(j in 1 : LW){
            Wall[[j]] <- .normalize(nextW[[j]])
            Wall[[j]] <- (Wall[[j]] + t(Wall[[j]]))/2;
        }
    }

    # Construct the combined affinity matrix by summing diffused matrices
    W <- matrix(0, nrow(Wall[[1]]), ncol(Wall[[1]]))
    for(i in 1:LW){
        W <- W + Wall[[i]]
    }

    W <- W/LW
    W <- .normalize(W)
    W <- (W + t(W)) / 2

    # Assign names to matrix
    dimnames(W) <- wall.names

    return(list(W=W, removed.pts=removed.pts))
}


################################################################################
############## Similarity measures #############################################
################################################################################

#' Scaled exponential euclidean distance
#'
#' @description It is a wrapper function to directly compute the scaled
#' exponential euclidean distance as implemented in SNFtool library.
#'
#' NOTE: from SNFtool "If the data is continuous, we recommend to use
#' the function dist2 as follows:"
#'
#' dist <- (dist2(as.matrix(Data1),as.matrix(Data1)))^(1/2);
#'
#' This is because the function dist2 actually computes the squared
#' euclidean distance (without the square root).
#'
#' @param M matrix/dataframe (samples x features).
#' @param kk integer. Number of nearest neighbours to sparsify similarity.
#' @param sigma numeric. Variance for local model.
#'
#' @return similarity matrix computed using scaled exponential euclidean
#'         distance.
#' @export
#'
#' @examples
#' # Create a matrix
#' set.seed(123);
#' M1 <- matrix(runif(50, min = 0, max = 1), nrow = 5); # 5 samples
#' rownames(M1) <- paste0("ID_", 1:nrow(M1));
#'
#' # Compute similarity matrix
#' sim <- scaled.exp.euclidean(M1, kk=3, sigma=0.5);
scaled.exp.euclidean <- function(M, kk=20, sigma=0.5){

    dist <- (SNFtool::dist2(as.matrix(M), as.matrix(M)))^(1/2);
    sim <- SNFtool::affinityMatrix(dist, K=kk, sigma=sigma);

    return(sim)
}



#' Scaled exponential chi-square distance
#'
#' @description It is a wrapper function to directly compute the scaled
#' exponential chi-square distance as implemented in SNFtool library.
#'
#' @param M matrix/dataframe (samples x features).
#' @param kk integer. Number of nearest neighbours to sparsify similarity.
#' @param sigma numeric. Variance for local model.
#'
#' @return similarity matrix computed using scaled exponential chi-square
#'         distance.
#' @export
#'
#' @examples
#'
#' # Create a matrix
#' set.seed(123);
#' M1 <- rbind(c(1,2,1,1,1,2,2), c(2,2,2,1,1,1,1), c(1,2,1,2,1,2,1),
#'            c(2,2,2,1,2,1,2));
#'
#' rownames(M1) <- paste0("ID_", 1:nrow(M1));
#'
#' # Compute similarity matrix
#' sim <- scaled.exp.chi2(M1, kk=3);
scaled.exp.chi2 <- function(M, kk=20, sigma=0.5){

    dist <- SNFtool::chiDist2(M);
    sim <- SNFtool::affinityMatrix(dist, K=kk, sigma=sigma);

    return(sim)
}


###############################################################################
############## UTILITIES ######################################################
###############################################################################

#' Get list of missing patients for each matrix
#'
#' @description This function takes in input a list of matrices/dataframes Mall
#' and returns for each one the list of patients that are missing (considering
#' the union of patients in all matrices). A patient is considered "missing"
#' if (I) it is present in the matrix but has |NA| > perc.na, (II) a patient is
#' not present in one matrix but it is present in at least one of the others.
#'
#' @param Mall list of named matrices/dataframes (samples x features).
#' @param perc.na percentage of NAs above which patient is considered missing.
#' @param miss.symbols vector of strings. If not NULL, the provided symbols in
#'                     dataframes is converted to NA.
#'
#' @return list of vectors. Each vector contains the indices of missing samples
#'         in a specific matrix/dataframe.
#' @export
#'
#' @examples
#'
#' # Create list of imput matrices
#' set.seed(123);
#' M1 <- matrix(runif(50, min = 0, max = 1), nrow = 5); # 5 samples
#' rownames(M1) <- paste0("ID_", 1:nrow(M1));
#' M1[1, 1:4] <- c(NA, NA, NA, NA);
#' M1[4, ] <- rep(NA, ncol(M1));
#'
#' M2 <- matrix(runif(80, min = 0, max = 1), nrow = 8); # 8 samples
#' rownames(M2) <- paste0("ID_", 3:10);
#' M2[2, ] <- rep(NA, ncol(M2));
#' M2[4, ] <- rep(NA, ncol(M2));
#' M2[8, 1] <- NA;
#'
#' M3 <- matrix(runif(80, min = 0, max = 1), nrow = 8); # 8 samples
#' rownames(M3) <- c(paste0("ID_", 1:5), paste0("ID_", 11:13));
#' M3[4, 1] <- NA;
#' M3[7, 1:8] <- rep(NA, 8);
#'
#' Mall <- list("M1"=M1, "M2"=M2, "M3"=M3);
#'
#' # Call function
#'
#' miss.pts <- get.miss.pts(Mall, perc.na=0.2, miss.symbols=NULL);
#'
#' # If missing values are encoded with some symbols (e.g. "?"),
#' # the function can take this situation into account
#'
#' M3[1, 3:4] <- "?"
#' Mall <- list("M1"=M1, "M2"=M2, "M3"=M3);
#'
#' miss.pts <- get.miss.pts(Mall, perc.na=0.2, miss.symbols=c("?"));
get.miss.pts <- function(Mall, perc.na=0.2, miss.symbols=NULL){

    # Check if matrices are named
    for (i in 1:length(Mall)){
        if (is.null(rownames(Mall[[i]]))){
            stop("get.miss.pts: matrices must be named")
        }
    }

    # Convert miss.symbols into NA
    if(!is.null(miss.symbols)){
        for(i in 1:length(Mall)){
            for(j in 1:length(miss.symbols)){
                Mall[[i]][Mall[[i]] == miss.symbols[j]] <- NA;
            }
        }
    }

    # Get all patients' names
    names <- lapply(Mall, rownames);
    names <- unique(unlist(names));

    # Get list of missing patients for each matrix
    miss.pts <- list();
    for(i in 1:length(Mall)){

        # Get patients not present in matrix
        miss.pts[[i]] <- setdiff(names, rownames(Mall[[i]]));

        # Get patients that are all NAs or NAs > perc.na
        n.na <- rowSums(is.na(Mall[[i]]));
        pts.many.na <- rownames(Mall[[i]])[(n.na/ncol(Mall[[i]])) > perc.na];
        miss.pts[[i]] <- c(miss.pts[[i]], pts.many.na);
    }

    return(miss.pts)
}




