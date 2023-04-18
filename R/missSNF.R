#' miss-SNF: integration of patients' networks considering missing cases
#'
#' @description Extension of the algorithm Similarity Network Fusion
#' (\url{https://doi.org/10.1038/nmeth.2810}) able to handle the absence (complete or
#' nearly complete) of a specific data source for a given patient by:
#' \loadmathjax
#' \itemize{
#' \item "reconstruct strategy" that can partially reconstruct missing
#'   data by using information from different sources (i.e. miss-SNF ONE).
#' \item "ignore strategy" which simply ignores missing data during the integration
#'   process (i.e. miss-SNF ZERO).
#' \item "equidistant strategy", similar to reconstruct but sets the similarity
#' of the partial samples with the others to a fixed value instead of zero.
#' \item "random strategy" that sets the similatiry randomly.
#' }
#'
#' @param Mall list of named matrices/dataframes (samples x features).
#' @param sims vector of strings. It is a vector containing the names of the
#'             similarity measures to apply to the matrices in Mall.
#'             "scaled.exp.euclidean" is the scaled exponential euclidean distance;
#'             "scaled.exp.chi2" is the scaled exponential chi-square distance
#' @param sims.arg list. List with the same length of "sims" where each elements
#'                 is a list containing additional arguments for each
#'                 similarity measure in the argument "sims". Set element to
#'                 NULL if you want to use default parameters for a specific
#'                 similarity measure (default).
#' @param mode string. If you want to partially reconstruct missing data use
#'             "reconstruct" or "one" (default), otherwise ignore them during integration
#'             using "ignore" or "zero". If you use "equidistant", the self-loop
#'             for partial samples is set to 0.5 while the similarity with other
#'             samples is set to the same value so that the sample is equidistant
#'             from all other samples but more similar to itself than others.
#'             The option "random" sets randomly the similarity of partial
#'             samples (let's consider this as a baseline).
#' @param perc.na percentage of NAs above which a patient is considered missing.
#' @param miss.symbols vector of strings. If not NULL, the provided symbols
#'                     in matrices are converted to NA.
#' @param K Number of neighbors in K-nearest neighbors part of the algorithm.
#' @param t Number of iterations for the diffusion process.
#' @param impute string. Kind of imputation method to apply in case of samples
#'               with few missing values. Options are: NULL, "mean", "median".
#'               NOTE: if impute=NULL, then perc.na has to be 0 (i.e.
#'               having even one NA will make the patient treated as missing).
#' @param d numeric. Set the diagonal of the matrix to "d" if
#' mode="reconstruct" or mode="one" (def d=1).
#' @param random.walk string. Use 1-step Random Walk to compute the local
#' similarity matrix S and/or p-step Random Walk (p>=2) to compute the global
#' similarity matrix P. random.walk=c("global", "local", "both", "none") and
#' defaults is random.walk="none".
#' @param p numeric. Number of steps for the p-step RW. Used only when
#' global similarity matrix is computed through p-step Random Walk.
#' @param seed numeric. Seed to get reproducible results. Needed only if
#' mode = "random".
#'
#' @return A list with two elements:
#' \itemize{
#' \item W : integrated similarity matrix. Note that the order of the patients
#' is different from the order of the original matrices.
#' \item removed.pts : vector with names of removed patients (i.e. patients
#' present only in one matrix of Mall and having too much NAs or, more in general,
#' if a patient is considered missing in all data sources).
#' }
#' @export
#' @importFrom stats runif
#'
#' @examples
#'
#' # Create list of input matrices
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
#' M3[7, ] <- rep(NA, ncol(M3));
#'
#' Mall <- list("M1"=M1, "M2"=M2, "M3"=M3);
#'
#' # Call miss.snf using "reconstruct strategy"
#' W.r <- miss.snf(Mall, sims=rep("scaled.exp.euclidean", 3),
#'                sims.arg=list(list(kk=2), list(kk=2), list(kk=2)),
#'                mode="reconstruct", K=3);
#'
#' # Call miss.snf using "ignore strategy"
#' W.i <- miss.snf(Mall, sims=rep("scaled.exp.euclidean", 3),
#'                sims.arg=list(list(kk=2), list(kk=2), list(kk=2)),
#'                mode="ignore", K=3);
miss.snf <- function(Mall, sims, sims.arg=vector("list", length(sims)),
                     mode="reconstruct", perc.na=0.2,
                     miss.symbols=NULL, K=20, t=20, impute="median",
                     d=1, random.walk="none", p=3, seed=NULL) {

    # Allow naming of mode as "one" or "zero"
    # Possible values for mode: one/reconstruct or zero/ignore
    if(mode == "one"){
        mode <- "reconstruct"
    }

    if(mode == "zero"){
        mode <- "ignore"
    }

    if(!(mode == "reconstruct" | mode == "ignore" | mode == "random" | mode == "equidistant")){

        stop("miss.snf: mode can be only one/reconstruct, zero/ignore, equidistant or random.")
    }

    # Check that a similarity measures is provided for each matrix
    if(length(Mall) != length(sims)){

        stop("miss.snf: number of matrices to integrate has to be equal to
             number of similarity measures provided.");
    }

    # Convert miss.symbols to NA
    if(!is.null(miss.symbols)){
        for(i in 1:length(Mall)){
            for(j in 1:length(miss.symbols)){
                Mall[[i]][Mall[[i]] == miss.symbols[j]] <- NA;
            }
        }
    }

    # Get list of missing patients for each matrix
    miss.pts <- get.miss.pts(Mall, perc.na=perc.na, miss.symbols=NULL);

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
    }else if(impute == "median" | impute == "mean"){
        for(i in 1:length(Mall)){
            Mall[[i]] <- impute_miss(Mall[[i]], perc.na=perc.na, method=impute,
                                     verbose=FALSE)
        }
    } else {
        stop("miss.snf: imputation method can be only median or mean.")
    }

    # Compute similarity matrices W (without considering missing patients)
    Wall <- list();
    for(i in 1:length(Mall)){

        idx.miss.pts <- which(rownames(Mall[[i]]) %in% miss.pts[[i]]); #idx patients with NAs
        sim.fun <- get(sims[[i]]);

        # Handle presence of matrices without patients to remove for NAs
        if(length(idx.miss.pts) != 0){
            dc <- Mall[[i]][-idx.miss.pts, ]
            #Wall[[i]] <- sim.fun(Mall[[i]][-idx.miss.pts, ], ...);
            Wall[[i]] <- do.call(sim.fun, c(list(dc), sims.arg[[i]]))
        } else{
            dc <- Mall[[i]]
            #Wall[[i]] <- sim.fun(Mall[[i]], ...);
            Wall[[i]] <- do.call(sim.fun, c(list(dc), sims.arg[[i]]))
        }
    }

    # NOTE: Here Wall represent a list of the similarity function (either the scaled
    # exponential similarity kernel or the chi square similarity kernel depending
    # on the arguments of the parameter sims). Matrices are not of the same size.


    # Align networks replacing missing patients with 0
    Wall_aligned <- NetInt::lalign.networks(fill=0, Wall);

    # Get names of patients removed from matrices because present in only
    # one matrix but having too much NAs
    removed.pts.tab <- table(unlist(miss.pts)) == length(Mall);
    removed.pts <- rownames(removed.pts.tab)[removed.pts.tab];

    # Set diagonal of missing elements to "d" if "reconstruct" mode is used
    # If "ignore" mode, missing patients are already set to vector of zeros.
    # If mode is "equidistant", the diagonal is set to 0.5 and the other elements
    # to 0.5/(|W|-1) for partial samples. If mode is "random", similarity for
    # partial samples is set to random values from uniform distribution.
    idx.miss.aligned <- list();
    if(mode == "random"){set.seed(seed)} #set seed for reproducibility if random

    for(i in 1:length(Wall_aligned)){

        idx <- which(!(rownames(Wall_aligned[[i]]) %in% rownames(Wall[[i]])));
        idx.miss.aligned[[i]] <- idx;
        if(mode == "reconstruct"){

            diag(Wall_aligned[[i]])[idx] <- d;

        } else if(mode == "equidistant"){

            eq <- 0.5/(ncol(Wall_aligned[[i]]) - 1)
            Wall_aligned[[i]][idx, ] <- eq
            Wall_aligned[[i]][, idx] <- eq

            diag(Wall_aligned[[i]])[idx] <- 0.5

        } else if(mode == "random"){

            # Create matrix of 0s with same dimensions of Wall_aligned[[i]]
            ran <- matrix(0, nrow=nrow(Wall_aligned[[i]]), ncol=ncol(Wall_aligned[[i]]))
            # Set columns for partial samples to random values from uniform
            # distribution (min and max of the distribution are taken from the
            # matrix). A little epsilon (.Machine$double.eps) is added to avoid
            # zeros among random elements.
            ran[, idx] <- runif(nrow(Wall_aligned[[i]])*length(idx),
                                min = min(Wall_aligned[[i]]),
                                max = max(Wall_aligned[[i]])) + .Machine$double.eps
            ran <- (ran+t(ran))/2 #make symmetric
            Wall_aligned[[i]][ran != 0] <- ran[ran != 0]
        }
    }

    # Use Wall_aligned in SNF code
    # NOTE: Wall_aligned represent the similarity kernel (either exponential or
    # chi square depending on sims) where for the missing patients i
    # (in reconstruct mode) have values Wall_aligned[][i,i] = 1
    Wall <- Wall_aligned;

    # FROM NOW ON STARTS SNF CODE (as implemented in SNFtool + changes
    # introduced on matrices P and S for missSNF)
    ############################################################################
    ############################################################################

    wall.names <- dimnames(Wall[[1]])

    LW <- length(Wall)

    #Normalize different networks to avoid scale problems.
    #### [J]: Wall after normalization contains "global" similarity matrices P
    #### [J]: Global similarity matrices can be computed by Wang's normalization
    #### or by p-step Random Walk Kernel.
    newW <- vector("list", LW)
    nextW <- vector("list", LW)
    for(i in 1:LW){

        if(random.walk == "global" | random.walk == "both"){

            Wall[[i]] <- RANKS::rw.kernel(Wall[[i]])
            Wall[[i]] <- RANKS::p.step.rw.kernel(Wall[[i]], p=p)
            Wall[[i]] <- NetPreProc::Prob.norm(Wall[[i]]) # Normalization with unnormalized lagrangian

        } else if(random.walk == "none" | random.walk == "local"){

            Wall[[i]] <- .normalize(Wall[[i]])

        } else {
            stop("miss.snf: random.walk can be only c('global', 'local',
                 'both', 'none')")
        }

        Wall[[i]] <- (Wall[[i]] + t(Wall[[i]]))/2

        ### [J]: set diagonal to d if mode="reconstruct", 0 otherwise for missing
        ### patients. The other elements for missing patients should remain zero
        ### even after normalization or p-step Random Walk. "equidistant" and
        ### "random" modes are implemented as before.

        if(mode == "equidistant"){

            idx <- idx.miss.aligned[[i]]

            eq <- 0.5/(ncol(Wall[[i]]) - 1)
            Wall[[i]][idx, ] <- eq
            Wall[[i]][, idx] <- eq

            diag(Wall[[i]])[idx] <- 0.5

        } else if(mode == "random"){

            idx <- idx.miss.aligned[[i]]
            ran <- matrix(0, nrow=nrow(Wall[[i]]), ncol=ncol(Wall[[i]]))
            ran[, idx] <- runif(nrow(Wall[[i]])*length(idx),
                                min = min(Wall[[i]]),
                                max = max(Wall[[i]])) + .Machine$double.eps
            ran <- (ran+t(ran))/2
            Wall[[i]][ran != 0] <- ran[ran != 0]


        } else if(mode == "reconstruct"){

            diag(Wall[[i]])[idx.miss.aligned[[i]]] <- d;

        } else {

            diag(Wall[[i]])[idx.miss.aligned[[i]]] <- 0;
        }
    }

    # NOTE: Wall here is the global similarity matrix P with Wall[][i,i]=1 in
    # reconstruct mode for the completely missing patients i

    ### Calculate the local transition matrix. (KNN step?)
    #### [J]: newW contains the local transition matrices S
    #### [J]: since missing patients has no neighbours by design
    #### (except for itself), in case of reconstruction strategy the weights
    #### are already set to zero and the diagonal to 1. In case of ignore
    #### strategy, the weights are already set to zero (division by zero
    #### avoided by the function .local.similarity.matrix()).
    #### [J]: local similarity matrix can be computed also by 1-step RW.
    #### When using 1-step RW, all weights are set to zero except for the
    #### diagonal.
    for(i in 1:LW){

        if(random.walk == "local" | random.walk == "both"){

            newW[[i]] <- RANKS::rw.kernel(Wall_aligned[[i]])
            newW[[i]] <- NetPreProc::Prob.norm(newW[[i]]) # Normalization with unnormalized lagrangian

        } else if (random.walk == "none" | random.walk == "global"){

            newW[[i]] <- .local.similarity.matrix(Wall_aligned[[i]], K);

        } else {
            stop("miss.snf: random.walk can be only c('global', 'local',
                 'both', 'none')")
        }

        # Set diagonal to d for mode reconstruct and 0 for mode ignore.
        # "equidistant" and "random" modes are computed as before.
        if(mode == "equidistant"){

            idx <- idx.miss.aligned[[i]]

            eq <- 0.5/(ncol(newW[[i]]) - 1)
            newW[[i]][idx, ] <- eq
            newW[[i]][, idx] <- eq

            diag(newW[[i]])[idx] <- 0.5

        } else if(mode == "random"){

            idx <- idx.miss.aligned[[i]]
            ran <- matrix(0, nrow=nrow(newW[[i]]), ncol=ncol(newW[[i]]))
            ran[, idx] <- runif(nrow(newW[[i]])*length(idx),
                                min = min(newW[[i]]),
                                max = max(newW[[i]])) + .Machine$double.eps
            ran <- (ran+t(ran))/2
            newW[[i]][ran != 0] <- ran[ran != 0]

        } else if(mode == "reconstruct"){

            diag(newW[[i]])[idx.miss.aligned[[i]]] <- d;

        } else {

            diag(newW[[i]])[idx.miss.aligned[[i]]] <- 0;
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

    ### [J]: if we use miss-SNF ignore, division considers only the
    ### actual number of available data sources
    if(mode != "ignore"){
        W <- W/LW
    } else{
        M <- simplify2array(Wall)
        M[M > 0] <- 1 # [J]: create an indicator array
        M <- rowSums(M, dims = 2) # [J]: count data sources
        W <- W * (1/M)
        W[is.nan(W)] <- 0 # [J]: replace NaN introduced with division by zero,
                          # possible when pts i and j are no present in the
                          # same data source, thus edge(i,j)=0 for each data source
    }


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
#' euclidean distance (euclidean distance without the square root).
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
#' @return list of vectors. Each vector contains the names of missing samples
#'         in a specific matrix/dataframe.
#' @export
#'
#' @examples
#'
#' # Create list of input matrices
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
#' M3[7, ] <- rep(NA, ncol(M3));
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


#' Impute patients with few NAs
#'
#' @description This function imputes missing samples using some specified
#' method (e.g. mean, median) only for patients/samples having a percentage
#' of NAs lower or equal to perc.na. In other words, not all missing data are
#' imputed but only the ones for which a patient has enough data to not be
#' considered completely missing.
#'
#' @param data matrix. Matrix (samples x features).
#' @param perc.na numeric. Percentage of missing features.
#' @param method string. An imputation method (possible options mean, median).
#' @param verbose boolean. Want to print messages? (def. TRUE)
#'
#' @return Imputed data matrix.
#' @export
impute_miss <- function(data, perc.na, method, verbose=TRUE) {

    # Retrieve imputation method
    impute_func <- get(method)

    # Get rows to impute (having few missing patients but not zero)
    n_na <- apply(is.na(data), 1, sum)
    idx_impute <- which(n_na/ncol(data) <= perc.na & n_na/ncol(data) != 0)

    if(length(idx_impute) == 0 & verbose == TRUE){
        message("impute_miss: no samples to impute are present.")
    }

    # Impute each column
    for (p in idx_impute){
        for(i in 1:ncol(data)){
            bool <- is.na(data[p, ])
            data[p, i] <- ifelse(bool[i], impute_func(data[, i], na.rm = TRUE), data[p, i])
        }
    }

    return(data)
}

