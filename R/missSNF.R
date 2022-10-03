# Implementation of missSNF, a modification of the algorithm Similarity Network 
# Fusion (https://doi.org/10.1038/nmeth.2810) to handle the absence 
# (complete or nearly complete) of a specific data modality.
# In particular, two strategies to implement missSNF are considered:
# - miss.snf.reconstruct() can partially reconstruct missing
#   data by using information from different sources
# - miss.snf.ignore() which simply ignores missing data during the integration
#   process




#' Compute similarity measures considering missing samples
#'
#' @param M matrix (samples x features).
#' @param miss.pts vector. Contains the names of missing patients for M.
#' @param sim function. Similarity measure to compute similarity matrix W.
#' @param mode string. It can be "reconstruct" or "ignore", depending on the 
#'             miss-SNF modality that we want to exploit.
#'             
#'
#' @return similarity matrix W.
#' @export
miss.sim <- function(M, miss.pts, sim, mode="reconstruct"){

    
    
    return(W)    
}



#' miss-SNF Reconstruct
#' 
#' @description extension of the algorithm Similarity Network Fusion
#' (https://doi.org/10.1038/nmeth.2810) able to handle the absence (complete or
#' nearly complete) of a specific data source for a given patient by partially
#' reconstructing the missing data using information available from the other
#' data sources.
#'
#' @param 
#'             
#' @param K 
#' @param t 
#'
#' @return 
#' @export
miss.snf.reconstruct <- function() {
    
}


###############################################################################
miss.snf.ignore <- function() {
    
    
}


################################################################################
############## Similarity measures #############################################
################################################################################

#' Scaled exponential euclidean distance
#' 
#' @description it is a wrapper function to directly compute the scaled 
#' exponential euclidean distance as implemented in SNFtool library.
#' # NOTE: from SNFtool "If the data is continuous, we recommend to use 
#' the function "dist2" as follows"
#' dist <- (dist2(as.matrix(Data1),as.matrix(Data1)))^(1/2);
#' NOTE: this is because the function dist2 actually computes the squared 
#'       euclidean distance (without the square root).
#'
#' @param M matrix (samples x features)
#' @param K integer. Number of nearest neighbours to sparsify similarity.
#' @param sigma numeric. Variance for local model.
#'
#' @return similarity matrix computed using scaled exponential euclidean 
#'         distance.
#' @export
scaled.exp.euclidean <- function(M, K=20, sigma=0.5){
    
    dist <- (dist2(M, M))^(1/2);
    sim <- affinityMatrix(dist, K=K, sigma=sigma);
    
    return(sim)
}



#' Scaled exponential chi-square distance
#' 
#' @description it is a wrapper function to directly compute the scaled 
#' exponential chi-square distance as implemented in SNFtool library.
#'
#' @param M matrix (samples x features)
#' @param K integer. Number of nearest neighbours to sparsify similarity.
#' @param sigma numeric. Variance for local model.
#'
#' @return similarity matrix computed using scaled exponential chi-square 
#'         distance.
#' @export
scaled.exp.chi2 <- function(M, K=20, sigma=0.5){
    
    dist <- chiDist2(M);
    sim <- affinityMatrix(dist, K=K, sigma=sigma);
    
    return(sim)
}


###############################################################################
############## UTILITIES ######################################################
###############################################################################

#' Get list of missing patients for each matrix
#' 
#' @description This function takes in input a list of matrices Mall and returns
#' for each matrix the list of patients that are missing (considering the union
#' of patients in all matrices). A patient is considered "missing" if (I) it is 
#' present in the matrix but has |NA| > perc.na, (II) a patient is not present
#' in one matrix but is present in at least one of the others.
#' 
#' @param Mall list of named matrices/dataframes (samples x features).
#' @param perc.na percentage of NAs above which patient is considered missing.
#' @param miss.symbol if not NULL, the provided symbol in matrices is converted 
#'                    to NA.
#'
#' @return list of vectors. Each vector contains the indices of missing samples
#'         in a specific matrix.
#' @export
get.miss.pts <- function(Mall, perc.na=0.2, miss.symbol=NULL){
    
    # Check if matrices are named
    for (i in 1:length(Mall)){
        if (is.null(rownames(Mall[[i]]))){
            stop("get.missing.pts: matrices must be named")
        }
    }
    
    # Convert miss.symbol into NA
    if(!is.null(miss.symbol)){
        for(i in 1:length(Mall)){
            Mall[[i]][Mall[[i]] == miss.symbol] <- NA;
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
        pts.many.na <- rownames(Mall[[i]][(n.na/ncol(Mall[[i]])) >= perc.na, ]);
        miss.pts[[i]] <- c(miss.pts[[i]], pts.many.na);
    }
    
    return(miss.pts)
}




