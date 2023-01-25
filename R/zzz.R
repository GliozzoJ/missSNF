missSNFStartupMessage <- function()
{
    # Startup message obtained as 
    # > figlet -f digital "miss-SNF"
    msg <- c(paste0(
        "  
+-+-+-+-+-+-+-+-+
|m|i|s|s|-|S|N|F|
+-+-+-+-+-+-+-+-+ \t version ", packageVersion("missSNF"), 
        "\n\nSimilarity Network Fusion (SNF) algorithm with missing data"))
    return(msg)
}

.onAttach <- function(lib, pkg)
{
    msg <- missSNFStartupMessage()
    if(!interactive())
        msg[1] <- paste("Package 'missSNF' version", packageVersion("missSNF"))
    packageStartupMessage(msg)      
    invisible()
}