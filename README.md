<div align="center">
<img src="./man/figures/miss-SNF_logo.svg" alt="miss-SNF logo" width="275"/>
</div>


## Installation

The "miss-SNF" package can be installed using devtools. Please
follow these steps:

0. Download or clone miss-SNF package from this repository.

```
# Code to clone repository
git clone https://github.com/AnacletoLAB/missSNF.git
```

You can download the .zip file containing this repository using the "< > Code" button.

1. Install devtools using:

```
install.packages("devtools");
```

2. Load devtools package:

```
library("devtools");
```

3. Then, you can install miss-SNF (some Bioconductor packages
has to be installed apart):

```
# Install BiocManager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor dependencies
BiocManager::install("RBGL")
BiocManager::install("graph")
BiocManager::install("limma")

# Install miss-SNF
install("./missSNF");
```

4. Now miss-SNF can be used as every R package:

```
# Load package
library(missSNF);

# View documentation
?miss.snf()
```
