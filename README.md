<img src="./man/figures/miss-SNF_logo.svg" alt="miss-SNF logo" width="275"/>

## Installation

The "miss-SNF" package can be installed using devtools. Please
follow these steps:

0. Download miss-SNF package from this repository.

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
BiocManager::install("RBGL")
BiocManager::install("graph")
BiocManager::install("limma")
install("./missSNF");
```

4. Now miss-SNF can be used as every R package:

```
# Load package
library(missSNF);

# View documentation
?miss.snf()
```
