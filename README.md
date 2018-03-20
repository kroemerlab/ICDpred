# ICDpred
ICDpred is an R package for predicting capacity of a molecule to induce ICD

## Requirements
The package can be installed on any distribution running the last version of R. Internet connection is required for installation and package use.
To install the package, you need to have [R](https://cran.r-project.org/) and [Java](https://www.java.com/fr/) installed.

## Packages dependencies
Some package dependencies must be fulfilled before installing the package. Open an R session and type:
```R
install.packages(pkgs=c('devtools','RCurl','rcdk'), repos = "http://cloud.r-project.org")
source("http://bioconductor.org/biocLite.R");biocLite(pkgs='ChemmineR', ask=F)
```

## Installation
In R console, type :
```R
devtools::install_github("kroemerlab/ICDpred")
```

##Usage
After installation, you have first to load the package by typing:
 ```R
library(ICDpred)
```
The *ICDscoring* function can now be used using the PubChem ID of your molecule(s):
 ```R
ICDscoring(CID=4212) #return mitoxantrone results
```
