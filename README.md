
## Overview

An implementation of the Parametric Empirical Bayes (PEB) algorithm for R, 
with methods for estimating the precision of the algorithm and for 
determining sample size requirements when estimating the required parameters.


## Getting Started

1. Install the release version of `devtools` from CRAN with `install.packages("devtools")`.

2. To install mfactor from github use: .


    ```R
    devtools::install_github("jdthorpe/PEB")
    ```
    note that in windows, you may have to build the github tools using: 

    ```R
    library(devtools)
    build_github_devtools()
    #### Restart R before continuing ####
    ```

    * Alternatively, downoload the github repository to a zip file and use:

        ```R
        install.packages("mfactor.zip", repos = NULL)

        # Remove the package after installation
        unlink("mfactor.zip")
        ```

3. Browse the vignettes:
	- [The Quick Start Guide](http://htmlpreview.github.io/?https://github.com/jdthorpe/PEB/blob/master/inst/doc/PEBquickStartGuide.html): `vignette("PEBquickStartGuide")`
	- [A Practical Guide to the PEB](http://htmlpreview.github.io/?https://github.com/jdthorpe/mfactor/blob/master/inst/doc/An-Introduction-to-Mfactors.html): `vignette("peb-parameters")`
	- [Theory Behind the Iterative Method for PEB Parameter Esatimation](http://htmlpreview.github.io/?https://github.com/jdthorpe/mfactor/blob/master/inst/doc/An-Introduction-to-Mfactors.pdf): `vignette("peb-iterative-method")`

