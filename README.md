# DynOmics #

The package provides functions to identify associations within one or between two time course 'omics' data and visualise the associations. Use : associateData to estimate the delays and identify associations of data sets containing time course 'omics' experiments; plot.associations: to visualise associated profiles.


### How to run dynOmics ###
In order to run dynOmics you require the latest version of R (>=3.2) and devtools. You can download the latest R version form [here](https://cran.r-project.org/bin/windows/base/) and follow the instructions to install R. Once installed type the following commands in your R console:

~~~~
install.packages('devtools')
library(devtools)
install_bitbucket('Jasmin87/dynOmics')
library(dynOmics)
~~~~

To see how you can use dynOmics and run examples type in your console:
~~~~
?dynOmics
?associateData
~~~~
### Who do I talk to? ###
The examples will hopefully answer most questions about dynOmics. However, additional questions can be directed to jasmin.straube[at]qimrberghofer.edu.au. We appreciate bug reports in the software or R functions and welcome any suggestions or comments for improvements.