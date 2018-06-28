#' Metabolite and Transcript Simulation Data
#'
#' Simulated data were received from Redestig et al., 2011. Metabolite and transcript levels were obtained using an impulse model
#' (Chechik and Koller, 2009). Functions were used to model five different metabolite patterns and for each metabolite 50 associated
#' transcript levels. Time lags were introduced in the range from -2 to 2 with the probability 0.1, 0.2, 0.4, 0.2, 0.1. Simulated
#' profiles have seven time points and normal distributed noise was introduced with mean zero and standard deviation 0.1. 
#'
#' \itemize{
#'   \item Metabolites. data matrix with 7 rows and 5 columns. Each row represents an experimental time sample, and each column a single metabolite.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Metabolites
#' @usage data(Metabolites)
#' @format  This data set contains the simulated expression of 5 metabolites for 7 time points. 
#' @references   Redestig,H. and Costa,I.G. Detection and interpretation of metabolite-transcript coresponses using combined profiling data. \emph{Bioinformatics} \bold{27}(13) (2011), pp. i357 65.
#' @source   The Metabolite Simulation Data is based on the the paper of Redestig \emph{et al.} (2011).
NULL
