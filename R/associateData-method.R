# Jasmin Straube, Queensland Facility of Advanced Bioinformatics
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU Moleculesral Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Moleculesral Public License for more details.
#
# You should have received a copy of the GNU Moleculesral Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#' Identify associations of trajectories within a data set or across two data sets
#' 
#' Function to estimate differences in expression initation of trajectories to identify associations between time course 'omics' data.
#' 
#' @importFrom snow makeCluster
#' @importFrom parallel stopCluster parLapply detectCores clusterExport
#' @importFrom stats cor cor.test fft na.omit p.adjust xtabs
#' @import methods 
#' @usage associateData(data1,data2,numCores)
#' @param data1 \code{data.frame} or \code{matrix} containing the time as rows and features as columns
#' @param data2 optional an additional \code{data.frame} or \code{matrix} containing the time as rows and features as columns
#' @param numCores alternative \code{numeric} value indicating the number of CPU cores to be used for parallelization. Default value is automatically estimated.
#' @details
#' associateData() takes as input two data sets of interest and performs a pairwise associations comparison between features using a fast Fourier transform approach to detect delays (also called 'associations') between the different features. Note that the argument `numCores` indicates the number of CPUs  and is detected by default in the function to perform parallelization. The final result is a table with a row for each pairwise comparison. The output presents the dynOmics estimated delay between two features, the p-value (`p`) and correlation coefficient (`cor`) from a Pearson's test, before and after the time profiles have been realigned according to the dynOmics estimated delay.
#' @return associateData returns an object of class \code{associations} containing the following components:
#' \itemize{
#' \item{Feature1}{ \code{character} the colnames or the index of data1.}
#' \item{Feature2}{ \code{character} the colnames or the index of data2.}
#' \item{delay}{ \code{numeric} estimated delay between feature1 and feature2.} 
#' \item{pBefore}{ \code{numeric} p-value of the test for association before applying the predicted time shift.}
#' \item{pAfter}{ \code{numeric} p-value of the test for association after applying the predicted time shift.}
#' \item{corBefore}{ \code{numeric} Pearson correlation before applying the predicted time shift.} 
#' \item{corAfter}{ \code{numeric} Pearson correlation after applying the predicted time shift.}
#' }
#' @references  Straube J., Bernard A., Huang B.E., Le Cao K.-A.(2015).  \emph{DynOmics - A new algorithm using fast Fourier transform to reveal dynamic molecule interactions } In preparation
#' @seealso \code{\link{summary.associations}}, \code{\link{plot.associations}}
#' @examples 
#' data(SmallExampleMetabTransc)
#' associations <- associateData(Metabolites[,1],Transcripts[,1:50])
#' summary(associations)
#' plot(associations,Metabolites,Transcripts,feature1=1)
#' @docType methods
#' @rdname associateData-methods

#setGeneric('associateData',function(data1,data2,numCores){standardGeneric('associateData')})
#setClassUnion("missingOrnumeric", c("missing", "numeric"))
#setClassUnion("matrixOrframe",c('matrix','data.frame'))
#setClassUnion("matrixOrframeorMissing",c('matrix','data.frame','missing'))
### @rdname associateData-methods
## @aliases associateData,matrixOrframe,matrixOrframeorMissing,missingOrnumeric-method
## @exportMethod associateData

#setMethod('associateData',c(data1="matrixOrframe",data2="matrixOrframeorMissing",numCores="missingOrnumeric"),function(data1,data2,numCores){
#  associateData(data1=data1,data2=data2,numCores=numCores)
#})

#' @export
associateData <- function(data1,data2,numCores){
  if(missing(numCores)){
    num.Cores <- detectCores()
  }else{
    num.Cores <- detectCores()
    if(num.Cores<numCores){
      warning(paste('The number of cores is bigger than the number of detected cores. Using the number of detected cores',num.Cores,'instead.'))
    }else{
      num.Cores <- numCores
    }
  }
  
  singleData <- 0
  if(missing(data2)){
    data2 <- data1
    singleData <- 1
  }
    
  
  if(sum(is.na(data1))>0|sum(is.na(data2))>0)
    stop("No missing data allowed.")
    
  
  nc <- ifelse(is.null(ncol(data1)),1,ncol(data1))
  nc2 <- ncol(data2)
  
  l.n <- ifelse(is.null(nrow(data1)),length(data1),nrow(data1))
  if(l.n!=nrow(data2))
    stop('Data must have the same time dimension.')
  


  #m <- matrix(NA,ncol=7,nrow=nc*nc2F)

  cl <- makeCluster(num.Cores,"SOCK")
  clusterExport(cl, list('data1','data2','fft','get.delay','convert.fft','nc','nc2','singleData'),envir=environment())

  ncF <- ifelse(singleData,nc-1,nc)
  df <- parLapply(cl, 1:ncF,fun = function(i){
# for(i in 1:nc){
    index <- 0
    if(singleData){
    comp <- c((i+1):nc2)
    }else{
      comp <- c(1:nc2)
    }
    m <- matrix(NA,ncol=7,nrow=length(comp))
    for(j in comp){
      index <- index+1
      if(nc==1){
        new.data1 <- scale(data1)
      }else{
        new.data1 <- scale(data1[,i])
      }
      
      new.data2 <- scale(data2[,j])
      
      cori<-pori <-delay2 <- pfft <- pori <-NA
      
      if(sum(is.na(new.data1))>0|sum(is.na(new.data2))>0)
        next

      del <- get.delay(new.data1,new.data2)

      m[index,] <- c(i,j,del)
    }
    return(m)
  })
  
  
  stopCluster(cl)
  sdf <- as.data.frame(do.call(rbind, lapply(df, unlist)))
  colnames(sdf) <- c('Feature1','Feature2','delay','pBefore','pAfter','corBefore','corAfter')
  if(!is.null(colnames(data1)))
    sdf$Feature1 <- colnames(data1)[sdf$Feature1]
  if(!is.null(colnames(data2)))
    sdf$Feature2 <- colnames(data2)[sdf$Feature2]
  class(sdf) <- c('associations','data.frame')
  return(sdf)
}


convert.fft <- function(cs, sample.rate=1) {
  cs <- cs / length(cs) # normalize
  distance.center <- Mod(cs)
  angle           <- (180*Arg(cs)/pi)
  data.frame(strength = distance.center,delay=angle)
}

get.delay<- function(x,y) {
  cord <- cor.test(x,y)
  if(abs(cord$estimate)!=1 & !is.na(cord)){
  l <- length(x)
  ori.fft<- convert.fft(fft(x))[2:(l/2),]

  cf <-convert.fft(fft(y))[2:(l/2),]
  #frequencey at maximum amplitude cycle
  freq <- which.max(ori.fft$strength)
  dels <-  ori.fft$delay[freq]- cf$delay[freq]
 
  dels <- ifelse(dels<0,dels+359,dels)
  #sequence length
  l2 <- l
  if(l%%2!=0)
    l2<-l+1
  l.freq <- l2/freq
  
  newest <- round(dels/(360/l.freq))
  if(dels>-1 & dels<90){
    newest <-  newest
  }
  #positive correlated with neg shift
  if(dels>269 & dels<360){
    newest <-  newest-l.freq
  }
  #negative correlated with neg shift
  if(dels>89 & dels<180){
    newest <-  (l.freq/2)-newest
  }
  #negative correlated with pos shift
  if(dels>180 & dels<270){
    newest <-  (l.freq/2)-newest
  }

  newest <- ifelse(dels=='180', 0,newest)
  newest <- round(newest)
  newest <- ifelse(abs(newest)==(l.freq),0,newest)
  newest <- ifelse((newest)>(l.freq/2),(newest)-(l.freq),newest)
  newest <- ifelse((newest)<(-l.freq/2),(newest)+(l.freq),newest)
  newest <- round(newest)
  delaytests <- c(newest, newest*-1)
  
  cors <- numeric(2)
  for(d in 1:2){
    dely <- delaytests[d]
    cors[d] <- ifelse(dely>0,cor(as.numeric(y[1:(l-abs(dely))]),as.numeric(x[(abs(dely)+1):l])),cor(as.numeric(x[1:(l-abs(dely))]),as.numeric(y[(abs(dely)+1):l])))
  }
  s <- sum(!is.na(cors))
  if(s==2){
    newest <- ifelse(abs(cors[1])>abs(cors[2]),newest,newest*-1)
  }else if(s==1){
    newest <- c(newest,newest*-1)[which(!is.na(cors))]
  }else{
    newest <-0
  }
  delaytests <- c(newest-1,newest,newest+1)
  delaytests[abs(delaytests)>=(l.freq) | is.na(delaytests)] <- 0
  delaytests <- unique(delaytests)
  cors <-  corsp <- numeric(length(delaytests))
  for(d in 1:length(delaytests)){
    dely <- delaytests[d]
    if(dely>0){
      c <- cor.test(as.numeric(y[1:(l-abs(dely))]),as.numeric(x[(abs(dely)+1):l]))
    }else{
      c <- cor.test(as.numeric(x[1:(l-abs(dely))]),as.numeric(y[(abs(dely)+1):l]))
    }
    cors[d] <-c$estimate
    corsp[d] <- c$p.value

  }
  ma <- which.max(abs(cors))
  cori <- cor.test(x,y,use="pairwise.complete.obs")
  delay2  <- delaytests[ma]
  pfft <- corsp[ma]
  cfft <- cors[ma]
  
  if(abs(cfft)<abs(cori$estimate)){
    delay2 <- 0
    cfft <- cori$estimate
    pfft <- cori$p.value
  }
  return(c(delay2,cori$p.value,pfft, cori$estimate,cfft))
  }else{
    
  }

  return(c(0,cord$p.value,cord$p.value, cord$estimate,cord$estimate))
}


