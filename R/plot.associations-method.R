# Jasmin Straube, Queensland Facility of Advanced Bioinformatics
# Part of this script was borrowed from the graphics and stats package.
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
#' Plot of \code{associations} objects
#' 
#' Plot showing the associated trajectories with or without estimated time shift.
#' 
#' @import gplots
#' @import ggplot2
#' @param x an object of class \code{associations}
#' @param data1 an object of class \code{matrix} or \code{data.frame}.
#' @param data2  an object of class \code{matrix} or \code{data.frame}.
#' @param time a vector of class \code{numeric} presenting the measured time points.
#' @param feature1 the reference feature to visualise, either the index or the name.
#' @param feature2 the associated feature to visualise, either the index or the name.
#' @param cutoff for the associated feature. If \code{fdr=TRUE} the false discovery rate (fdr) corrected p-value (default \code{cutoff=0.05}). If \code{fdr=FALSE} the absolute Pearson Correlation cutoff (default \code{cutoff=0.9}).
#' @param fdr (default TRUE) indicating if the false discovery rate of the corrected p-values from the \code{associations} object should be used as cutoff to visualize associated profiles. If FALSE the absolute Peason correlation is used as cutoff. 
#' @param absCor (default FALSE) if \code{fdr=FALSE} you can choose to visualise associations invariant for positive or negative correlation.
#' @param withShift (default FALSE) indicating if the associated feature should be plotted with the time shift.
#' @param \ldots ignored
#' @details The function allows to visualise features with and without realignement (or shift) of the time profiles according to the estimated delays using associateData() function from the dynOmics package. Features to be visualised can be filtered either using FDR corrected p-values or a correlation threshold.
#' @return plot showing the associated data as calculated by associateData()
#' @seealso \code{\link{associateData}}, \code{\link{summary.associations}}
#' @examples 
#' data(SmallExampleMetabTransc)
#' associations <- associateData(Metabolites[,1:2],Transcripts[,1:100])
#' #if you only define feature1 or feature2 if will plot all associations
#' plot(associations,Metabolites,Transcripts,feature1=1,withShift = TRUE)
#' #if you define feature1 and feature2 it will only plot these two profiles
#' plot(associations,Metabolites,Transcripts,feature1="Metabolite 1",feature2="Transcript 2")
#' @method plot associations
#' @export
plot.associations <- function(x,data1,data2, time, feature1,feature2, cutoff,fdr=T,absCor=T,withShift=F,...){
  l.n <- ifelse(is.null(ncol(data1)),1,ncol(data1))
  l.r <- ifelse(is.null(nrow(data1)),length(data1),nrow(data1))
  if(missing(time)){
    
    time <- 1:l.r
  }
  
  if(missing(feature2) & missing(feature1))
    stop("You need to define feature1 or feature2")
  
  if(!missing(feature1))
  feature1 <- ifelse(is.numeric(feature1),x$Feature1[feature1],feature1)
  if(!missing(feature2))
  feature2 <- ifelse(is.numeric(feature2),x$Feature2[feature2],feature2)
  
  if(missing(cutoff) & fdr)
    cutoff <- 0.05
  
  if(missing(data2))
    data2 <- data1
  
  if(missing(cutoff) & !fdr)
    cutoff <- 0.9
  
  if(min(cutoff)< -1|max(cutoff)>1)
    stop('Select the cutoff within range -1 and 1. ')

  if(!missing(feature2) & !missing(feature1)){

    index <- which(x$Feature1==feature1 & x$Feature2==feature2)
    if(length(index)==0)
      stop('No association found.')
    
    indexT <- x$Feature2[index]

    data <- data2[,indexT]
    dc <- length(indexT)
    

    if(l.n==1){
      data1 <- scale(unlist(as.vector(data1)))
    }else{
      data1 <- scale(unlist(as.vector(data1[,feature1])))
    }
    
    newdf2 <- data.frame(value=data1,time=time,Feature=paste('Feature',feature1))
  }
  
  if(missing(feature2)){
    if(fdr){
      index <- which(x$Feature1==feature1)[p.adjust(x$pAfter[x$Feature1==feature1],method='BH')<cutoff]
    }else{
      if(absCor){
        index <- which(x$Feature1==feature1 & abs(x$corAfter)>=cutoff)
      }else{
        if(cutoff <0){
          index <- which(x$Feature1==feature1 & x$corAfter<=cutoff)
        }else{
          index <- which(x$Feature1==feature1 & x$corAfter>=cutoff)
        }
      }
      
    }
    if(length(index)==0)
      stop('No association found.')

    indexT <- x$Feature2[index]
    
    data <- data2[,indexT]
    dc <- length(indexT)
    
    print(l.n)
    
    if(l.n==1){
      data1 <- scale(unlist(as.vector(data1)))
    }else{
      data1 <- scale(unlist(as.vector(data1[,feature1])))
    }
    
    newdf2 <- data.frame(value=data1,time=time,Feature=paste('Feature',feature1))
  }
  
  if(missing(feature1)){
    if(fdr){
      index <- which(x$Feature2==feature2)[p.adjust(x$pAfter[x$Feature2==feature2],method='BH')<cutoff]
    }else{
      if(absCor){
        index <- which(x$Feature2==feature2 & abs(x$corAfter)>=cutoff)
      }else{
        if(cutoff <0){
      index <- which(x$Feature2==feature2 & x$corAfter<=cutoff)
        }else{
          index <- which(x$Feature2==feature2 & x$corAfter>=cutoff)
    }
      }
    }

    if(length(index)==0)
      stop('No association found.')
    
    indexT <- x$Feature1[index]
    print(l.n)
    
    if(l.n==1){
      data1 <- scale(unlist(as.vector(data1)))
    }else{
      data <- scale(unlist(as.vector(data1[,indexT])))
    }

    dc <- length(indexT)
    newdf2 <- data.frame(value=scale(unlist(as.vector(data2[,feature2]))),time=time,Feature=paste('Feature',feature2))
  }

  
  numt <- length(indexT)
  if(withShift){
    li <- length(time)
 
    index <- na.omit(index)
    delays <- x$delay[index]
    
    wholeL <- li*numt - sum(abs(delays))
    datad <- numeric(wholeL)
    timed <- numeric(wholeL)
    cors <- numeric(wholeL)
    inds <- numeric(wholeL)
    corsi <-  x$corAfter[index]

   for(i in 1:numt){
    times <- NA
    datan <- NA
    if(numt>1){
    dat <- scale(data[,i])
    }else{
    dat <- scale(data)  
    }
     if(delays[i]==0){
       datan <- dat
       times <- time
     }else if(delays[i]>0){
         first <- (delays[i])+1
         datan <- dat[1:(li-delays[i])]
         times <- time[first:li] 
       
       }else{
         first <- abs(delays[i])+1
         datan <- dat[first:li]
         times <- time[1:(li+delays[i])] 
       }
      ti <- length(times)
     range <- ((i-1)*ti+1):(i*ti)
     datad[range] <-datan
     timed[range] <-times
     cors[range] <- rep(corsi[i],each=ti)
     inds[range] <- rep(i, each=ti)
     
   }


    newdf <- data.frame(value=unlist(as.vector(datad)),time=timed,ind=inds,cor=cors)
    
  }else{
  if(numt>1){
    data <- apply(data,2,scale)
  }else{
    data <- scale(data)
  }
  newdf <- data.frame(value=unlist(as.vector(data)),time=rep(time,dc),ind=rep(1:dc,each=length(time)),cor=rep(x$corAfter[index],each=length(time)))
  }
  fe <- ifelse(missing(feature1),feature2,feature1)
  value <- ind <- Feature <- NULL
  l.pos <- length(x$corAfter[index][x$corAfter[index]>0])
  cor.neg <-signif(mean(x$corAfter[index][x$corAfter[index]<0],na.rm=T),2)
  l.neg <- length(x$corAfter[index][x$corAfter[index]<0])
  cor.pos <-signif(mean(x$corAfter[index][x$corAfter[index]>0],na.rm=T),2)
  title <- paste("Associations with feature",fe,ifelse(withShift,"with shift",""),"\n Av cor pos=",cor.pos,'(#',l.pos,"),neg=" ,cor.neg ,'(#',l.neg,')')
  
   ggplot(newdf) + geom_line(aes(x = time, y = value,  group = ind,col=cor),data = newdf) +ggtitle(title) +labs(x='Time ',y='Intensity')+theme(axis.title.x = element_text(size=14),axis.text.x = element_text(size=10),axis.title.y=element_text(size=14),axis.text.y=element_text(size=14)) + scale_color_continuous(limits=c(-1,1),low = "blue",high = "orange") + theme(panel.background = element_rect(fill = 'white', colour = 'white'),legend.title=element_blank())  + geom_line(aes(x = time, y = value,linetype=Feature),size=1.1,data = newdf2)+xlim(range(newdf2$time))#+scale_linetype_identity()
  }
  
