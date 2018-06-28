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
#' Summary of a \code{associations} Object
#' 
#' Summarises the \code{associations} object returned by the \code{\link{associateData}} method. 
#' @param object An object of class \code{associations} .
#' @param ... Additional arguments which are passed to \code{summary}.
#' @return summary of the \code{associations} object.
#' @examples
#' \dontrun{
#' data(Metabolites)
#' data(Transcripts)
#' associations <- associateData(Metabolites[,1],Transcripts[,c(1:50)])
#' summary(associations)
#' }
#' @method summary associations
#' @export
summary.associations <-function(object, ...){
            cat('Identifying associations between time profiles using fast Fourier transform. \n ')
            cat(' \n ')
            cat('Significant associations <0.05 found with shift and without shift:\n')
            print(xtabs(~(object$pAfter<0.05)+(object$pBefore<0.05)))
            cat(' \n ')
            cat('Predicted delays: \n ')
            print(table(object$delay))
          }

