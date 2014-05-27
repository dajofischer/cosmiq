\name{quantify_combined}

\alias{quantify_combined}

\title{Generate report with combined ion intensities}

\description{

     \code{quantify_combined} calculates a data table of mz intensities
     from a sum of ions within a selected mz window and from all MS scans.
}

\usage{

    quantify_combined(xs,xy, rtcombine)

}

\arguments{

  \item{xs}{xcmsSet object}

  \item{xy}{table including mz location parameters}

  \item{rtcombine}{Numerical, with two entries. If no retention 
  time information is to be used, \code{rtcombine} provides the 
  retention time window where ion intensities will be summed for 
  quantification. If \code{rtcombine = c(0,0)}, 
  All scans will be used}

}

\details{

 This function is only used, when no retention time information 
 should be used, for example with direct infusion experiments. 

}

\author{

   David Fischer 2013
   
}
\examples{
      \dontrun{
    #  see package vignette
      }
}

