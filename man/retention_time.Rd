\name{retention_time}

\alias{retention_time}

\title{Detection of chromatographic peak locations from extracted 
ion chromatograms (EICs)}

\description{

     For each EIC in the EIC matrix, \code{retention_time} localizes 
chromatographic peaks using \code{peakdetection}.
}

\usage{

    retention_time(xy, xs, eicmatrix, RTscales, RTSNR.Th, 
RTSNR.area, mintr)

}

\arguments{

  \item{xy}{table including mz location parameters}

  \item{xs}{xcmsSet object}

  \item{eicmatrix}{Matrix containing extracted ion chromatograms (EICs)}

  \item{RTscales}{Peak Detection of Chromatographic peaks on the 
combined extracted ion chromatogram: scales for continuous wavelet 
transformation (CWT), see also the \code{peakDetectionCWT} function 
of the \code{MassSpecWavelet} package.}
  
  \item{RTSNR.Th}{Peak Detection of Chromatographic peaks on the 
combined extracted ion chromatogram: Signal to noise ratio 
threshold}
  
  \item{RTSNR.area}{Peak Detection of Chromatographic peaks on the 
combined extracted ion chromatogram: Area around the peak for noise 
determination. Indicates number of surrounding peaks on the first 
CWT scale. \code{default = 20}}

  \item{mintr}{Minimal peak width intensity treshold (in 
percentage), for which two overlapping peaks are considered as 
separated (\code{default = 0.5})}

}

\details{Based on the combined extracted ion chromatograms, there is 
another peak detection step to be performed. The same algorithm as 
described for the peak picking of mass signals (continuous wavelet 
transformation) is also used for peak picking in the retention time 
domain.

}

\author{

   David Fischer 2013
   
}
\examples{
      \dontrun{
    #  see package vignette
      }
}

