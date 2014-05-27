\name{eicmatrix}

\alias{eicmatrix}

\title{Generate matrix of combined extraced ion chromatograms (EICs)}

\description{

     for each selected mass window, \code{eicmatrix} calculates EICs of
     every raw file and combines them together. It is recommended to
     correct the retention time first using \code{retcor.obiwarp}
}

\usage{

    eicmatrix(xs, xy, center)

}

\arguments{


  \item{xs}{xcmsSet object}

  \item{xy}{table including mz location parameters}

  \item{center}{file number which is used as a template for retention
  time correction}

}

\details{
    
    For each detected mass, an extracted ion chromatogram (EIC) is 
calculated. In order to determine the elution time for each detected mass, 
the EICs of every mass are combined between all acquired runs.

    Make sure that xcms-retention time correction (centwave) was 
applied to the dataset.  The output will be a matrix of EIC's

}

\author{

   David Fischer 2013
   
}


\examples{
      \dontrun{
    #  see package vignette section 
    # 'Generation and combination of extracted ion chromatograms'
    xs<-create_datamatrix(xs=xs, rxy=rxy)
      }
}
