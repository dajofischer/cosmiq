\name{combine_spectra}

\alias{combine_spectra}

\title{Combine mass spectra of each scan and each file into a single master 
spectrum}

\description{

     \code{combine_spectra} imports each raw file using \code{xcmsRaw} and 
assigns each ion to a previously defined mass vector, which is created using 
bin size parameter \code{mzbin}. This process is repeated for each raw file.
}

\usage{

    combine_spectra(xs, mzbin=0.003, linear=FALSE, continuum=FALSE)    

}

\arguments{

  \item{xs}{xcmsSet object}

  \item{mzbin}{Bin size for the generation of mass vector}

  \item{linear}{logical. If TRUE, linear vector will be generated with 
  \code{mzbin} increments. If FALSE, mass vector will be generated 
  using a non-linear function. This option is recommended for TOF-type 
  mass detectors}

   \item{continuum}{boolean flag. default value is FALSE.}

}

\details{
    This processing step calculates a combined mass spectrum. 
    Mass spectra not only from all scans of a single LCMS run 
    alone are combined but from all acquired datasets. 
    As a result, signal to noise ratio increases for each 
    additional LCMS run. 
}

\author{

   David Fischer 2013
   
}


\examples{

cdfpath <- file.path(find.package("faahKO"), "cdf")

my.input.files <- dir(c(paste(cdfpath, "WT", sep='/'),
    paste(cdfpath, "KO", sep='/')), full.names=TRUE)


# create xcmsSet object
xs <- new("xcmsSet")
xs@filepaths <- my.input.files

x<-combine_spectra(xs=xs, mzbin=0.25,
    linear=TRUE, continuum=FALSE)

plot(x$mz, x$intensity, type='l',
    xlab='m/Z', ylab='ion intensity')


}
