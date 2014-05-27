\name{peakdetection}

\alias{peakdetection}

\title{An algorithm for the detection of peak locations and boundaries 
      in mass spectra or ion chromatograms}

\description{

     \code{peakdetection} uses Continuous wavelet transformation 
(CWT) to determine optimal peak location. A modified algorithm of Du 
et al. (2006) is used to localize peak positions.}

\usage{

    peakdetection(scales, y, x, SNR.Th, SNR.area, mintr)    

}

\arguments{

  \item{scales}{vector with the scales to perform CWT}

  \item{y}{vector of ion intensities}

  \item{x}{vector of mz bins}

  \item{SNR.Th}{Signal-to-noise threshold}
  
  \item{SNR.area}{Window size for noise estimation}

  \item{mintr}{Minimal peak width intensity treshold (in 
percentage), for which two overlapping peaks are considered as 
separated. default is set to 0.5.}

}

\details{

 A peak detection algorithm based on continuous wavelet 
transformation (CWT) is used for this step (modified from Du et al., 
2006). Peak detection based on CWT has the advantage that a sliding 
scale of wavelets instead of a single filter function with fixed 
wavelength is used. This allows for a flexible and automatic 
approximation of the peak width. As a result it is possible to 
locate both narrow and broad peaks within a given dynamic range.

}

\author{

   David Fischer 2013
   
}


\references{

    Du, P., Kibbe, W. A., & Lin, S. M. (2006). 
    Improved peak detection in mass spectrum by incorporating 
continuous wavelet transform-based pattern matching. Bioinformatics, 
22(17), 2059-2065. doi:10.1093/bioinformatics/btl355

}

\examples{
    
cdfpath <- file.path(find.package("faahKO"), "cdf")

my.input.files <- dir(c(paste(cdfpath, "WT", sep='/'),
    paste(cdfpath, "KO", sep='/')), full.names=TRUE)


# create xcmsSet object
xs <- new("xcmsSet")
xs@filepaths <- my.input.files

op<-par(mfrow=c(3,1))

x<-combine_spectra(xs=xs, mzbin=0.25,
    linear=TRUE, continuum=FALSE)

plot(x$mz, x$intensity, 
    type='h', 
    main='original',
    xlab='m/Z', ylab='ion intensity')

xy <- peakdetection(x=x$mz, y=x$intensity,
    scales=1:10, SNR.Th=1.0, SNR.area=20, mintr=0.5)

id.peakcenter<-xy[,4]

plot(x$mz[id.peakcenter], x$intensity[id.peakcenter], type='h', 
    main='filtered')


plot(x$mz, x$intensity, type='l',
    xlim=c(400,450),
    main='zoom',
    log='y',
    xlab='m/Z', ylab='ion intensity (log scale)')

points(x$mz[id.peakcenter], x$intensity[id.peakcenter],col='red', type='h')
    
}
