\name{cosmiq}

\alias{cosmiq}

\title{cosmiq - main wrapper function}

\description{
    This is the main wrapper function for the package cosmiq. Every 
processing step of \code{cosmiq} will be calculated during this function, 
including mass spectra combination, detection of relevant masses, generation 
and combination of extracted ion chromatograms, detection of chromatographic 
peaks, Localisation and quantification of detected peaks. 

     \code{combine_spectra} imports each raw file using \code{xcmsRaw} and 
assigns each ion to a previously defined mass vector, which is created using 
bin size parameter \code{mzbin}. This process is repeated for each raw file.
}

\usage{

    cosmiq (files, RTinfo=TRUE, mzbin=0.003, 
        center=0, linear=FALSE, profStep=1, retcorrect=TRUE, 
        continuum=FALSE, rtcombine=c(0,0), scales=c(1:10), 
        SNR.Th=10, SNR.area=20, RTscales=c(1:10, seq(12,32,2)), 
        RTSNR.Th=20, RTSNR.area=20, mintr=0.5)

}

\arguments{

  \item{files}{String vector containing exact location of each file}
  
  \item{RTinfo}{Logical indicating whether retention time should be 
used or not. If FALSE, the quantitative information will be 
retrieved as a sum of ion intensities whithin a retention time 
window. This retention time window }
  
  \item{rtcombine}{Numerical, with two entries. If no retention time 
information is to be used, \code{rtcombine} provides the retention 
time window where ion intensities will be summed for quantification. 
If \code{rtcombine = c(0,0)}, All scans will be used}
  
  \item{mzbin}{Bin size for the generation of mass vector}
  
  \item{center}{Indicates number of file which will be used for 
retention time alignment. If zero, file will be automatically 
selected based on maximum summed ion intensity of the combined mass 
spectrum.}
  
  \item{linear}{logical. If TRUE, linear vector will be generated 
with \code{mzbin} increments. If FALSE, mass vector will be 
generated using a non-linear function. This option is recommended 
for TOF-type mass detectors}
  
  \item{profStep}{step size (in m/z) to use for profile generation 
from the raw data files.}
  
  \item{retcorrect}{Logical, should retention time correction be 
used? \code{default = TRUE}}
  
  \item{continuum}{Logical, is continuum data used? \code{default = 
FALSE}}
  
  \item{scales}{Peak Detection of Mass Peaks in the combined mass 
spectrum: scales for continuous wavelet transformation (CWT). }
  
  \item{SNR.Th}{Peak Detection of Mass Peaks in the combined mass 
spectrum: Signal to noise ratio threshold}
  
  \item{SNR.area}{Peak Detection of Mass Peaks in the combined mass 
spectrum: Area around the peak for noise determination. Indicates 
number of surrounding peaks on the first CWT scale. \code{default = 
20}}
  
  \item{RTscales}{Peak Detection of Chromatographic peaks on the 
combined extracted ion chromatogram: scales for continuous wavelet 
transformation (CWT), see also the \code{peakDetectionCWT} function 
of the \code{MassSpecWavelet} package. }
  
  \item{RTSNR.Th}{Peak Detection of Chromatographic peaks on the 
combined extracted ion chromatogram: Signal to noise ratio 
threshold}
  
  \item{RTSNR.area}{Peak Detection of Chromatographic peaks on the 
combined extracted ion chromatogram: Area around the peak for noise 
determination. Indicates number of surrounding peaks on the first 
CWT scale. \code{default = 20}}
  
\item{mintr}{Minimal peak width intensity treshold (in percentage), 
for which two overlapping peaks are considered as separated. The 
default calue is 0.5.}

}

\details{Current data processing tools focus on a peak detection 
strategy where initially each LCMS run is treated independently from 
each other. Subsequent data processing for the alignment of 
different samples is then calculated on reduced peak tables. This 
function involves the merging of all LCMS datasets of a given 
experiment as a first step of raw data processing. The merged LCMS 
dataset contains an overlay and the sum of ion intensities from all 
LCMS runs. Peak detection is then performed only on this merged 
dataset and quantification of the signals in each sample is guided 
by the peak location in the merged data.

See also \code{xcms::retcor.obiwarp} and \code{cosmiq::create_datamatrix} which 
executes each step of the cosmiq function.


For running examples please read the package vignette or \code{?cosmiq::create_datamatrix}.

}

\author{

   Christan Panse, David Fischer 2014
   
}

\examples{

    \dontrun{
    # the following lines consume 2m17.877s  
    # on a Intel(R) Xeon(R) CPU X5650  @ 2.67GHz
    library(faahKO)

    cdfpath <- file.path(find.package("faahKO"), "cdf")
    my.input.files <- dir(c(paste(cdfpath, "WT", sep='/'), 
        paste(cdfpath, "KO", sep='/')), full.names=TRUE)

     x <- cosmiq(files=my.input.files, mzbin=0.25, SNR.Th=0, linear=TRUE)

     image(t(x$eicmatrix^(1/4)), 
        col=rev(gray(1:20/20)),
        xlab='rt', 
        ylab='m/z', axes=FALSE )

     head(xcms::peakTable(x$xcmsSet))
    }

}
