\name{create_datamatrix}

\alias{create_datamatrix}

\title{Quantifying mz/RT intensities using peak locations from master map}

\description{

     \code{create_datamatrix} identifies mz/RT peak boundaries in 
each raw file using the information from a master mass spectrum and 
master EIC. For each mz/RT location, the peak volume is calculated 
and stored in a report table.
}

\usage{

    create_datamatrix(xs,rxy)

}

\arguments{

  \item{xs}{xcmsSet object}

  \item{rxy}{matrix containing mz and RT boundaries for each 
detected peak}
}

\details{With the information about their position in the combined 
datasets, each individual mz/RT feature is located in the raw data.

}

\author{

   David Fischer 2013
   
}

\examples{
    \dontrun{
    cdfpath <- file.path(find.package("faahKO"), "cdf")
    my.input.files <- dir(c(paste(cdfpath, "WT", sep='/'),
        paste(cdfpath, "KO", sep='/')), full.names=TRUE)

    #
    # create xcmsSet object
    # todo
    xs <- new("xcmsSet")
    # consider only two files!!!
    xs@filepaths <- my.input.files[1:2]

    class<-as.data.frame(c(rep("KO",2),rep("WT", 0)))
    rownames(class)<-basename(my.input.files[1:2])
    xs@phenoData<-class

    x<-combine_spectra(xs=xs, mzbin=0.25,
        linear=TRUE, continuum=FALSE)

    plot(x$mz, x$intensity, type='l',
        xlab='m/Z', ylab='ion intensity')

    xy <- peakdetection(x=x$mz, y=x$intensity, scales=1:10,
    SNR.Th=0.0, SNR.area=20, mintr=0.5)

    id.peakcenter<-xy[,4]

    plot(x$mz, x$intensity, type='l',
        xlim=c(440,460),
            xlab='m/Z', ylab='ion intensity')

    points(x$mz[id.peakcenter], x$intensity[id.peakcenter],
        col='red', type='h')


    # create dummy object
    xs@peaks<-matrix(c(rep(1, length(my.input.files) * 6),
        1:length(my.input.files)), ncol=7)

    colnames(xs@peaks) <- c("mz", "mzmin", "mzmax", "rt",
        "rtmin", "rtmax", "sample")

    xs<-xcms::retcor(xs, method="obiwarp", profStep=1,
        distFunc="cor", center=1)



    eicmat<-eicmatrix(xs=xs, xy=xy, center=1)

    #  process a reduced mz range for a better package build performance
    (eicmat.mz.range<-range(which(475 < xy[,1] & xy[,1] < 485)))

    eicmat.filter <- eicmat[eicmat.mz.range[1]:eicmat.mz.range[2],]
    xy.filter <- xy[eicmat.mz.range[1]:eicmat.mz.range[2],]

    #
    # determine the new range and plot the mz versus RT map
    (rt.range <- range(as.double(colnames(eicmat.filter))))
    (mz.range<-range(as.double(row.names(eicmat.filter))))

    image(log(t(eicmat.filter))/log(2), col=rev(gray(1:20/20)),
        xlab='rt [in seconds]', ylab='m/z', axes=FALSE,
        main='overlay of 12 samples using faahKO')

    axis(1, seq(0,1, length=6), round(seq(rt.range[1], rt.range[2], length=6)))
    axis(2, seq(0,1, length=4), seq(mz.range[1], mz.range[2], length=4))

    #
    # determine the chromatographic peaks
    rxy<-retention_time(xs=xs, RTscales=c(1:10, seq(12,32, by=2)),
        xy=xy.filter,
        eicmatrix=eicmat.filter,
        RTSNR.Th=120, RTSNR.area=20)

    rxy.rt <- (rxy[,4] - rt.range[1])/diff(rt.range)
    rxy.mz <- (rxy[,1] - mz.range[1])/diff(mz.range)

    points(rxy.rt, rxy.mz, pch="X", lwd=2, col="red")


    xs<-create_datamatrix(xs=xs, rxy=rxy)

    peaktable <- xcms::peakTable(xs)

    head(peaktable)
    }
}
