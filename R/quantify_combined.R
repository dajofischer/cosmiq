#R

quantify_combined <- 

function (xs, xy, rtcombine = c(0, 0)) 
{
    files <- xs@filepaths
    a <- xcmsRaw(files[1])
    if (rtcombine[1] == 0) {
        rtcombine[1] <- a@scantime[1]
    }
    if (rtcombine[2] == 0) {
        rtcombine[2] <- a@scantime[length(a@scantime)]
    }
    groupdata <- matrix(c(xy[, 1:3], rep(mean(rtcombine), nrow(xy)), 
        rep(rtcombine[1], nrow(xy)), rep(rtcombine[2], nrow(xy)), 
        rep(length(files), nrow(xy))), ncol = 7)
    message("Quantifying masses ...")
    for (j in 1:length(files)) {
        dat <- matrix(rep(0, nrow(xy) * 1), ncol = 1)
        if (j != 1) {
            a <- xcmsRaw(files[j])
        }
        for (k in 1:size(xy, 1)) {
            temp <- rawEIC(a, mzrange <- xy[k, 2:3], rtrange = rtcombine)
            dat[k, 1] <- sum(temp$intensity)
        }
        dat <- cbind(groupdata[, 1:6], dat)
        dat <- cbind(dat, rep(j, nrow(xy)))
        if (j == 1) {
            peakdata <- dat
        }
        else {
            peakdata <- rbind(peakdata, dat)
        }
        message(paste(j, "of", length(files), "done."))
    }
    colnames(groupdata) <- c("mzmed", "mzmin", "mzmax", "rtmed", 
        "rtmin", "rtmax", "npeaks")
    colnames(peakdata) <- c("mz", "mzmin", "mzmax", "rt", "rtmin", 
        "rtmax", "into", "sample")
    groupidx <- matrix(seq(1, nrow(xy) * length(files)), nrow(xy), 
        length(files))
    groupidx <- split(groupidx, row(groupidx))
    xs@groupidx <- groupidx
    xs@peaks <- peakdata
    xs@groups <- groupdata
    return(xs)
}

