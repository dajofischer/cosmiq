#R

# taken from xcms:::phenoDataFromPaths
.xcmsphenoDataFromPaths<-function (paths) 
{
    sclass <- gsub("^\\.$", "sample", dirname(paths))
    lev <- strsplit(sclass, "/")
    levlen <- sapply(lev, length)
    if (length(lev) > 1 && !all(levlen[1] == levlen)) 
        stop("Directory tree must be level")
    pdata <- as.data.frame(matrix(unlist(lev), nrow = length(lev), 
        byrow = TRUE))
    redundant <- apply(pdata, 2, function(col) length(unique(col)) == 
        1)
    if (!any(!redundant)) {
        redundant[length(redundant)] <- FALSE
    }
    pdata <- pdata[, !redundant, drop = FALSE]
    if (ncol(pdata) == 1) {
        scomp <- strsplit(substr(sclass, 1, min(nchar(sclass))), 
            "")
        scomp <- matrix(c(scomp, recursive = TRUE), ncol = length(scomp))
        i <- 1
        while (all(scomp[i, 1] == scomp[i, -1]) && i < nrow(scomp)) i <- i + 
            1
        i <- min(i, tail(c(0, which(scomp[1:i, 1] == .Platform$file.sep)), 
            n = 1) + 1)
        if (i > 1 && i <= nrow(scomp)) 
            sclass <- substr(sclass, i, max(nchar(sclass)))
        pdata <- data.frame(class = sclass)
    }
    rownames(pdata) <- gsub("\\.[^.]*$", "", basename(paths))
    pdata
}

cosmiq<-
function (files, RTinfo = TRUE, mzbin = 0.003, center = 0, linear = FALSE, 
    profStep = 1, retcorrect = TRUE, continuum = FALSE, rtcombine = c(0, 
        0), scales = c(1:10), SNR.Th = 10, SNR.area = 20, RTscales = c(1:10, 
        seq(12, 32, 2)), RTSNR.Th = 20, RTSNR.area = 20, mintr = 0.5) 
{
    xs <- new("xcmsSet")
    xs@filepaths <- files
    xs@phenoData <- .xcmsphenoDataFromPaths(files)
    x <- combine_spectra(xs = xs, mzbin = mzbin, linear = linear, 
        continuum = continuum)
    message("Detecting mz peaks on master spectrum ...")
    xy <- peakdetection(x = x[[1]], y = x[[2]], scales = scales, 
        SNR.Th = SNR.Th, SNR.area = SNR.area, mintr = mintr)
    if (RTinfo == FALSE) {
        xs <- quantify_combined(xs = xs, xy = xy, rtcombine = rtcombine)
        peakparams <- list(RTinfo = RTinfo, mzbin = mzbin, linear = linear, 
            continuum = continuum, rtcombine = rtcombine, scales = scales, 
            SNR.Th = SNR.Th, SNR.area = SNR.area)
        output <- list(xcmsSet = xs, masterspec = x, peakparams = peakparams, 
            mzpeaks = xy)
    }
    else {
        if (length(files) > 1 & retcorrect == TRUE) {
            xs@peaks <- matrix(c(rep(1, length(files) * 6), 1:length(files)), 
                ncol = 7)
            colnames(xs@peaks) <- c("mz", "mzmin", "mzmax", "rt", 
                "rtmin", "rtmax", "sample")
            if (center == 0) {
                maxintensity <- rep(0, length(files))
                for (i in 1:length(files)) {
                  temp <- xcmsRaw(files[i])
                  maxintensity[i] <- sum(temp@env$intensity)
                }
                center <- which.max(maxintensity)
            }
            xs <- xcms::retcor(xs, method = "obiwarp", profStep = profStep, 
                distFunc = "cor", center = center)
        }
        else {
            for (i in 1:length(files)) {
                temp <- xcmsRaw(files[i])
                xs@rt$corrected[[i]] <- temp@scantime
                xs@rt$raw[[i]] <- temp@scantime
            }
            center <- 1
        }
        eicmat <- eicmatrix(xs = xs, xy = xy, center = center)
        rxy <- retention_time(xs = xs, RTscales = RTscales, xy = xy, 
            eicmatrix = eicmat, RTSNR.Th = RTSNR.Th, RTSNR.area = RTSNR.area)
        xs <- create_datamatrix(xs = xs, rxy = rxy)
        peakparams <- list(RTinfo = RTinfo, mzbin = mzbin, center = center, 
            linear = linear, profStep = profStep, retcorrect = retcorrect, 
            continuum = continuum, rtcombine = rtcombine, scales = scales, 
            SNR.Th = SNR.Th, SNR.area = SNR.area, RTscales = RTscales, 
            RTSNR.Th = RTSNR.Th, RTSNR.area = RTSNR.area)
        output <- list(xcmsSet = xs, eicmatrix = eicmat, masterspec = x, 
            peakparams = peakparams, mzpeaks = xy)
    }
    return(output)
}

