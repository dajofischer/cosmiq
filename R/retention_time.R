#R


retention_time <- 
function (xy, xs, eicmatrix, RTscales = c(1:10, seq(12, 32, 2)), 
    RTSNR.Th = 10, RTSNR.area = 20, mintr = 0.5) 
{
    message("Detecting chromatographic peaks from EIC matrix ...")

    rxy <- zeros(1, 10)
    colnames(rxy) <- c("mz", "mzmin", "mzmax", "rt", "rtmin", 
        "rtmax", "xcenter", "xleft", "xright", "scale")
    for (k in 1:size(eicmatrix, 1)) {
        temp <- peakdetection(x = as.numeric(colnames(eicmatrix)), 
            y = eicmatrix[k, ], scales = RTscales, SNR.Th = RTSNR.Th, 
            SNR.area = RTSNR.area, mintr = mintr)
        if (isempty(temp) == FALSE) {
            mzcor <- t(matrix(rep(xy[k, 1:3]), 3, size(temp, 
                1)))
            rtcor <- temp[, c(1:6, 9)]
            if (is.matrix(rtcor) == FALSE) {
                rtcor <- t(as.matrix(rtcor))
            }
            temp <- cbind(mzcor, rtcor)
            colnames(temp) <- colnames(rxy)[1:10]
            rxy <- rbind(rxy, temp)
        }
        message(paste(k, "of ", size(eicmatrix, 1), "Peaknum = ", 
            size(rxy, 1)))
    }
    rxy <- rxy[2:size(rxy, 1), ]
    return(rxy)
}


