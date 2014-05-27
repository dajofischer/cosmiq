#R

eicmatrix <-
function (xs, xy, center) 
{
    message("Creating EIC matrix ...")
    files <- xs@filepaths
    xs.length <- zeros(length(xs@rt$raw), 1)
    for (i in 1:length(xs@rt$raw)) {
        xs.length[i] <- length(xs@rt$raw[[i]])
    }
    eicmatrix <- zeros(nrow(xy), max(xs.length))
    for (i in 1:length(files)) {
        a <- xcmsRaw(files[i])
        temp.l <- 1:xs.length[i]
        rtcor <- xs@rt$corrected[[i]]
        rtraw <- xs@rt$raw[[i]]
        rtcor[1] <- rtraw[1]
        rtcor[xs.length[i]] <- rtraw[xs.length[i]]
        for (k in 1:size(xy, 1)) {
            temp <- rawEIC(a, mzrange <- xy[k, 2:3])
            eicmatrix[k, temp.l] <- eicmatrix[k, temp.l] + interp1(rtcor, 
                temp$intensity, rtraw, method = "linear")
        }
        message(paste(i, "of", length(files), "done."))
    }
    header <- xs@rt$raw[[center]]
    if (length(header) < max(xs.length)) {
        b <- which(xs.length == max(xs.length))
        b <- b[1]
        headerfill <- xs@rt$raw[[b]]
        header <- c(header, headerfill[(length(header) + 1):length(headerfill)])
    }
    colnames(eicmatrix) <- header
    rownames(eicmatrix) <- xy[, 1]
    return(eicmatrix)
}


