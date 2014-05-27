#R



create_datamatrix <-
function (xs, rxy) 
{
    message("Quantifying mz/RT features ...")
    files <- xs@filepaths
    xs.length <- zeros(length(xs@rt$raw), 1)
    for (i in 1:length(xs@rt$raw)) {
        xs.length[i] <- length(xs@rt$raw[[i]])
    }
    for (i in 1:length(files)) {
        a <- xcmsRaw(files[i])
        dat <- zeros(size(rxy, 1), 3)
        colnames(dat) <- c("into", "maxo", "sample")
        dat[, 3] <- i
        dat <- cbind(rxy[, 1:6], dat)
        rtcor <- xs@rt$corrected[[i]]
        rtraw <- xs@rt$raw[[i]]
        rtcor[1] <- rtraw[1]
        rtcor[xs.length[i]] <- rtraw[xs.length[i]]
        for (k in 1:size(rxy, 1)) {
            if (rxy[k, 9] > xs.length[i]) {
                rxy[k, 9] <- xs.length[i]
            }
            if (rxy[k, 7] > xs.length[i]) {
                rxy[k, 7] <- xs.length[i]
            }
            temp <- rawEIC(a, mzrange <- rxy[k, 2:3])
            temp <- interp1(rtcor, temp$intensity, rtraw, method = "linear")
            tempcw <- MassSpecWavelet::cwt(temp, scales = rxy[k, 
                10], wavelet = "mexh")
            loc1 <- which(diff(sign(diff(tempcw))) == -2) + 1
            loc1.min <- abs(loc1 - rxy[k, 7])
            if (isempty(loc1.min) == 0) {
                c <- loc1[which(min(loc1.min) == loc1.min)]
                if (length(c) > 1) {
                  tempmax <- tempcw[c]
                  c <- c[max(tempmax) == tempmax]
                }
                if (length(c) > 1) {
                  c <- c[1]
                }
                retnew <- c - rxy[k, 7] + rxy[k, 7:9]
                if (sum(retnew < 1) > 0) {
                  retnew[retnew < 1] <- 1
                }
                if (sum(retnew > xs.length[i]) > 0) {
                  retnew[retnew > xs.length[i]] <- xs.length[i]
                }
                rtlim <- interp1(rtcor, rtraw, rtraw[retnew], 
                  method = "linear")
                pwid <- (rtlim[3] - rtlim[2])/(retnew[3] - retnew[2])
                dat[k, 7] <- pwid * sum(temp[retnew[2]:retnew[3]])
                dat[k, 8] <- max(temp[retnew[2]:retnew[3]])
                dat[k, 4:6] <- rtlim
            }
        }
        if (i == 1) {
            xs@peaks <- dat
        }
        else {
            xs@peaks <- rbind(xs@peaks, dat)
        }
        message(paste(i, "of", length(files), "done."))
    }
    groupidx <- matrix(seq(1, size(rxy, 1) * length(files)), 
        size(rxy, 1), length(files))
    groupidx <- split(groupidx, row(groupidx))
    xs@groupidx <- groupidx
    xsgroups <- zeros(size(rxy, 1), 1)
    xsgroups[, ] <- length(files)
    xsgroups <- cbind(rxy[, 1:6], xsgroups)
    colnames(xsgroups) <- c("mzmed", "mzmin", "mzmax", "rtmed", 
        "rtmin", "rtmax", "npeaks")
    xs@groups <- xsgroups
    return(xs)
}

