#R
combine_spectra <- 
function (xs, mzbin = 0.003, linear = FALSE, continuum = FALSE) 
{
    my.files <- xs@filepaths
    a <- xcmsRaw(my.files[1])
    mzrange <- a@mzrange
    if (linear == TRUE) {
        range <- seq(mzrange[1], mzrange[2], mzbin)
    }
    else {
        p1 = 2/3
        p2 = 1/3
        xr = seq(0, 1, mzbin/(mzrange[2] - mzrange[1]))
        range = (p1 * (xr^2)) + (p2 * xr)
        range = mzrange[1] + range * (mzrange[2] - mzrange[1])
    }
    if (continuum == TRUE) {
        range <- sort(unique(a@env$mz))
        range <- range[c(1:(length(range) - 1))] 
            + ((range[c(2:length(range))] 
            - range[c(1:(length(range) - 1))])/2)
        range <- c(mzrange[1], range, mzrange[2])
    }
    y <- as.vector(matrix(0, 1, (length(range) - 1)))
    message("Combining spectra ...")
    for (j in 1:length(my.files)) {
        if (j != 1) {
            a <- xcmsRaw(my.files[j])
        }
        A <- a@env$mz >= mzrange[1] & a@env$mz <= mzrange[2]
        a@env$mz <- a@env$mz[A]
        a@env$intensity <- a@env$intensity[A]
        b <- sort(a@env$mz, index.return = TRUE)
        b <- b$ix
        c <- hist(a@env$mz[b], plot = FALSE, range)
        index <- 1
        for (i in 1:length(c$counts)) {
            if (c$counts[i] > 0) {
                y[i] <- y[i] + sum(a@env$intensity[b[index:(index + 
                  c$counts[i])]])
                index <- index + c$counts[i]
            }
        }
        message(paste(j, "of", length(my.files), "done."))
    }
    A <- is.na(y) == FALSE
    y <- y[A]
    c$mi <- c$mi[A]
    x <- list(mz = c$mi, intensity = y)
    return(x)
}


