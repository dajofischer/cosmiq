#R

.xcmsDescendMin<-function (y, istart = which.max(y)) 
{
    if (!is.double(y)) 
        y <- as.double(y)

    unlist(.C(cosmiq_DescendMin, y, length(y), 
	as.integer(istart - 1), 
	ilower = integer(1), 
	iupper = integer(1))[4:5]) + 1
}


peakdetection <-
function (scales = c(1:10, seq(12, 32, 2)), y, x, SNR.Th = 20, 
    SNR.area = 20, mintr = 0.5) 
{
    message("Detecting mz peaks on master spectrum ...")
    cw <- MassSpecWavelet::cwt(y, scales = scales, wavelet = "mexh")
    cw <- t(cw)
    wsize <- size(cw)
    A <- cw >= cw[c(1, 1:wsize[1] - 1), ] & cw >= cw[c(1, 1:wsize[1] - 
        1), c(1, 1:wsize[2] - 1)] & cw >= cw[, c(1, 1:wsize[2] - 
        1)] & cw >= cw[c(2:wsize[1], wsize[1]), c(2:wsize[2], 
        wsize[2])] & cw >= cw[c(2:wsize[1], wsize[1]), c(1, 1:wsize[2] - 
        1)] & cw >= cw[c(1, 1:wsize[1] - 1), c(2:wsize[2], wsize[2])] & 
        cw >= cw[, c(2:wsize[2], wsize[2])] & cw >= cw[c(2:wsize[1], 
        wsize[1]), ]
    loc1 <- which(diff(sign(diff(cw[1, ]))) == -2) + 1
    A[1, loc1] <- TRUE
    loc1 <- matrix(c(loc1, cw[1, loc1], rep(0, length(loc1) * 
        2)), ncol = 4)
    xy <- meshgrid(1:wsize[2], y = 1:wsize[1])
    xy <- matrix(c(xy[[1]][A], xy[[2]][A], cw[A], rep(0, 2 * 
        sum(A))), ncol = 5)
    for (i in 1:nrow(loc1)) {
        loc1[i, 3] <- i - SNR.area
        loc1[i, 4] <- i + SNR.area
    }
    loc1[loc1[, 3] < 1, 3] <- 1
    loc1[loc1[, 4] > nrow(loc1), 4] <- nrow(loc1)
    xyacc <- rbind(matrix(c(xy[, 1], rep(0, 2 * nrow(xy))), ncol = 3), 
        matrix(c(loc1[, 1], rep(1, 2 * nrow(loc1))), ncol = 3))
    b <- sort(xyacc[, 1], index.return = TRUE)
    xyacc <- xyacc[b[[2]], ]
    xyacc[, 3] <- cumsum(xyacc[, 3])
    xyacc[xyacc[, 3] == 0, 3] <- 1
    xyacc <- xyacc[xyacc[, 2] == 0, 3]
    for (i in 1:nrow(xy)) {
        temp <- loc1[loc1[xyacc[i], 3]:loc1[xyacc[i], 4], 2]
        temp <- temp[temp > 1]
        if (isempty(temp) == 0) {
            xy[i, 5] <- xy[i, 3]/median(temp)
        }
    }
    if (SNR.Th > 0) {
        if (sum(xy[, 5] > SNR.Th) != 1) {
            xy <- xy[xy[, 5] > SNR.Th, ]
        }
        else if (sum(xy[, 5] > SNR.Th) == 1) {
            xy <- t(as.matrix(xy[xy[, 5] > SNR.Th, ]))
        }
    }
    if (isempty(xy) == FALSE) {
        xy <- matrix(c(c(xy), c(zeros(size(xy, 1), 3))), size(xy, 
            1), 8)
        for (i in 1:nrow(xy)) {
            xy[i, 6:7] <- .xcmsDescendMin(cw[xy[i, 2], ], istart = xy[i, 
                1])
            xy[i, 8] <- sum(y[c(xy[i, 1] - scales[1]):c(xy[i, 
                1] + scales[1])])
        }
        b <- sort.int(xy[, 2], decreasing = TRUE, index.return = TRUE)
        b <- b[[2]]
        xy <- xy[b, ]
        if (is.matrix(xy) == FALSE) {
            xy <- t(as.matrix(xy))
        }
        for (i in seq(1, nrow(xy))) {
            A <- which(xy[, 7] >= xy[i, 6] & xy[, 6] <= xy[i, 
                7] & xy[, 1] != 0)
            if (sum(A) > 1) {
                if (sum(y[xy[i, 6]:xy[i, 7]]) != 0) {
                  v <- c(min(xy[A, 6]):max(xy[A, 7]))
                  yv <- y[v] - min(y[v])
                  yvmin = yv < max(yv) * mintr
                  B <- cumsum(yvmin)
                  C <- hist(B[yvmin == 0], plot = FALSE, 0:max(B))
                  cwest <- round(max(C$counts)/2)
                  cwest <- abs(cwest - scales)
                  cwest <- which(cwest == min(cwest))
                  cwest <- cwest[1]
                  newmax <- v[which(diff(sign(diff(cw[cwest, 
                    v]))) == -2) + 1]
                  if (isempty(newmax)) {
                    A <- A[2:length(A)]
                    xy[A, ] <- 0
                  }
                  else {
                    xy[A, ] <- 0
                    xynew <- zeros(length(newmax), 8)
                    for (j in 1:length(newmax)) {
                      xynew[j, 1] <- newmax[j]
                      xynew[j, 2] <- cwest
                      xynew[j, 3] <- cw[cwest, newmax[j]]
                      xynew[j, 6:7] <- .xcmsDescendMin(cw[cwest, 
                        ], istart = newmax[j])
                      temp <- abs(xynew[j, 1] - loc1[, 1])
                      temp <- which(min(temp) == temp)
                      temp <- temp[1]
                      vsn <- (temp - SNR.area):(temp + SNR.area)
                      vsn <- vsn[vsn > 0 & vsn < nrow(loc1)]
                      vsn <- loc1[vsn, 2]
                      vsn <- vsn[vsn > 1]
                      if (isempty(vsn) == 0) {
                        xynew[j, 5] <- xynew[j, 3]/median(vsn)
                      }
                      if (SNR.Th > 0) {
                        if (xynew[j, 5] <= SNR.Th) {
                          xynew[j, ] <- 0
                        }
                      }
                    }
                    xy <- rbind(xy, xynew)
                  }
                }
                else {
                }
            }
        }
        if (sum(xy[, 1] != 0) != 1) {
            xy <- xy[xy[, 1] != 0, ]
        }
        else if (sum(xy[, 1] != 0) == 1) {
            xy <- t(as.matrix(xy[xy[, 1] != 0, ]))
        }
        b <- sort.int(xy[, 1], decreasing = FALSE, index.return = TRUE)
        b <- b[[2]]
        xy <- xy[b, ]
        if (is.matrix(xy) == FALSE) {
            xy <- t(as.matrix(xy))
        }
        xy <- matrix(c(x[xy[, 1]], x[xy[, 6]], x[xy[, 7]], xy[, 
            1], xy[, 6], xy[, 7], xy[, 3], xy[, 5], scales[xy[, 
            2]]), size(xy, 1), 9)
        colnames(xy) <- c("center", "left", "right", "xcenter", 
            "xleft", "xright", "signal cw", "SN", "scale")
    }
    return(xy)
}

