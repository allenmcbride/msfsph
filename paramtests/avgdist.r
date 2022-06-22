# avgdist.r
# Allen McBride
# June 21, 2022
#
# This file contains functions related to comparing two results. The
# first two functions compute an average distance between two results
# and are used in paramPlot.r. avgDist() compares results pixel-by-pixel
# and is intended for comparing agent density with idealized density,
# in which case the simulator outputs a kernel density estimate for each
# grid cell. avgDistAgents() compares results at each agent location with
# idealized results at the same location; it is intended for evaluating
# the accuracy of "abtract" (non-density) fields. The remaining display*()
# functions are intended for interactive use.

avgDist <- function(fnAgents, fnFDM) {
    cat("Average distance between: ", fnAgents, fnFDM, ": ")
    agents <- scan(fnAgents, quiet=TRUE)
    fdm <- scan(fnFDM, quiet=TRUE, skip=4)
    if (length(agents) != length(fdm))
        stop(paste0("Length of ", fnAgents, ": ", length(agents), "; length of ", fnFDM, ": ", length(fdm)))
    validPix <- (agents > 0)
    dist <- sum(abs(fdm[validPix] - agents[validPix])) / length(validPix)
    cat(dist, "\n")
    return(dist)
}

avgDistAgents <- function(fnAgents, fnFDM) {
    cat("Average distance between: ", fnAgents, fnFDM, ": ")
    agents <- read.table(fnAgents, col.names=c("x", "y", "val"))
    fdmPixX <- scan(fnFDM, n=1, skip=0)
    fdmPixY <- scan(fnFDM, n=1, skip=1)
    fdmPhysX <- scan(fnFDM, n=1, skip=2)
    fdmPhysY <- scan(fnFDM, n=1, skip=3)
    fdmVec <- scan(fnFDM, skip=4)
    if (length(fdmVec) != fdmPixX * fdmPixY)
        stop(paste("Expected", fdmPixX * fdmPixY, "values in file", fnFDM, "; found", length(fdmVec)))
    rows <- floor(agents$y * fdmPixY / fdmPhysY)
    cols <- floor(agents$x * fdmPixX / fdmPhysX)
    indices <- rows * fdmPixX + cols + 1 #R is one-indexed
    dists <- abs(agents$val - fdmVec[indices])
    return(mean(dists))
}

displaySquareField <- function(fn, cutoff = 0.0, ratio = 1, skip = 0) {
    dev.new()
    field <- scan(fn, quiet=TRUE, skip=skip)
    if (cutoff == 0.0) cutoff <- max(field)
    # In next line, byrow = TRUE would seem correct because the field
    # files are row-major, but image() displays the matrix transposed
    # so it turns out wrong that way.
    fieldMatrix <- matrix(field, nrow = sqrt(ratio * length(field)))
    image(fieldMatrix, breaks = seq(0.0, cutoff, length.out = 13), useRaster = TRUE, main = fn)
}

displayHistogram <- function(fn, nBins, cutoff = 0.0) {
    dev.new()
    field <- scan(fn, quiet=TRUE)
    if (cutoff == 0.0) cutoff <- max(field)
    field <- field[field < cutoff]
    hist(field, breaks = seq(0.0, cutoff, length.out = nBins), main = fn)
}

displayAgentField <- function(fnAgents, fnFDM) {
    agents <- read.table(fnAgents, col.names=c("x", "y", "val"))
    fdmPhysX <- scan(fnFDM, n=1, skip=2)
    fdmPhysY <- scan(fnFDM, n=1, skip=3)
    plot(x=0, xlim=c(0, fdmPhysX), ylim=c(0, fdmPhysY), type='n')
    points(agents$x, agents$y, pch = 19, col = rgb(0.0, 0.0, (agents$val / max(agents$val))))
}
