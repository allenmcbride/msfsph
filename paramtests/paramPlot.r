# paramPlot.r
# Allen McBride
# June 21, 2022
#
# The two main functions here are paramPlot() and paramPlotOne(). Both of
# these functions make plots from a replicate set of results, comparing
# each to idealized results. A set of replicate experimental results is
# specified by a list of filenames. The paramFilenames() function can
# be called by the same code that calls paramPlot*() to help find these
# filenames. By contrast, in cases where the idea of control replicates
# applies, these are specified by the beginning of the set of filenames,
# and paramPlot*() searches for all relevant files. In each case,
# simulation results (experimental or control) are compared with idealized
# results by calling functions from avgdist.r. These comparisons are saved
# in an .rds file for later re-use to speed the process of updating plots.
#
# paramPlot() is for multi-way (factorial) experiments and paramPlotOne()
# is for cases with a single varying parameter. The two functions
# contain duplicate code, but too much state would need to be passed
# to any helper functions to factor out the common code in an elegant
# way. A more object-oriented approach may be better in the future.

library(ggplot2)
library(viridis)
source("avgdist.r")

# This function is meant for interactive use; it is a hack for situations
# in which comparison functions expect idealized results to be made with
# FDM simulations, but actual idealized results are made with agent-based
# simulations. This function takes a set of replicates of the latter
# and averages the grids together.
prepareFakeFDM <- function(controlName, subdir = "") {
    filenames <- paramFilenames(controlName = controlName, subdir = subdir)
    print(filenames)
    data <- NULL
    for (fn in filenames) {
        if (is.null(data))
            data <- scan(fn, quiet=TRUE)
        else
            data <- data + scan(fn, quiet=TRUE)
    }
    data <- data / length(filenames)
    data <- c(0, 0, 0, 0, data) #Real FDM files have dimension info in first four values
    write(data, file = paste0(controlName, "-fakefdm"), ncolumns = 1)
}

paramFilenames <- function(paramList = c(), fieldName = "gaussianDens", controlName = "", extraName = NULL, subdir = "") {
    paramPattern <- if (length(paramList) > 0) paste(paste0(paramList, "-.*"), collapse = "") else ""
    extraPattern <- if (!is.null(extraName)) paste0("-", extraName) else ""
    fullPattern <- paste0("^", controlName, paramPattern, extraPattern, "-.*-", fieldName, "\\.bz2$")
    path <- paste0("../fields/", subdir)
    print(paste("paramList:", paramList))
    print(paste("filename pattern:", fullPattern))
    return(list.files(path = path, full.names = TRUE, pattern = fullPattern))
}

confidenceRadius <- function(data) {
    return(qt(0.975, length(data) - 1) * sd(data) / sqrt(length(data)))
}

# This function and the next make prettier axis labels for logarithmic x-axes.
log2labelOne <- function(val, fancyLogLabel, signifLabel) { 
    if (val < 1 && fancyLogLabel) 
        return(paste0("$2^{", round(log2(val)), "}$"))
    else
        if (signifLabel)
            return(signif(val, 3))
        else
            return(val)
}

log2labelfunc <- function(fancyLogLabel, signifLabel) {
    return(function(vals) {
               return(sapply(vals, log2labelOne, fancyLogLabel, signifLabel))
           })
}

findComparisons <- function(savefilename, requested, distfunc) {
    savefile <- paste0(savefilename, ".rds")

    # Create save file if needed
    if (!file.exists(savefile)) {
        savedComparisons <- data.frame(filename = character(), 
                                       comparisonFilename = character(), 
                                       mtimeSim = as.POSIXct(integer()), 
                                       mtimeComp = as.POSIXct(integer()), 
                                       comparison = numeric())
        saveRDS(savedComparisons, savefile)
    }

    # Recover saved info and figure out what is usable and what needs
    # updating Unlike the "requested" dataframe, the saved dataframe
    # will have NAs in the paramMain column, so we can use that to
    # distinguish any saved comparisons that aren't requested for this
    # particular plot.
    all <- merge(requested, readRDS(savefile), all = TRUE)
    rowNeedsComparison <- !is.na(all$paramMain) &
        (is.na(all$comparison) | 
         all$mtimeSim < file.mtime(all$filename) | 
         all$mtimeComp < file.mtime(all$comparisonFilename))

    # Make comparisons where needed, update timestamps
    if (any(rowNeedsComparison)) {
        all$comparison[rowNeedsComparison] <- mapply(distfunc, all$filename[rowNeedsComparison], all$comparisonFilename[rowNeedsComparison])
    }
    all$mtimeSim[rowNeedsComparison] <- file.mtime(all$filename[rowNeedsComparison])
    all$mtimeComp[rowNeedsComparison] <- file.mtime(all$comparisonFilename[rowNeedsComparison])

    # Re-save relevant columns of data
    saveRDS(all[c("filename", "comparisonFilename", "mtimeSim", "mtimeComp", "comparison")], savefile)

    return(all[!is.na(all$paramMain), ])
}

paramPlot <- function(plotname, savefilename, filenames, paramMain, paramBy, comparisonFilename, xLabel, byLabel, controlName = NULL, controlLabel = "", height = 3.0, widthToHeight = 1.7, fieldName = "gaussianDens", extraName = NULL, fancyLogLabel = FALSE, signifLabel = FALSE, marginOffset = 60, skip = c(TRUE), errorwidth = 0.015, order = NULL) {

    distfunc <- ifelse(fieldName == "gaussianDens", avgDist, avgDistAgents)

    if (!is.null(controlName)) {
        controlFilenames <- paramFilenames(fieldName = fieldName, controlName = controlName, extraName = extraName)
        controlRequested <- data.frame(filename = controlFilenames,
                                comparisonFilename = comparisonFilename,
                                paramMain = "dummy")
        controls <- findComparisons(controlName, controlRequested, distfunc)
        controlMean <- mean(controls$comparison)
        controlLine <- geom_hline(yintercept = controlMean, linetype = "dashed", size = 0.5 * height / 3.0)
    } else {
        controlLine = NULL
    }

    mygroup = if (is.null(order)) as.factor(paramBy) else factor(paramBy, levels = order)
    requested <- data.frame(filename = filenames,
                            comparisonFilename = comparisonFilename,
                            paramMain = paramMain,
                            paramBy = mygroup)

    comparisons <- findComparisons(savefilename, requested, distfunc)

    nRepsFrame <- aggregate(data.frame(comparisons$comparison), by = list(comparisons$paramMain), length)
    names(nRepsFrame) <- c("Level", "Number of replicates")
    print(nRepsFrame)

    # Compute means and standard errors, aggregating over replicates
    aggregated <- aggregate(data.frame(comparisons$comparison), by = list(comparisons$paramMain, comparisons$paramBy), mean)
    names(aggregated) <- c("paramMain", "paramBy", "mean")
    aggregated$confidenceRadius <- aggregate(comparisons$comparison, list(comparisons$paramMain, comparisons$paramBy), confidenceRadius)$x

    # Plot comparisons
    yLabel <- "Distance from ideal"
    mybreaks = sort(unique(aggregated$paramMain))[skip]
    mywidth = errorwidth * nlevels(as.factor(paramMain))
    ggplot(aggregated, aes(x = paramMain, y = mean, group = paramBy)) +
        scale_x_continuous(trans = "log2", breaks = mybreaks, labels = log2labelfunc(fancyLogLabel, signifLabel)) +
        ylim(0, NA) +
        labs(x = xLabel, y = yLabel) + 
        geom_line(aes(color = paramBy), size = 0.4 * height / 3.0) +
        geom_point(aes(color = paramBy), size = 0.5 * height / 3.0) +
        controlLine +
        geom_errorbar(aes(ymin = mean - confidenceRadius, ymax = mean + confidenceRadius, width = mywidth, color = paramBy), size = 0.4 * height / 3.0) + 
        scale_color_viridis(option = "turbo", discrete = TRUE, name = byLabel) +
        theme(legend.position="top") +
        theme(legend.title = element_text(margin = margin(r = -marginOffset)))
    ggsave(paste0("paramPlot-", plotname, ".svg"), width = widthToHeight * height, height = height)
    return(aggregated)
}

paramPlotOne <- function(plotname, savefilename, filenames, paramMain, comparisonFilename, xLabel, xlog = FALSE, bar = FALSE, controlName = NULL, controlLabel = "", height = 3.0, widthToHeight = 1.7, fieldName = "gaussianDens", extraName = NULL, fancyLogLabel = FALSE, signifLabel = FALSE, ticklabels = NULL, skip = c(TRUE)) {

    distfunc <- ifelse(fieldName == "gaussianDens", avgDist, avgDistAgents)
    
    if (!is.null(controlName)) {
        controlFilenames <- paramFilenames(fieldName = fieldName, controlName = controlName, extraName = extraName)
        controlRequested <- data.frame(filename = controlFilenames,
                                comparisonFilename = comparisonFilename,
                                paramMain = "dummy")
        controls <- findComparisons(controlName, controlRequested, distfunc)
        controlMean <- mean(controls$comparison)
        controlLine <- geom_hline(yintercept = controlMean, linetype = "dashed", size = 0.5 * height / 3.0)
    } else {
        controlLine = NULL
    }

    requested <- data.frame(filename = filenames,
                            comparisonFilename = comparisonFilename,
                            paramMain = paramMain)

    comparisons <- findComparisons(savefilename, requested, distfunc)

    nRepsFrame <- aggregate(data.frame(comparisons$comparison), by = list(comparisons$paramMain), length)
    names(nRepsFrame) <- c("Level", "Number of replicates")
    print(nRepsFrame)

    # Compute means and standard errors, aggregating over replicates
    aggregated <- aggregate(data.frame(comparisons$comparison), by = list(comparisons$paramMain), mean)
    names(aggregated) <- c("paramMain", "mean")
    aggregated$confidenceRadius <- aggregate(comparisons$comparison, list(comparisons$paramMain), confidenceRadius)$x

    # Plot comparisons
    yLabel <- "Distance from ideal"

    paramPlot <- 
        if (bar) {
            if (!is.null(controlName)) {
                controlConf <- confidenceRadius(controls$comparison)
                aggregated <- rbind(list(controlLabel, controlMean, controlConf), aggregated)
            }
            if (identical(aggregated$paramMain, c("false", "true")))
                aggregated <- aggregated[c(2, 1), ]
            ggplot(aggregated, aes(x = as.factor(paramMain), y = mean)) +
                scale_x_discrete(limits = aggregated$paramMain, labels = if(is.null(ticklabels)) waiver() else ticklabels) +
                ylim(0, NA) +
                labs(x = xLabel, y = yLabel) + 
                geom_col() +
                geom_errorbar(aes(ymin = mean - confidenceRadius, ymax = mean + confidenceRadius, width = .015 * nlevels(as.factor(paramMain))))
        } else {
            mybreaks = sort(unique(aggregated$paramMain))[skip]
            ggplot(aggregated, aes(x = paramMain, y = mean)) +
                scale_x_continuous(trans = if (xlog) "log2" else "identity", 
                                   breaks = mybreaks,
                                   labels = if (xlog) log2labelfunc(fancyLogLabel, signifLabel) else signif(mybreaks, 2)) +
                ylim(0, NA) +
                labs(x = xLabel, y = yLabel) + 
                geom_line(size = 0.4 * height / 3.0) +
                geom_point(size = 0.3 * height / 3.0) +
                controlLine +
                geom_errorbar(aes(ymin = mean - confidenceRadius, ymax = mean + confidenceRadius, width = .015 * nlevels(as.factor(paramMain))), size = 0.4 * height / 3.0)
        }

    ggsave(filename = paste0("paramPlot-", plotname, ".svg"), plot = paramPlot, width = widthToHeight * height, height = height)
    return(aggregated)
}
