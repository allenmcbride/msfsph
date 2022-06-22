source("paramPlot.r")

findMainVals <- function(filenames) {
    return(as.numeric(sapply(strsplit(filenames, "-"), "[[", 2)) / 3600)
}

allThree <- function(plotnameSuffix, comparisonFilename) {
    paramsList <- list(c("sensormttf"), c("propulsionmttf"), c("productionmttf"))
    savefilenameList <- list("sensormttf", "propulsionmttf", "productionmttf")
    filenamesList <- sapply(paramsList, paramFilenames, subdir = "DIFFUSION_TRANSLATION")
    mainValsList <- sapply(filenamesList, findMainVals)
    mainLabelBase <- "MTTF (hours) for"
    mainLabelSpecific <- list("sensors", "propulsion", "production")
    mainLabelList <- mapply(paste, list(mainLabelBase), mainLabelSpecific)
    plotNamesList <- list("sensorMTTF", "propulsionMTTF", "productionMTTF")

    mapply(paramPlotOne,
           plotname = paste0(plotNamesList, plotnameSuffix),
           savefilename = savefilenameList,
           filenames = filenamesList,
           paramMain = mainValsList,
           xLabel = mainLabelList,
           MoreArgs = list(comparisonFilename = comparisonFilename,
                           xlog = TRUE,
                           controlName = "diffusion_translation-d40",
                           height = 1.5,
                           skip = c(TRUE, FALSE))
           )

    warnings()
}

mapply(allThree,
       plotnameSuffix = c("-realfdm", "-fakefdm"),
       comparisonFilename = c(
                              "../fields/half-initp25-d40-idealFinal-fields-0.bz2",
                              "../fields/diffusionBindhalf0-d40-fakefdm"
       )
)
