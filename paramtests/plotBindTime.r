source("paramPlot.r")

params <- c("bindhalf", "half")
savefilename <- paste(params, collapse = "-")
filenames <- paramFilenames(paramList = params, subdir = "DIFFUSION_TRANSLATION")
bindhalfs <- as.numeric(sapply(strsplit(filenames, "-"), "[[", 2))
halfs <- as.numeric(sapply(strsplit(filenames, "-"), "[[", 4))
bindhalfToHalf <- bindhalfs / halfs
compfnRealfdm <- "../fields/half-initp25-d40-idealFinal-fields-0.bz2"
compfnFakefdm <- "../fields/diffusionBindhalf0-d40-fakefdm"
bindhalfLabel <- "Half-life of morphogen bound to sensors (minutes)"
halfLabel <- "Half-life of free morphogen (minutes)"
bindhalfToHalfLabel <- "Ratio of half-life of morphogen bound to sensor to half-life of free morphogen"

#If I give a symbol to morphogen half-life like \tau_{1/2}, I could maybe have room for a legend on the right. Though maybe other plots will need it on top anyway.

mapply(paramPlot,
       plotname = c("bindhalfByHalf-realfdm", "bindhalfByHalf-fakefdm"),
       comparisonFilename = c(compfnRealfdm, compfnFakefdm),
       MoreArgs = list(
                       savefilename = savefilename,
                       filenames = filenames,
                       paramMain = bindhalfs, 
                       paramBy = halfs, 
                       xLabel = bindhalfLabel,
                       byLabel = halfLabel,
                       controlName = "diffusion_translationBindhalf0-d40",
                       height = 2.7
       )
)

#paramPlot(plotname = "relBindhalfByHalf",
#          savefilename = savefilename,
#          filenames = filenames,
#          paramMain = bindhalfToHalf, 
#          paramBy = halfs, 
#          comparisonFilename = comparisonFilename,
#          xLabel = bindhalfToHalfLabel,
#          byLabel = halfLabel,
#          xlog = TRUE,
#          widthToHeight = 3)
