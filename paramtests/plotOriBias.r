source("paramPlot.r")

params <- c("oribias")
savefilename <- paste(params, collapse = "-")
filenames <- paramFilenames(params, subdir = "DIFFUSION_TRANSLATION")
#main <- as.numeric(sapply(strsplit(filenames, "-"), "[[", 2))
mainFirst <- sapply(strsplit(filenames, "-"), "[[", 2)
mainSecond <- sapply(strsplit(filenames, "-"), "[[", 3)
mainString <- ifelse(substr(mainFirst, nchar(mainFirst), nchar(mainFirst)) == "e",
                     paste0(mainFirst, "-", mainSecond),
                     mainFirst)
main <- as.numeric(mainString)
comparisonFilenameRealfdm <- "../fields/half-initp25-d40-idealFinal-fields-0.bz2"
comparisonFilenameFakefdm <- "../fields/diffusionBindhalf0-d40-fakefdm"
mainLabel <- "Orientation bias $\\sigma$ (radians)"

mapply(paramPlotOne,
       plotname = c("oriBias-realfdm", "oriBias-fakefdm"),
       comparisonFilename = c(comparisonFilenameRealfdm, comparisonFilenameFakefdm),
       MoreArgs = list(
                       savefilename = savefilename,
                       filenames = filenames,
                       paramMain = main, 
                       xLabel = mainLabel,
                       xlog = TRUE,
                       controlName = "diffusion_translation-d40",
                       height = 1.5,
                       fancyLogLabel = TRUE,
                       skip = c(FALSE, TRUE, FALSE)
       )
)
