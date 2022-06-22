source("paramPlot.r")

params <- c("bindrate")
savefilename <- paste(params, collapse = "-")
filenames <- paramFilenames(params, subdir = "DIFFUSION_TRANSLATION")
main <- as.numeric(sapply(strsplit(filenames, "-"), "[[", 2))
compfnRealfdm <- "../fields/half-initp25-d40-idealFinal-fields-0.bz2"
compfnFakefdm <- "../fields/diffusionBindhalf0-d40-fakefdm"
#mainLabel <- expression(Association~rate~"("*M^-1*s^-1*")")
mainLabel <- "Association rate ($\\si{\\per\\Molar\\per\\second}$)"

mapply(paramPlotOne,
       plotname = c("bindRate-realfdm", "bindRate-fakefdm"),
       comparisonFilename = c(compfnRealfdm, compfnFakefdm),
       MoreArgs = list(
                       savefilename = savefilename,
                       filenames = filenames,
                       paramMain = main, 
                       xLabel = mainLabel,
                       xlog = TRUE,
                       controlName = "diffusion_translation-d40",
                       height = 2.1,
                       fancyLogLabel = TRUE,
                       skip = c(FALSE, TRUE)
       )
)
