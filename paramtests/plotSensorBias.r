source("paramPlot.r")

params <- c("sensbias")
savefilename <- paste(params, collapse = "-")
filenames <- paramFilenames(params, subdir = "DIFFUSION_TRANSLATION")
biases <- as.numeric(sapply(strsplit(filenames, "-"), "[[", 2))
compfnRealfdm <- "../fields/half-initp25-d40-idealFinal-fields-0.bz2"
compfnFakefdm <- "../fields/diffusionBindhalf0-d40-fakefdm"
biasLabel <- "Sensor bias lognormal parameter"

mapply(paramPlotOne,
       plotname = c("sensorBias-realfdm", "sensorBias-fakefdm"),
       comparisonFilename = c(compfnRealfdm, compfnFakefdm),
       MoreArgs = list(
                       savefilename = savefilename,
                       filenames = filenames,
                       paramMain = biases, 
                       xLabel = biasLabel,
                       xlog = TRUE,
                       controlName = "diffusion_translation-d40",
                       height = 2.1,
                       fancyLogLabel = TRUE,
                       skip = c(FALSE, TRUE)
       )
)
