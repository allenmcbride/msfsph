source("paramPlot.r")

params <- c("speedbias")
savefilename <- paste(params, collapse = "-")
filenames <- paramFilenames(paramList = params, subdir = "DIFFUSION_TRANSLATION")
main <- as.numeric(sapply(strsplit(filenames, "-"), "[[", 2))
#filenames <- filenames[main > 0.01]
#main <- main[main > 0.01]
compfnRealfdm <- "../fields/half-initp25-d40-idealFinal-fields-0.bz2"
compfnFakefdm <- "../fields/diffusionBindhalf0-d40-fakefdm"
mainLabel <- "Speed bias $\\sigma$ (unitless)"

print(
      mapply(paramPlotOne,
             plotname = c("speedBias-realfdm", "speedBias-fakefdm"),
             comparisonFilename = c(compfnRealfdm, compfnFakefdm),
             MoreArgs = list(
                             savefilename = savefilename,
                             filenames = filenames,
                             paramMain = main, 
                             xLabel = mainLabel,
                             xlog = TRUE,
                             controlName = "diffusion_translation-d40",
                             controlLabel = "ctrl",
                             height = 1.5,
                             fancyLogLabel = TRUE,
                             skip = c(TRUE, FALSE)
             )
      )
)
