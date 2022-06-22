source("paramPlot.r")

params <- c("motionnoise")
filenames <- paramFilenames(params)
savefilename <- paste(params, collapse = "-")
main <- as.numeric(sapply(strsplit(filenames, "-"), "[[", 2))
#filenames <- filenames[main > 0.01]
#main <- main[main > 0.01]
comparisonFilename <- "../fields/growth-d40-idealFinal-fields-0.bz2"
mainLabel <- "Additive motion noise standard deviation (mm/s)"

paramPlotOne(plotname = "motionNoise",
             savefilename = savefilename,
             filenames = filenames,
             paramMain = main, 
             comparisonFilename = comparisonFilename,
             xLabel = mainLabel,
             xlog = TRUE,
             controlName = "growth-d40",
             height = 2.1,
             fancyLogLabel = TRUE,
             skip = c(FALSE, TRUE)
)
