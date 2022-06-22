source("paramPlot.r")

params <- c("leakhalf")
filenames <- paramFilenames(paramList = params, subdir = "GROWTH")
savefilename <- paste(params, collapse = "-")
main <- as.numeric(sapply(strsplit(filenames, "-"), "[[", 2)) / 60
comparisonFilename <- "../fields/growth-d40-idealFinal-fields-0.bz2"
mainLabel <- "Half-life of leaky integral (hours)"

print(
      paramPlotOne(plotname = "leakhalf",
                   savefilename = savefilename,
                   filenames = filenames,
                   paramMain = main, 
                   comparisonFilename = comparisonFilename,
                   xLabel = mainLabel,
                   xlog = TRUE,
                   height = 1.5,
                   skip = c(TRUE, FALSE)
      )
)

print(
      paramPlotOne(plotname = "leakhalfsub",
                   savefilename = savefilename,
                   filenames = filenames[main >= 4],
                   paramMain = main[main >= 4], 
                   comparisonFilename = comparisonFilename,
                   xLabel = mainLabel,
                   xlog = TRUE,
                   height = 1.5,
                   skip = c(TRUE, FALSE)
      )
)
