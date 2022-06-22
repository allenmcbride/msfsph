source("paramPlot.r")

fieldName <- "agentField-0"
params <- c("laplcorr")
savefilename <- paste(params, collapse = "-")
filenames <- paramFilenames(paramList = params, fieldName = fieldName)

print(
      paramPlotOne(plotname = "laplcorr",
                   savefilename = savefilename,
                   filenames = filenames,
                   paramMain = sapply(strsplit(filenames, "-"), "[[", 2),
                   comparisonFilename = "../fields/laplcorr-idealFinal-fields-0.bz2",
                   #xLabel = "Heated disc in bath, with and without Laplacian correction",
                   xLabel = NULL,
                   xlog = FALSE,
                   bar = TRUE,
                   height = 2.1,
                   widthToHeight = 1,
                   fieldName = fieldName,
                   ticklabels = c("With", "Without")
      )
)
