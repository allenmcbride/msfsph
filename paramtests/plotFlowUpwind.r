source("paramPlot.r")

fieldName <- "agentField-0"
params <- c("swimming")
savefilename <- paste(params, collapse = "-")
filenames <- paramFilenames(paramList = params, fieldName = fieldName, extraName = "adv80", subdir = "ADVECTION")
main <- sapply(strsplit(filenames, "-"), "[[", 2)
comparisonFilename <- "../fields/upwind-d40-idealFinal-fields-0.bz2"

print(
      paramPlotOne(plotname = "flowUpwind",
                   savefilename = savefilename,
                   filenames = filenames,
                   paramMain = main, 
                   comparisonFilename = comparisonFilename,
                   xLabel = NULL,
                   bar = TRUE,
                   height = 2.1,
                   widthToHeight = 1,
                   fieldName = fieldName,
                   ticklabels = c("Swimming", "Crawling")
      )
)
