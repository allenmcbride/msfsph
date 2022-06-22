source("paramPlot.r")

params <- c("swimming")
savefilename <- paste(params, collapse = "-")
filenames <- paramFilenames(paramList = params, extraName = "grow", subdir = "GROWTH")
main <- sapply(strsplit(filenames, "-"), "[[", 2)
comparisonFilename <- "../fields/growth-d160-idealFinal-fields-0.bz2"

print(
      paramPlotOne(plotname = "flowGrow",
                   savefilename = savefilename,
                   filenames = filenames,
                   paramMain = main, 
                   comparisonFilename = comparisonFilename,
                   xLabel = NULL,
                   bar = TRUE,
                   height = 2.1,
                   widthToHeight = 1,
                   ticklabels = c("Swimming", "Crawling")
      )
)
