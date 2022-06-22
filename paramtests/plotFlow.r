source("paramPlot.r")

params <- c("swimming")
savefilename <- paste(params, collapse = "-")
filenames <- paramFilenames(paramList = params, extraName = "diff", subdir = "DIFFUSION")
main <- sapply(strsplit(filenames, "-"), "[[", 2)
comparisonFilename <- "../fields/half-initp25-d80-lores-idealFinal-fields-0.bz2"

print(
      paramPlotOne(plotname = "flow",
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
