source("paramPlot.r")

params <- c("targetCorr")
filenames <- paramFilenames(paramList = params, subdir = "GROWTH")
savefilename <- paste(params, collapse = "-")
main <- sapply(strsplit(filenames, "-"), "[[", 2)
comparisonFilename <- "../fields/growth-d40-idealFinal-fields-0.bz2"
#mainLabel <- "Recruitment with and without target density correction"
#mainLabel <- "Target density correction"
#mainLabel <- " "

print(
      paramPlotOne(plotname = "recruitmentTarget",
                   savefilename = savefilename,
                   filenames = filenames,
                   paramMain = main, 
                   comparisonFilename = comparisonFilename,
                   xLabel = NULL,
                   xlog = FALSE,
                   bar = TRUE,
                   height = 2.1,
                   widthToHeight = 1,
                   controlName = "growth-d40",
                   controlLabel = "Corrected",
                   ticklabels = c("With", "Without")
      )
)
