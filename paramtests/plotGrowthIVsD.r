source("paramPlot.r")

filenames <- paramFilenames(controlName = "leakhalf-3840", subdir = "GROWTH")
main <- rep("Proportional", length(filenames))
comparisonFilename <- "../fields/growth-d40-idealFinal-fields-0.bz2"
mainLabel <- NULL

print(
      paramPlotOne(plotname = "growthIVsD",
                   savefilename = "growthIVsD",
                   filenames = filenames,
                   paramMain = main, 
                   comparisonFilename = comparisonFilename,
                   xLabel = mainLabel,
                   xlog = FALSE,
                   bar = TRUE,
                   height = 2.1,
                   widthToHeight = 1,
                   controlName = "growth-d40",
                   controlLabel = "Integral"
      )
)
