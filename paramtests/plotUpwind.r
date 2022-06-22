source("paramPlot.r")

fieldName <- "agentField-0"
name <- "upwind"
savefilename <- name
filenames <- paramFilenames(paramList = c(name), fieldName = fieldName, subdir = "ADVECTION")

print(
      paramPlotOne(plotname = name,
                   savefilename = savefilename,
                   filenames = filenames,
                   paramMain = sapply(strsplit(filenames, "-"), "[[", 2),
                   comparisonFilename = paste0("../fields/", name, "-d40-idealFinal-fields-0.bz2"),
                   #xLabel = "Maintaining fixed pattern while agents move, with and without upwind dot product",
                   xLabel = NULL,
                   xlog = FALSE,
                   bar = TRUE,
                   height = 2.1,
                   widthToHeight = 1,
                   fieldName = fieldName,
                   ticklabels = c("With", "Without")
      )
)
