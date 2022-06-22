source("paramPlot.r")

params <- c("ecc")
savefilename <- paste(params, collapse = "-")
filenames <- paramFilenames(paramList = params, subdir = "DIFFUSION")
main <- sapply(strsplit(filenames, "-"), "[[", 2)
compfnVals <- paste0("../fields/half-initp25-mass-", c("", "ecc0.8-", "ecc0.96-"), "idealFinal-fields-0.bz2")
compfnMap <- c("0.0" = compfnVals[1], "0.8" = compfnVals[2], "0.96" = compfnVals[3])
compfns <- compfnMap[main]

mainLabel <- "Eccentricity"

print(
      paramPlotOne(plotname = "ecc",
                   savefilename = savefilename,
                   filenames = filenames,
                   paramMain = main, 
                   comparisonFilename = compfns,
                   xLabel = mainLabel,
                   bar = TRUE,
                   height = 2.1,
                   widthToHeight = 1.5
      )
)
