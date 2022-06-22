source("paramPlot.r")

params <- c("half", "diff")
filenames <- paramFilenames(paramList = params)
savefilename <- paste(params, collapse = "-")
halfs <- as.numeric(sapply(strsplit(filenames, "-"), "[[", 2))
diffs <- as.numeric(sapply(strsplit(filenames, "-"), "[[", 4))
#hs <- sqrt(diffs * 60.0 * halfs * (1000.0 ^ 2.0) / log(2.0))
#hs <- round(sqrt(diffs * 60.0 * halfs / log(2.0)))
hs <- sqrt(diffs * 60.0 * halfs / log(2.0))
comparisonFilename <- "../fields/half-initp25-d40-idealFinal-fields-0.bz2"
#diffLabel <- expression(Diffusivity~"("*mu*m^2/s*")")
diffLabel <- "Diffusivity ($\\si{\\square\\micro\\meter\\per\\second}$)"
halfLabel <- "Half-life (minutes)"
hLabel <- "Smoothing length (microns)"

paramPlot(plotname = "hByDiff",
          savefilename <- paste(params, collapse = "-"),
          filenames = filenames,
          paramMain = hs, 
          paramBy = diffs, 
          comparisonFilename = comparisonFilename,
          xLabel = hLabel,
          height = 2.7,
          byLabel = diffLabel,
          errorwidth = 0.008)

paramPlot(plotname = "hByHalf",
          savefilename <- paste(params, collapse = "-"),
          filenames = filenames,
          paramMain = hs, 
          paramBy = halfs, 
          comparisonFilename = comparisonFilename,
          xLabel = hLabel,
          byLabel = halfLabel,
          height = 2.7,
          signifLabel = TRUE,
          marginOffset = 20,
          errorwidth = 0.008)

paramPlot(plotname = "halfByDiff",
          savefilename <- paste(params, collapse = "-"),
          filenames = filenames,
          paramMain = halfs, 
          paramBy = diffs, 
          comparisonFilename = comparisonFilename,
          xLabel = halfLabel,
          byLabel = diffLabel,
          height = 2.7,
          marginOffset = 200)

#paramPlot(plotname = "diffByHalf",
#          savefilename <- paste(params, collapse = "-"),
#          filenames = filenames,
#          paramMain = diffs, 
#          paramBy = halfs, 
#          comparisonFilename = comparisonFilename,
#          xLabel = diffLabel,
#          byLabel = halfLabel,
#          xlog = TRUE)

#paramPlot(plotname = "diffByH",
#          savefilename <- paste(params, collapse = "-"),
#          filenames = filenames,
#          paramMain = diffs, 
#          paramBy = hs, 
#          comparisonFilename = comparisonFilename,
#          xLabel = diffLabel,
#          byLabel = hLabel,
#          xlog = TRUE)

#paramPlot(plotname = "halfByH",
#          savefilename <- paste(params, collapse = "-"),
#          filenames = filenames,
#          paramMain = halfs, 
#          paramBy = hs, 
#          comparisonFilename = comparisonFilename,
#          xLabel = halfLabel,
#          byLabel = hLabel,
#          xlog = TRUE)

