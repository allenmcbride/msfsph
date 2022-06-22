source("paramPlot.r")

params <- c("half", "diam")
savefilename <- paste(params, collapse = "-")
#filenames <- paramFilenames(paramList = params, subdir = "DIFFUSION")
filenames <- paramFilenames(paramList = params, subdir = ".")
halfs <- as.numeric(sapply(strsplit(filenames, "-"), "[[", 2))
diams <- as.numeric(sapply(strsplit(filenames, "-"), "[[", 4))
diff <- 100.0
hs <- sqrt(diff * 60.0 * halfs / log(2.0))
relHs <- round(hs / diams, 1)
comparisonFilename <- paste0("../fields/half-initp25-d", diams, "-idealFinal-fields-0.bz2")
hLabel <- "Smoothing length ($\\si{\\micro\\meter}$)"
relHLabel <- "Ratio of smoothing length to diameter"
diamLabel <- "Diameter ($\\si{\\micro\\meter}$)"

print(
      paramPlot(plotname = "hByDiam",
                savefilename = savefilename,
                filenames = filenames,
                paramMain = hs, 
                paramBy = diams, 
                comparisonFilename = comparisonFilename,
                xLabel = hLabel,
                byLabel = diamLabel,
                signifLabel = TRUE,
                marginOffset = 110)
)

#paramPlot(plotname = "diamByH",
#          savefilename = savefilename,
#          filenames = filenames,
#          paramMain = diams, 
#          paramBy = hs, 
#          comparisonFilename = comparisonFilename,
#          xLabel = diamLabel,
#          byLabel = hLabel,
#          xlog = FALSE,
#          signifLabel = TRUE)

#paramPlot(plotname = "relHByDiam",
#          savefilename = savefilename,
#          filenames = filenames,
#          paramMain = relHs, 
#          paramBy = diams, 
#          comparisonFilename = comparisonFilename,
#          xLabel = relHLabel,
#          byLabel = diamLabel,
#          xlog = FALSE)
#
#paramPlot(plotname = "diamByRelH",
#          savefilename = savefilename,
#          filenames = filenames,
#          paramMain = diams, 
#          paramBy = relHs, 
#          comparisonFilename = comparisonFilename,
#          xLabel = diamLabel,
#          byLabel = relHLabel,
#          xlog = FALSE)
