source("paramPlot.r")

params <- c("half", "speedup")
savefilename <- paste(params, collapse = "-")
filenames <- paramFilenames(params, subdir = "DIFFUSION_TRANSLATION")
halfs <- as.numeric(sapply(strsplit(filenames, "-"), "[[", 2))
diffs <- round(372.15416687891657557517^2 * (log(2)/(60*halfs)))
speedups <- as.numeric(sapply(strsplit(filenames, "-"), "[[", 4))
#relSpeedups <- signif(speedups / diffs, 3)
relSpeedups <- speedups / diffs
compfnRealfdm <- "../fields/half-initp25-d40-idealFinal-fields-0.bz2"
compfnFakefdm <- "../fields/diffusionBindhalf0-d40-fakefdm"
speedupLabel <- "Speed multiplier"
relSpeedupLabel <- "Speed multiplier (unitless) / diffusivity ($\\si{\\square\\micro\\meter\\per\\second}$)"
halfLabel <- "Half-life (minutes)"
diffLabel <- "Diffusivity ($\\si{\\square\\micro\\meter\\per\\second}$)"

mapply(paramPlot,
       plotname = c("speedByDiff-realfdm", "speedByDiff-fakefdm"),
       comparisonFilename = c(compfnRealfdm, compfnFakefdm),
       MoreArgs = list(
                       savefilename = savefilename,
                       filenames = filenames,
                       paramMain = speedups, 
                       paramBy = diffs, 
                       xLabel = speedupLabel,
                       byLabel = diffLabel,
                       height = 2.7
       )
)

mapply(paramPlot,
       plotname = c("relSpeedByDiff-realfdm", "relSpeedByDiff-fakefdm"),
       comparisonFilename = c(compfnRealfdm, compfnFakefdm),
       MoreArgs = list(
                       savefilename = savefilename,
                       filenames = filenames,
                       paramMain = relSpeedups, 
                       paramBy = diffs, 
                       xLabel = relSpeedupLabel,
                       byLabel = diffLabel,
                       height = 2.7,
                       signifLabel = TRUE,
                       marginOffset = 200,
                       skip = c(FALSE, TRUE)
       )
)

#mapply(paramPlot,
#       plotname = c("speedByHalf-realfdm", "speedByHalf-fakefdm"),
#       comparisonFilename = c(compfnRealfdm, compfnFakefdm),
#       MoreArgs = list(
#                       savefilename = savefilename,
#                       filenames = filenames,
#                       paramMain = speedups, 
#                       paramBy = halfs, 
#                       xLabel = speedupLabel,
#                       byLabel = halfLabel,
#                       xlog = TRUE
#       )
#)

#paramPlot(plotname = "halfBySpeed",
#          savefilename = savefilename,
#          filenames = filenames,
#          paramMain = halfs, 
#          paramBy = speedups, 
#          comparisonFilename = comparisonFilename,
#          xLabel = halfLabel,
#          byLabel = speedupLabel,
#          xlog = TRUE)
