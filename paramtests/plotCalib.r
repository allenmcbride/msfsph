source("paramPlot.r")

params <- c("speedup")
filenamesCalib <- paramFilenames(paramList = c("half-16-speedup"), subdir = "DIFFUSION_TRANSLATION")
filenamesNoVal <- paramFilenames(paramList = params, extraName = "noCalibVal", subdir = "DIFFUSION_TRANSLATION")
filenamesNoGrad <- paramFilenames(paramList = params, extraName = "noCalibGrad", subdir = "DIFFUSION_TRANSLATION")
filenamesNoCalib <- paramFilenames(paramList = params, extraName = "nocalib", subdir = "DIFFUSION_TRANSLATION")
filenames <- c(filenamesCalib, filenamesNoVal, filenamesNoGrad, filenamesNoCalib)
speedupCalib <- as.numeric(sapply(strsplit(filenamesCalib, "-"), "[[", 4))
speedupNoVal <- as.numeric(sapply(strsplit(filenamesNoVal, "-"), "[[", 2))
speedupNoGrad <- as.numeric(sapply(strsplit(filenamesNoGrad, "-"), "[[", 2))
speedupNoCalib <- as.numeric(sapply(strsplit(filenamesNoCalib, "-"), "[[", 2))
speedup <- c(speedupCalib, speedupNoVal, speedupNoGrad, speedupNoCalib)
type <- c(rep("Full calib.", length(filenamesCalib)),
          rep("Grad. only", length(filenamesNoVal)),
          rep("Conc. only", length(filenamesNoGrad)),
          rep("No calib.", length(filenamesNoCalib)))
savefilename <- "calib"
compfnRealfdm <- "../fields/half-initp25-d40-idealFinal-fields-0.bz2"
compfnFakefdm <- "../fields/diffusionBindhalf0-d40-fakefdm"
speedupLabel <- "Speed multiplier"
calibLabel <- ""

print(
      mapply(paramPlot,
             plotname = c("calib-realfdm", "calib-fakefdm"),
             comparisonFilename = c(compfnRealfdm, compfnFakefdm),
             MoreArgs = list(
                             savefilename = savefilename,
                             filenames = filenames,
                             paramMain = speedup, 
                             paramBy = type,
                             xLabel = speedupLabel,
                             byLabel = calibLabel,
                             height = 2.7,
                             order = c("Full calib.", "Grad. only", "Conc. only", "No calib.")
             )
      )
)
