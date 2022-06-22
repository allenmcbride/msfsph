#!/usr/bin/env python3

# This script should be run after compiling without correcting for
# calibrated value. (I could also consider testing with ADVECTION instead
# to get at the effect on Laplacian.)

import paramTest

program = "DIFFUSION_TRANSLATION"
params = [
        paramTest.Param(name = 'speedup', levels = [2**i for i in range(-3, 7)])
        ]

otherargs = [
        "--initdens", "0.25", 
        "--tlim", "28800", 
        "--maxspeed", "0.01", 
        "--bindhalf", "0",
        "--nonewcalib",
        ]
replicates = 5

paramTest.run("agents-" + program + "-noCalibVal", params, otherargs, replicates, program, extraNameString = "noCalibVal", dryRun = False)
