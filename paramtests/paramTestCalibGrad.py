#!/usr/bin/env python3

# This script should be run after compiling without correcting for
# calibrated gradient.

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

paramTest.run("agents-" + program + "-noCalibGrad", params, otherargs, replicates, program, extraNameString = "noCalibGrad", dryRun = False)
