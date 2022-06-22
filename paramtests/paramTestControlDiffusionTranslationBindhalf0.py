#!/usr/bin/env python3

import paramTest

program = "DIFFUSION_TRANSLATION"

params = [
        paramTest.Param(name = 'diffusion_translationBindhalf0-d40', levels = None),
        ]

otherargs = [
        "--initdens", "0.25", 
        "--tlim", "28800", 
        "--bindhalf", "0",
        "--nonewcalib",
        ]
replicates = 60 

paramTest.run("agents-" + program, params, otherargs, replicates, ".", dryRun = False, oneShot = True)
