#!/usr/bin/env python3

import paramTest

params = [
        paramTest.Param(name = 'diffusionBindhalf0-d40', levels = None),
        ]

otherargs = [
        "--initdens", "0.25", 
        "--tlim", "28800", 
        "--bindhalf", "0",
        "--nonewcalib",
        ]
replicates = 60 

paramTest.run("agents-DIFFUSION", params, otherargs, replicates, dryRun = False, oneShot = True)
