#!/usr/bin/env python3

import paramTest

program = "DIFFUSION_TRANSLATION"

params = [
        paramTest.Param(name = 'diffusion_translation-d40', levels = None),
        ]

otherargs = [
        "--initdens", "0.25", 
        "--tlim", "28800", 
        "--nonewcalib",
        ]
replicates = 20 

paramTest.run("agents-" + program, params, otherargs, replicates, ".", dryRun = False, oneShot = True)
