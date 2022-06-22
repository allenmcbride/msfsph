#!/usr/bin/env python3

import paramTest

params = [
        paramTest.Param(name = 'diffusion-d40', levels = None),
        ]

otherargs = [
        "--initdens", "0.25", 
        "--tlim", "28800", 
        "--nonewcalib",
        ]
replicates = 20 

paramTest.run("agents-DIFFUSION", params, otherargs, replicates, dryRun = False)
