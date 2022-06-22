#!/usr/bin/env python3

import paramTest

program = "DIFFUSION_TRANSLATION"
params = [
        paramTest.Param(name = 'bindrate', levels = [2**i for i in range(-3, 9)]),
        ]

otherargs = [
        "--initdens", "0.25", 
        "--tlim", "28800", 
        "--nonewcalib",
        ]
replicates = 5

paramTest.run("agents-" + program, params, otherargs, replicates, program, dryRun = False)
