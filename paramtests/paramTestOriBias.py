#!/usr/bin/env python3

import paramTest

program = "DIFFUSION_TRANSLATION"
params = [
        paramTest.Param(name = 'oribias', levels = [(2**i) for i in range(-19, -1)]),
        ]

otherargs = [
        "--initdens", "0.25", 
        "--tlim", "28800", 
        "--nonewcalib",
        ]
replicates = 10

paramTest.run("agents-" + program, params, otherargs, replicates, program, dryRun = False)
