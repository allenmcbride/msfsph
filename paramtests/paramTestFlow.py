#!/usr/bin/env python3

import paramTest

program = "DIFFUSION"

params = [
        paramTest.Param(name = 'swimming', levels = ["false", "true"]),
        ]

otherargs = [
        "--initdens", "0.25", 
        "--diam", "160",
        "--tlim", "28800", 
        "--hperdiam", "6",
        "--nonewcalib",
        ]
replicates = 10

paramTest.run("agents-" + program, params, otherargs, replicates, program, extraNameString = "diff", dryRun = False)
