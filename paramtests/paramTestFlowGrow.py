#!/usr/bin/env python3

import paramTest

program = "GROWTH"

params = [
        paramTest.Param(name = 'swimming', levels = ["false", "true"]),
        ]

otherargs = [
        "--initdens", "0.1", 
        "--diam", "160",
        "--tlim", "64800", 
        "--hperdiam", "6",
        "--half", "16,16,32",
        "--diff", "100,100,100",
        "--xphys", "6.0",
        "--nonewcalib",
        ]
replicates = 10

paramTest.run("agents-" + program, params, otherargs, replicates, program, extraNameString = "grow", dryRun = False)
