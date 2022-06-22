#!/usr/bin/env python3

import paramTest

program = "GROWTH"
params = [
        paramTest.Param(name = f"targetCorr-false", levels = None),
        ]

otherargs = [
        "--initdens", "0.1",
        "--tlim", "64800",
        "--half", "2,2,4",
        "--diff", "100,100,100",
        "--bindhalf", "1",
        "--xphys", "6.0",
        "--nonewcalib",
        ]
replicates = 20

paramTest.run("agents-" + program + "-NOTARGCORR", params, otherargs, replicates, program)
