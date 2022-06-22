#!/usr/bin/env python3

import paramTest
import math

params = [
        paramTest.Param(name = 'motionnoise', levels = [2**i for i in range(-10, 6)]),
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
replicates = 5

paramTest.run("agents-GROWTH", params, otherargs, replicates, dryRun = False)
