#!/usr/bin/env python3

import paramTest

program = "GROWTH"
params = [
        #paramTest.Param(name = 'leakhalf', levels = [3600 * 2**i for i in range(-2, 9)]),
        paramTest.Param(name = 'leakhalf', levels = [60 * 2**i for i in range(-2, 12)]),
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

paramTest.run("agents-" + program + "-derivctrl", params, otherargs, replicates, program)
