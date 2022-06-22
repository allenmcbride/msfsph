#!/usr/bin/env python3

import paramTest

program = "DIFFUSION"

params = [
        #paramTest.Param(name = 'ecc', levels = [0, 0.8, 0.96]),
        paramTest.Param(name = 'ecc', levels = [0.96]),
        ]

otherargs = [
        "--initdens", "0.25", 
        "--mass", "0.0012566370614359175",
        "--tlim", "28800", 
        "--nonewcalib",
        ]
replicates = 20

paramTest.run("agents-" + program, params, otherargs, replicates, program, dryRun = False, oneShot = True)
