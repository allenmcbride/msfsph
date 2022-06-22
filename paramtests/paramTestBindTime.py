#!/usr/bin/env python3

import paramTest

program = "DIFFUSION_TRANSLATION"

params = [
        paramTest.Param(name = 'bindhalf', levels = [2**i for i in range(-2, 10)]),
        paramTest.Param(name = 'half', levels = [2**i for i in range(1, 7)]),
        ]

otherargs = [
        "--initdens", "0.25", 
        "--tlim", "28800", 
        "--smoothingh", "372.15416687891657557517",
        "--nonewcalib",
        ]
replicates = 5

paramTest.run("agents-" + program, params, otherargs, replicates, program, dryRun = False)
