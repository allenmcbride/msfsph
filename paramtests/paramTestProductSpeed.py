#!/usr/bin/env python3

import paramTest

program = "DIFFUSION_TRANSLATION"
params = [
        #paramTest.Param(name = 'half', levels = [2**i for i in range(-2, 11)]),
        paramTest.Param(name = 'half', levels = [2**i for i in range(1, 7)]),
        #paramTest.Param(name = 'half', levels = [2**i for i in range(3, 9)]),
        paramTest.Param(name = 'speedup', levels = [2**i for i in range(-3, 7)])
        ]

otherargs = [
        "--initdens", "0.25", 
        "--tlim", "28800", 
        "--smoothingh", "372.15416687891657557517",
        #"--smoothingh", "1488.61666751566630230067",
        "--maxspeed", "0.002", 
        "--bindhalf", "0",
        "--nonewcalib",
        ]
replicates = 5

paramTest.run("agents-" + program, params, otherargs, replicates, program, dryRun = False)
