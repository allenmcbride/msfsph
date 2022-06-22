#!/usr/bin/env python3

import paramTest

program = "ADVECTION"

paramsT = [ paramTest.Param(name = f"upwind-true", levels = None), ]

paramsF = [ paramTest.Param(name = f"upwind-false", levels = None), ]

otherargs = [
        "--initdens", "0.1", 
        "--tlim", "3600", 
        "--nonewcalib",
        ]
replicates = 5
fieldName = "agentField-0"

paramTest.run("agents-" + program, paramsT, otherargs, replicates, program, fieldName)
paramTest.run("agents-" + program + "-noUpwind", paramsF, otherargs, replicates, program, fieldName)
