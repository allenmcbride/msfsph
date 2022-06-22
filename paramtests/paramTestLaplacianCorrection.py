#!/usr/bin/env python3

import paramTest

paramsT = [ paramTest.Param(name = f"laplcorr-true", levels = None), ]

paramsF = [ paramTest.Param(name = f"laplcorr-false", levels = None), ]

otherargs = [
        "--initdens", "0.1", 
        "--tlim", "3600", 
        "--nonewcalib",
        ]
replicates = 5
fieldName = "agentField-0"

paramTest.run("agents-PARABOLOID", paramsT, otherargs, replicates, fieldName)
paramTest.run("agents-PARABOLOID-noLaplCorr", paramsF, otherargs, replicates, fieldName)
