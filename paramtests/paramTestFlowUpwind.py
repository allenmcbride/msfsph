#!/usr/bin/env python3

import paramTest

program = "ADVECTION"

params = [
        paramTest.Param(name = 'swimming', levels = ["false", "true"]),
        ]

otherargs = [
        "--initdens", "0.1", 
        "--diam", "80",
        "--tlim", "3600", 
        "--hperdiam", "6",
        "--nonewcalib",
        ]
replicates = 10
fieldName = "agentField-0"

paramTest.run("agents-" + program, params, otherargs, replicates, program, extraNameString = "adv80", dryRun = False, oneShot = True)
