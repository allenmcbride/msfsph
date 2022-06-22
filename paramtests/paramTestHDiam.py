#!/usr/bin/env python3

import paramTest
import functools
import paramcheck

params = [
        paramTest.Param(name = 'half', levels = [2**i for i in range(-2, 11)]),
        paramTest.Param(name = 'diam', levels = [10 * 2**i for i in range(0, 2)])
        ]

otherargs = [
        "--initdens", "0.25", 
        "--tlim", "28800", 
        "--minxpix", "750",
        "--nonewcalib",
        ]
replicates = 10

def checker(paramSet):
    return paramcheck.parametersOkay(True, paramSet[1], 0.0, 100.0, paramSet[0])

paramTest.run("agents-DIFFUSION", params, otherargs, replicates, checkparams = checker, dryRun = False, oneShot = True)
