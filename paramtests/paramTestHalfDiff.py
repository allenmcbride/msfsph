#!/usr/bin/env python3

import paramTest
import functools
import paramcheck

params = [
        paramTest.Param(name = 'half', levels = [2**i for i in range(-2, 11)]),
        paramTest.Param(name = 'diff', levels = [25, 50, 100, 200, 400, 800])
        ]

otherargs = [
        "--initdens", "0.25", 
        "--tlim", "28800", 
        "--nonewcalib",
        ]
replicates = 5

def checker(paramSet):
    return paramcheck.parametersOkay(True, 40, 0.0, paramSet[1], paramSet[0])

paramTest.run("agents-DIFFUSION", params, otherargs, replicates, checkparams = checker)
