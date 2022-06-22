#!/usr/bin/env python3

#This one is redundant with paramTestGrowthLeakRate.py

import paramTest

params = [
        paramTest.Param(name = f"deriv-true", levels = None),
        ]

otherargs = [
        "--initdens", "0.1", 
        "--tlim", "64800",
        "--half", "2,2,4",
        "--bindhalf", "1",
        "--xphys", "6.0",
        "--nonewcalib",
        ]
replicates = 5

paramTest.run("agents-GROWTH-derivCtrl", params, otherargs, replicates)
