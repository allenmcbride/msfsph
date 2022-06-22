#!/usr/bin/env python3

import paramTest
import subprocess

def run(binary, picdir, otherargs):
    paramTest.runHelper(binary, picdir, ["--ideal"] + otherargs, ".", False)

#run("agents-DIFFUSION", "half-initp25-d80", ["--tlim", "28800", "--diff", "100", "--diam", "80", "--initdens", "0.25"])
#run("agents-DIFFUSION", "half-initp25-d40", ["--tlim", "28800", "--diff", "100", "--diam", "40", "--initdens", "0.25"])
#run("agents-DIFFUSION", "half-initp25-d20", ["--tlim", "28800", "--diff", "100", "--diam", "20", "--initdens", "0.25"])
#run("agents-DIFFUSION", "half-initp25-d10", ["--tlim", "28800", "--diff", "100", "--diam", "10", "--initdens", "0.25"])
#run("agents-DIFFUSION", "half-initp25-d40-ecc0.8", ["--tlim", "28800", "--diff", "100", "--diam", "40", "--ecc", "0.8", "--initdens", "0.25"])
#run("agents-DIFFUSION", "half-initp25-d40-ecc0.96", ["--tlim", "28800", "--diff", "100", "--diam", "40", "--ecc", "0.96", "--initdens", "0.25"])
#run("agents-DIFFUSION", "half-initp25-mass-ecc0.8", ["--tlim", "28800", "--diff", "100", "--mass", "0.0012566370614359175", "--ecc", "0.8", "--initdens", "0.25"])
#run("agents-DIFFUSION", "half-initp25-mass-ecc0.96", ["--tlim", "28800", "--diff", "100", "--mass", "0.0012566370614359175", "--ecc", "0.96", "--initdens", "0.25"])

#run("agents-PARABOLOID", "laplcorr", ["--tlim", "3600", "--diff", "1000", "--diam", "40", "--initdens", "0.1"])

#run("agents-ADVECTION", "upwind-d40", ["--diam", "40"])

#run("agents-GROWTH", "growth-d40", ["--tlim", "64800", "--diam", "40", "--initdens", "0.1", "--xphys", "6.0"])
run("agents-GROWTH", "growth-d160", ["--tlim", "64800", "--diam", "160", "--initdens", "0.1", "--xphys", "6.0", "--hperdiam", "6"])
