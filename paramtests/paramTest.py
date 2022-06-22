# paramTest.py
# Allen McBride
# June 22, 2022
#
# run() looks at available results files to see how many replicates are
# still needed, and launches simulations as needed until enough results
# are found. Because results are searched after each new replicate,
# this script can be run on several machines at once with a shared
# filesystem. the oneShot option skips checking for results and simply
# runs the desired number of replicates; it is useful when the number of
# available machines is much larger than the number of replicates desired.

import filelock
import glob
import itertools
import pathlib
import random
import subprocess
import sys
import typing

class Param(typing.NamedTuple):
    name: str
    levels: list

def countReplicates(baseName, fieldName, subdir):
    nRep = len(glob.glob("fields/" + subdir + "/" + baseName + "-*" + "-" + fieldName + ".bz2"))
    print("Number of completed replicates found:", nRep)
    return nRep

def nextNumber(baseName):
    repsdir = pathlib.Path("replicates")
    repsdir.mkdir(exist_ok = True)
    filepath = repsdir / baseName
    lock = filelock.FileLock(str(filepath) + ".lock")
    with lock:
        if filepath.exists():
            with filepath.open() as file:
                numbersUsed = [int(line) for line in file]
            nextNumber = max(numbersUsed) + 1
        else:
            nextNumber = 1
        with filepath.open(mode = "a") as file:
            file.write(f"{nextNumber}\n")
    print("Next replicate label:", nextNumber)
    return nextNumber

def nameAndArgs(paramNames, paramSet, extraNameString):
    paramArgs = [] 
    picdir = ""
    for pName, pLevel in zip(paramNames, paramSet):
        paramArgs.extend([f"--{pName}={str(pLevel)}"])
        picdir += f"{pName}-{pLevel}-"
    if extraNameString is not None:
        picdir += extraNameString + "-"
    return picdir, paramArgs

def runHelper(cmdName, picdir, allOtherArgs, subdir, dryRun = False):
    logdir = pathlib.Path("logs")
    logdir.mkdir(exist_ok = True)
    logfilename = f"{picdir}.log"
    filepath = logdir / logfilename
    with filepath.open(mode = "w") as logfile:
        args = [
            "nice", "-n", "19",
            "time",
            "./" + cmdName, 
            "--picdir", picdir, 
            "--subdir", subdir,
            ] + allOtherArgs
        print("To run:", " ".join(args), flush = True)
        if not dryRun:
            subprocess.run(args, stdout = logfile, stderr = subprocess.STDOUT, check = True)

def run(cmdName, params, otherargs, nReplicatesNeeded, subdir = '.', fieldName = "gaussianDens", extraNameString = None, checkparams = lambda paramSet: True, dryRun = False, oneShot = False):
    if not isinstance(cmdName, str):
        sys.exit("First argument to run() should be name of a command to run")

    if params[0].levels is None:
        if len(params) > 1:
            sys.exit("Error: For control runs, first argument should be a list of one Param object with levels member set to None.")
        baseName = params[0].name
        if oneShot:
            picdir = f"{baseName}-{nextNumber(baseName)}"
            runHelper(cmdName, picdir, otherargs, subdir, dryRun)
        elif dryRun:
            runHelper(cmdName, baseName, otherargs, subdir, dryRun)
        else:
            moreNeeded = True
            while(moreNeeded):
                nReplicatesFound = countReplicates(baseName, fieldName, subdir)
                if nReplicatesFound >= nReplicatesNeeded:
                    print("Found all needed replicates for:", baseName)
                    moreNeeded = False
                else:
                    print("Still need", nReplicatesNeeded - nReplicatesFound, "replicates for:", baseName)
                    picdir = f"{baseName}-{nextNumber(baseName)}"
                    runHelper(cmdName, picdir, otherargs, subdir, dryRun)

    else:
        paramSets = [p for p in itertools.product(*[p.levels for p in params]) if checkparams(p)]
        paramNames = [p.name for p in params]
        if oneShot:
            for paramSet in paramSets:
                picdir, paramArgs = nameAndArgs(paramNames, paramSet, extraNameString)
                baseName = picdir.rstrip("-")
                picdir += str(nextNumber(baseName))
                runHelper(cmdName, picdir, paramArgs + otherargs, subdir, dryRun)
        elif dryRun:
            for paramSet in paramSets:
                picdir, paramArgs = nameAndArgs(paramNames, paramSet, extraNameString)
                runHelper(cmdName, picdir, paramArgs + otherargs, subdir, dryRun)
        else:
            while len(paramSets) > 0:
                print("Number of parameter sets remaining:", len(paramSets))
                paramSet = random.choice(paramSets)
                picdir, paramArgs = nameAndArgs(paramNames, paramSet, extraNameString)
                baseName = picdir.rstrip("-")
                nReplicatesFound = countReplicates(baseName, fieldName, subdir)
                if nReplicatesFound >= nReplicatesNeeded:
                    print("Found all needed replicates for:", paramSet)
                    paramSets.remove(paramSet)
                else:
                    print("Still need", nReplicatesNeeded - nReplicatesFound, "replicates for:", paramSet)
                    picdir += str(nextNumber(baseName))
                    runHelper(cmdName, picdir, paramArgs + otherargs, subdir, dryRun)
