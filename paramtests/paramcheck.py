# paramcheck.py
# Allen McBride
# June 22, 2022
#
# The parametersOkay() function is just to filter out parameter
# combinations that would take too long to calibrate for.

import math

def decay(half):
    return math.log(2.0) / (half * 60.0)

def smoothingh(diff, half):
    return math.sqrt(diff / decay(half))

def parametersOkay(passthrough, diam, ecc, diff, half):
    if smoothingh(diff, half) / diam < 1:
        return False
    const = 0.0
    if passthrough:
        if ecc == 0.0:
            timeConst = 32
        elif ecc == 0.8:
            timeConst = 10
        elif ecc == 0.96:
            timeConst = 10
    else:
        if ecc == 0.0:
            timeConst = 8
        elif ecc == 0.8:
            timeConst = 5
        elif ecc == 0.96:
            timeConst = 5
    timeOkay = smoothingh(diff, half) / diam < timeConst
    return timeOkay

def parametersOkayC1(passthrough, diam, ecc, diff, half):
    if smoothingh(diff, half) / diam < 1:
        return False
    const = 0.0
    if passthrough:
        timeConst = 0
    else:
        if ecc == 0.0:
            timeConst = 0
        elif ecc == 0.8:
            timeConst = 8
        elif ecc == 0.96:
            timeConst = 8
    timeOkay = smoothingh(diff, half) / diam < timeConst
    return timeOkay
