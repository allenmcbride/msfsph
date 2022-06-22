June 22, 2022

This directory contains scripts for launching experiments to test the effects of changing the simulator's parameters as well as for
visualizing the results. 
The key scripts are `paramTest.py` for launching and `paramPlot.r` for visualization. `paramcheck.r` and `avgdist.r` contain some helper functions for the former scripts. 
The rest of the Python scripts (e.g., `paramTestHalfDiff.py`) launch experiments for particular simulation parameters. 
The rest of the R scripts (e.g., `plotHalfDiff.r`) visualize the results. 
Each Python script has a corresponding R script, more or less, though sometimes one visualization will apply to multiple experiments (e.g., MTTF). 
Also, `plotEdge.r` is a stand-alone script to visualize the results of an edge test hard-coded into the simulator.
