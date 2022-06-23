June 22, 2022

**msfsph** is a simulator supporting the research for my Ph.D. dissertation, [Toward Scalable Morphogenetic Engineering: Natural Computing in SPH Swarm Control](https://trace.tennessee.edu/utk_graddiss/7111/). It simulates agents moving in an aqueous medium, communicating through the secretion of morphogen-like chemicals that diffuse and degrade in their shared environment. The code for these agents is compiled from the Morphgen language, designed by Dr. Bruce J. MacLennan in part to describe artificial morphogenetic processes in the continuum limit using partial differential equations.

msfsph depends on [SDL2](https://www.libsdl.org/), [Eigen](https://eigen.tuxfamily.org/index.php), and [cxxopts](https://github.com/jarro2783/cxxopts). cxxopts is a header-only library, and the cxxopts.hpp file should be placed in the src directory of this project.
