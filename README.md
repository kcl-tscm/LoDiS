# Low Dimensional System Molecular Dynamics (LoDiS)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

![Lodis Logo](/images/lodislogo.png)

## Contents
* [Background](#background)
* [Installation](#installation)
* [Usage](#usage)
* [Outputs](#outputs)
* [Example](#example)
* [Contributors](#contributors)
* [References](#references)


## Background
The LoDiS package is a 0D classical molecular dynamics Fortran 90 - Python hybrid software designed to simulate processes for finite-size systems
between 10-10000 atoms. Incorporated tools allow for investigations into growth, coalescence, quenching, phase transition, 
canonical esemble (NVT), microcanonical esemble (NVE) and metadynamics with the choice of either a gas phase or ligand environment, and the addition of
a MgO substrate.  

Supported nanosystems include:
* Mono- and bi-metallic clusters (metal-metal interactions are modelled by the Rosato-Guillope-Legrande potential) [1]
* Noble gases (Lennard-Jones potential) [2]
* Carbon-based systems (Pacheco-Girifalco potential) [3]

For more general information and publications visit [Baletto group website](https://balettogroup.weebly.com/lodis.html)

## Installation
Clone the repository into a local directory:
```
https://github.com/kcl-tscm/LoDiS.git
```

Modify the Makefile in the LoDiS_GIT/base directory to run your local Fortran compiler and its libraries (openmp required).

Compile all the .f90 files by running the Makefile:
```
cd LoDiS/LODIS_GIT
make -f Makefile
```

## Usage
Before configuring the simulation, several files must be prepared beforehand:  
* .xyz file with the initial nanocluster atom positions 
* .pot file with the potential parameters
* .pot file with the MgO substrate parameters (only when MgO substrate is present)
* .xyz file with second cluster atom positions (only during coalescence)


LoDiS comes equipped with two input schemes: an interactive interface and an input.in file. The former is a more user friendly method for those unfamiliar
with the program, providing visual aids and error checks before simulation. In addition, post-simulation analysis unique to the scheme can be automatically run.
To run the interface and by extension LoDiS, use the following command on terminal:
```
python /PATH/LoDiS_GIT/py_interface/lodis_gui.py
```
/PATH is either the relative from the **input.in** file directory (i.e the current directory) to the LoDiS_GIT directory, or the absolute path.
After allocaring values to all the necessary parameters and confirming, an **input.in** file will be created and LoDiS will begin the simulation.

The second method involves copying the complete **input.in** file from the **/LoDiS_GIT/base** directory into the output directory and modifying the parameter
values within it. The .yxz and .pot files paths are included in the **input.in** and should be given as the relative path for the best result.
After altering the **input.in** file, begin the simulation by running the executable generated by compiling the Makefile, which should be called **LODIS_all**: 
```
/PATH/LoDiS_GIT/base/LODIS_all <input.in> output.out
```
Background information regarding each procedure as well tutorials on inputs and outputs are presented in the manual.
Supporting descriptions about the input parameters can be found in the [LoDiS Documentation](https://github.com/kcl-tscm/LoDiS/wiki/LoDiS-Documentation).

## Outputs
Running a simulation regardless of the procedure will output the following six files:
* energy.out - the caloric data over the course of the simulation
* movie.xyz - the trajectory data.
* error.out - a file listing any errors that may have occurred over the course of the run period or may have halted the process altogether
* definition.out - a binary file
* pr.out - the final positions of the components
* output.out - a file containing a runthorugh of the process, inluding parameters and errors.

Depending on the procedure, additional files may be generated such as **meta.out** for Metadynamics and **coalescing.out** for Coalescence.   


## Example
A rapid 100K/ns melting simulation of an Ag 147 atom icosahedron is readily available with the use of the provided files in the **example_input_files** directory.

Edit the paths to the for Ag147.xyz and Ag_Ag.pot files before running the simulation:
```
filepos      = '~/Documents/LoDiS/input_example_files/Ag147.xyz',             ! Initial atom positions file, ONLY .xyz format
  
filepot      = '~/Documents/LoDiS/input_example_files/Ag_Ag.pot',             ! Potential parameters file, ONLY .pot format
```

NOTE: neither Ag_309.xyz nor Ag_Ag.MgO.pot are used in this particular run.

## Contributors
* Francesca Baletto (francesca.baletto@kcl.ac.uk)
* Raphael Pinto-Miles (raphpmx@gmail.com)
* Kevin Rossi (k1992@hotmail.it)
* Vagner Rigo (vagnerrigo@gmail.com)
* Riccardo Ferrando (ferrando@fisicaunige.it)
* Christine Mottet (mottet@cinam.univ-mrs.fr)



## References
[1] V. Rosato, M. Guillope and B. Legrand, *Philosophical Magazine A* **59**, 321 (1989)

[2] J. E. Jones, *Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences* **106**, 463 (1924)

[3] L. A. Girifalco, *The Journal of Physical Chemistry* **95**, 5370 (1991)

[4] K. Rossi, *Journal of Physics: Condensed Matter* **29**, 145402 (2017)

[5] K. Rossi and F. Baletto, *Physical Chemistry Chemical Physics* **19**, 11057 (2017)
