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
The LoDiS package is a semi-empirical classical molecular dynamics software to investigate processes for finite-size systems
between 10-10000 atoms. Incorporated tools allow for investigations into growth/coalescence, quenching, phase transition, 
canonical NVT and metadynamics with or without the presence of an MgO substrate. 

Supported nanosystems include:
* Mono- and bi-metallic clusters (metal-metal interactions are modelled by the Rosato-Guillope-Legrand potential) [1]
* Noble gases (Lennard-Jones potential) [2]
* Carbon-based systems (Pacheco-Girifalco potential) [3]

LoDiS allows for two choices of coordination number calculations; either by Fermi distribution formalisation [3] or 
via analytical and polynomial formalisation [4].

For more general information and publications visit [Baletto group website](https://balettogroup.weebly.com/lodis.html)

## Installation
Clone the repository into a local directory:
```
https://github.com/kcl-tscm/LoDiS.git
```

Open the Makefile and select the correct version of fortran (listed under gfortran/ifort) and its libraries on your local computer by removing the necessary '#' symbols.

Compile all the .f90 files by running the Makefile:
```
cd LoDiS/LODIS_GIT
make -f Makefile
```

## Usage
Copy the **input.in** file into a chosen output directory.
The input variables for an MD simulation, including which type of process is to be simulated, are determined by modifying the **input.in** file.
Supporting information about the input parameters can be found in the [LoDiS Documentation](https://github.com/kcl-tscm/LoDiS/wiki/LoDiS-Documentation).

Other files required for simulation include:  
* .xyz file with the initial nanocluster atom positions 
* .pot file with the potential parameters
* .pot file with the MgO substrate parameters (only when MgO substrate is present)
* .xyz file with coalescence cluster positions (only during coalescence/growth)

List the files to be read by the software in the **input.in** with either the relative path or ideally the absolute path.
Example files for reference are provided in the **example_input_files** directory.

To run the simulation, use the following command line in terminal:
```
./PATH/LODIS_all <input.in> output.out
```
/PATH is the relative path from the **input.in** file directory (i.e the current directory) to the LoDiS_GIT directory.
Upon initializing the simulation, LoDiS will run in the background till completion.

## Example
A rapid 100K/ns melting simulation of an Ag 147 atom icosahedron is readily available with the use of the provided files in the **example_input_files** directory.

Edit the paths to the for Ag147.xyz and Ag_Ag.pot files before running the simulation:
```
filepos      = '~/Documents/LoDiS/input_example_files/Ag147.xyz',             ! Initial atom positions file, ONLY .xyz format
  
filepot      = '~/Documents/LoDiS/input_example_files/Ag_Ag.pot',             ! Potential parameters file, ONLY .pot format
```

NOTE: neither PT38TO.xyz nor Ag_Ag.MgO.pot are used in this particular run.

Running the simulation will output the following six files:
* energy.out
* movie.xyz
* error.out
* definition.out
* pr.out
* output.out

The caloric and trajectory data are written into the energy.out and movie.xyz files respectively.

## Contributors
* Francesca Baletto (francesca.baletto@kcl.ac.uk)
* Raphael Pinto-Miles (raphpmx@gmail.com)
* Kevin Rossi (k1992@hotmail.it)

## References
[1] V. Rosato, M. Guillope and B. Legrand, *Philosophical Magazine A* **59**, 321 (1989)

[2] J. E. Jones, *Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences* **106**, 463 (1924)

[3] L. A. Girifalco, *The Journal of Physical Chemistry* **95**, 5370 (1991)

[4] K. Rossi, *Journal of Physics: Condensed Matter* **29**, 145402 (2017)

[5] K. Rossi and F. Baletto, *Physical Chemistry Chemical Physics* **19**, 11057 (2017)
