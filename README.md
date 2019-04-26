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
The LoDiS package is a 0D classical molecular dynamics software designed to simulate processes for finite-size systems
between 10-10000 atoms. Incorporated tools allow for investigations into growth, coalescence, quenching, phase transition, 
canonical esemble (NVT), microcanonical esemble (NVE) and metadynamics with the choice of a vacuum or ligand environment, and the presence of
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

Open the Makefile and select the correct version of fortran (listed under gfortran/ifort) and its libraries on your local computer by removing the necessary '#' symbols.

Compile all the .f90 files by running the Makefile:
```
cd LoDiS/LODIS_GIT
make -f Makefile
```

## Usage
Copy the **input.in** file into a chosen output directory.
The input variables for an MD simulation, including which type of process is to be simulated, are determined by modifying the **input.in** file. Background
information regarding each procedure as well tutorials on inputs and outputs are presented in the manual.
Supporting descriptions about the input parameters can be found in the [LoDiS Documentation](https://github.com/kcl-tscm/LoDiS/wiki/LoDiS-Documentation).

Additional files required for simulation include:  
* .xyz file with the initial nanocluster atom positions 
* .pot file with the potential parameters
* .pot file with the MgO substrate parameters (only when MgO substrate is present)
* .xyz file with second cluster atom positions (only during coalescence)


List the files to be read by the sofware in the **input.in** with either the relative path or the absolute path.
Example files for reference are provided in the **example_input_files** directory.

To run the simulation, use the following command line in the terminal:
```
PATH/LODIS_all <input.in> output.out
```
/PATH is either the relative from the **input.in** file directory (i.e the current directory) to the LoDiS_GIT directory, or the absolute path.
Upon initializing the simulation, LoDiS will run in the background till completion.

## Outputs
Running a simulation regardless of the procedure will output the following six files:
* energy.out - the caloric data over the course of the simulation
* movie.xyz - the trajectory data.
* error.out - a file listing any errors that may have occurred over the course of the run period or may have halted the process altogether
* definition.out - a binary file
* pr.out - the final positions of the components
* output.out - a file containing a runthorugh of the process, inluding parameters and errors.

An additional 4 file types depending on the process can be generated. The first two, '**meta.out**' and '**metahistory.out**', are produced 
after a metadynamics run whilst undergoing a NVE/NVT/quenching will give a '**velocity.out**' file. The last file type '**coalescing.out**' is outputted
a coalescence.   


## Example
A rapid 100K/ns melting simulation of an Ag 147 atom icosahedron is readily available with the use of the provided files in the **example_input_files** directory.

Edit the paths to the for Ag147.xyz and Ag_Ag.pot files before running the simulation:
```
filepos      = '~/Documents/LoDiS/input_example_files/Ag147.xyz',             ! Initial atom positions file, ONLY .xyz format
  
filepot      = '~/Documents/LoDiS/input_example_files/Ag_Ag.pot',             ! Potential parameters file, ONLY .pot format
```

NOTE: neither Ag_309.xyz nor Ag_Ag.MgO.pot are used in this particular run.

## The LoDiS Team
* Francesca Baletto (francesca.baletto@kcl.ac.uk)
* Raphael Pinto-Miles (raphpmx@gmail.com)
* Kevin Rossi (k1992@hotmail.it)
* Vagner Rigo (vagnerrigo@gmail.com)
* Riccardo Ferrando (ferrando@fisicaunige.it)
* Christine Mottet (mottet@cinam.univ-mrs.fr)
* Henrik Moerkved (henrik.morkved@kcl.ac.uk)
* Wei Zhao (wei.1.zhao@kcl.ac.uk)
* Robert Jones (robert.m.jones@kcl.ac.uk)


## References
[1] V. Rosato, M. Guillope and B. Legrand, *Philosophical Magazine A* **59**, 321 (1989)

[2] J. E. Jones, *Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences* **106**, 463 (1924)

[3] L. A. Girifalco, *The Journal of Physical Chemistry* **95**, 5370 (1991)

[4] K. Rossi, *Journal of Physics: Condensed Matter* **29**, 145402 (2017)

[5] K. Rossi and F. Baletto, *Physical Chemistry Chemical Physics* **19**, 11057 (2017)

[6] L. Verlet, *Phys. Rev.* **159**, 98 (1967)

[7] H.C.Andersen,*The Journal of Chemical Physics* **72**, 2384 (1980)

[8] F. Cleri and V. Rosato, *Phys. Rev. B* **48**, 22 (1993)

[9] F. Baletto, R. Ferrando, A. Fortunelli, F. Montalenti, and C. Mottet, *The Journal of Chemical Physics* **116**, 3856 (2002)

[10] J. E. Jones and S. Chapman, *Proceedings of the Royal Society of London. Series A, Containing Papers of a Mathematical and Physical Character* **106**, 463 (1924)

[11] R. Cortes-Huerto, J. Goniakowski, and C. Noguera, *The Journal of Chemical Physics* **138**, 244706 (2013)

[12] I. Atanasov, G. Barcaro, F. Negreiros, A. Fortunelli, and R. Johnston, *The Journal of chemical Physics* **138**, 224703 (2013).

[13] W. Vervisch, C. Mottet, and J. Goniakowski, *Phys. Rev. B* **65**, 245411 (2002)

[14] F. Baletto, C. Mottet, and R. Ferrando, *Phys. Rev. Lett.* **84**, 5544 (2000)

[15] F. Baletto, *Journal of Physics: Condensed Matter* **31**, 113001 (2019)

[16] A. Laio and M. Parrinello, *Proceedings of the National Academy of Sciences* **99**, 12562 (2002)

[17] G. A. Tribello, J. Cuny, H. Eshet, and M. Parrinello, *The Journal of Chemical Physics* **135**, 114109 (2011)

[18] G. Santarossa, A. Vargas, M. Iannuzzi, and A. Baiker, *Phys. Rev. B* **81**, 174205 (2010)
