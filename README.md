# Low Dimensional System Molecular Dynamics (LoDiS)

![LoDiS Logo](../../LoDiS_MD/documentation/Mini_Version_2016/lodis_logo.png)
Format: ![Alt Text](url)

## CONTENTS
* [Background](#background)
* [Documentation](#documentation)
* [Installation](#installation)
* [Usage](#usage)
* [Outputs](#outputs)
* [Example](#example)
* [References](#references)

## Background
The LoDiS package is a semi-empirical classical molecular dynamics software to investigate processes for finite-size systems
between 10-1000 atoms. Incorporated tools allow for investigations into growth/coalescence, quenching, phase transition, 
free energy and metadynamics with or without an MgO substrate. Supported nanosystems include:

* Mono- and bimetallic clusters (metal-metal interactions are modelled by the Rosato-Guillope-Legrande potential)
* Noble gases (Lennard-Jones potential)
* Carbon-based sytems (Pacheco-Girifalco potential)

## Installation
Clone the repository into a local directory:
```
https://github.com/kcl-tscm/LoDiS.git
```

Open and modify the Makefile to run the correct version of fortran (gfortran/ifort) and its libraries on your local computer.
Compile all the .f90 files by running the Makefile:
```
cd LoDiS/LODIS_GIT
make -f Makefile
```

## Usage
Copy the **input.in** file into a chosen output directory.
Inputs for an MD simulation including which type of process to be simulated are determined by modifying the **input.in** file.
Supporting annotations within the file provide information on each variable. To run the simulation, use the following command line in terminal:
```
./PATH/LODIS_all <input.in> output.out
```
/PATH is the relative path from the **input.in** file directory (i.e current directory) to the LoDiS_GIT directory.
Upon initializing the simulation, LoDiS will run in the background till completion

## Example
A rapid 100K/ns melting simulation of an Ag 147 atom iscosahedron is readily available with the use of the provived files in the **input_example_files** directory.
The chosen melting rate is determined by the variable relation:
```
deltat/(npas*tstep)
```

**Generally tstep shouldn't be changed unless required**

**filepos must be in an .xyz file format**

**filepot, mgo_pot and filepos2 must be .pot file format**

NOTE: filepos2 is not used in this particular run

Running the simulation will output the following six files:
* energy.out
* movie.xyz
* error.out
* definition.out
* pr.out
* output.out

The caloric and trajectory data are written into the energy.out and movie.xyz files respectively.

## References