# Low Dimensional System Molecular Dynamics (LoDiS)

## CONTENTS
* [Background](##background)
* [Documentation](##documentation)
* [Installation](##installation)
* [Usage](##usage)
* [Outputs](##outputs)
* [Example](##example)

## Background
The LoDiS package 

## Documentation

## Installation
Clone the repository into a local directory:
'''
https://github.com/kcl-tscm/LoDiS.git
'''

Open and modify the Makefile to run the correct version of fortran on you local computer.
Compile all the .f90 files by running the Makefile:
'''
cd LoDis/LODIS_GIT
make -f Makefile
'''

## Usage
Inputs for an MD simulation including the type of process are determined by modifying the **input.in** file.
Supporting annotations within the file provide information on eath of the input variables. To run the simulation
input the command line into terminal:
'''
./PATH/LODIS_all <input.in> output.out
'''
**/PATH/ is the relative path from the input.in file to the LoDiS directory.**
NOTE: Upon initializing the simulation, LoDiS will run in the background till completion.

## Example

