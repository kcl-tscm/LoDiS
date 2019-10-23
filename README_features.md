# New contents for Low Dimensional System Molecular Dynamics (LoDiS)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

![Lodis Logo](/images/lodislogo.png)

## Contents
* [Aim](#Aim)
* [Contributor](#contributor)
* [Description](#description)
* [Outputs](#outputs)
* [Example](#example)
* [References](#references)


## Aim
The LoDiS package is in a continuous evolution. Recent updates that might need to be tested are available here.

Currently we are working on coalescence, vibrational spectra, pressure, and nanogenomics

For more general information and publications visit [Baletto group website](https://www.balettogroup.org/lodis.html)

Common rules are the follwoing:
Before configuring the simulation, several files must be prepared beforehand:  
* .xyz file with the initial nanocluster atom positions 
* .pot file with the potential parameters
* .pot file with the MgO substrate parameters (only when MgO substrate is present)
* .xyz file with second cluster atom positions (only during coalescence)

Every new tool should be added as a logical variable deafult value as false that can be defined in the input **input.in** 

## Contributor
the person in charge of any changes shoudl be added in the specific part of the code and their name added in the **contributors** list.

## Description
Describe briefly what subroutines have been changed and a list of your changes. Also add a description of the new input and which variables have been introduced.

## Outputs
Provide a brief description of the new output. Please specify the units if needed.

## Example
Please provide an example of possible output

## References
Add any useful reference.
