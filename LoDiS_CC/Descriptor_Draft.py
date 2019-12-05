# This project is going to be a collection of all other codes for the
#Classification and characterisation of bi/monometallic structures
#The general jist is that we'll be reading in movie.xyz & energy.out
#Then producing statistical plots correlating descriptors


import time
import matplotlib.pyplot as plt
import numpy as np

import pickle
import sys
import os

#The following imports are the dependencies found in the \
#parent directory.

import Movie_Read as Read
import Kernels
from Kernels import KB_Dist as KB
import Distances as Dist


print("Welcome to this LoDiS post-processing scheme."
      "This script takes energy.out and movie.xyz files as arguments"
      "unless otherwise specified by name in the following input"
      "requests.")

#Below is the general scheme by which one takes the trajectory input.

Energy, Trajectory, Elements, filepath = Read.Manual_Input()


#A little bit of proof-of-concept code to plot the first
#few PDDFs using this scheme.

if __name__ == '__main__':

    Positions=[];Distances=[];PDDF=[]

    for x in range(5):
        Positions.append(np.column_stack((Elements[x],Trajectory[x])))
        Distances.append(Dist.Distance.Euc_Dist(Positions[x]))   #All possible pairings
        PDDF.append(Kernels.Kernels.KDE_Uniform(Distances[x],0.25))
        plt.plot(PDDF[x][0],PDDF[x][1])
        plt.show()
        
