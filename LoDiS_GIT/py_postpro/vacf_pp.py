import os
import sys
import cmath as cm
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


def acf_ff(tstep, atom_tot):
    # Read VACF file fort.56
    if os.path.isfile('fort.56'):
        vfile = open('fort.56', 'r+')
    else:
        #print('VACF file *fort.56* not detected within this directory')     #Used in gcn function
        print('Halting VACF analysis')
        return 'none'

    vread = vfile.readlines()
    vread = [i.split() for i in vread]
    acf_t = [float(l[1]) for l in vread]                      #VACF values at each timestep
    print('Finished reading VACF file')

    # Units conversion for VACF
    npas = len(acf_t) +1                      #total number of steps
    del_f = 1/(npas*tstep)                    #Frequency increment

    for m in range(npas -1):
        acf_t[m] = (acf_t[m]*npas*acf_t[0])/((npas - m))        #units = [m**2/s**2]

    #Calculate the VACF in frequency domain
    """
    f = m/(npas*tstep), m = 0,1,2,....npas-1
    t = p*tstep,       p = 1,2, ... npas-1
    ft = mp/npas
    """
    print('Carrying out Fourier transform to obtain vibrational spectral density (VDOS)')
    f_i = np.fft.fftfreq(npas, tstep)                  #Obtain frequencies from the discrete timesteps
    acf_f = 2*np.abs(np.fft.fft(acf_t))                #Fast Fourier Transform of the VACF, double due to imaginary amplitude

    #Plot f vs I
    g = plt.figure(1)
    h = 4.1357*10**(-15)                          # eV * s
    plot_f_i = f_i[:int((npas)/2)]/(10**12)           #Hz --> 10 ^12Hz
    plot_mev = h*f_i[:int((npas)/2)]*(10**3)           #Units in meV
    plot_acf_f = acf_f[:int((npas)/2)]
    fig = plt.figure()
    plt.plot(plot_mev,plot_acf_f )  #plot only in the positive domain [:n/2]
    #plt.xlabel('Frequency (10^12 Hz)')
    plt.xlabel('Phonon Energy (eV)')
    plt.ylabel('Intensity (a.u.)')
    #plt.xlim(0,max(plot_f_i))
    plt.ylim(min(plot_acf_f), max(plot_acf_f))
    plt.show()
    vfile.close()



