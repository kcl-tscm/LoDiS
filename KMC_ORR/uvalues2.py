# These are the variables going into our code that apply to all the systems and are not dependent on the geometry or chemical ordering

inv_kb=11604.51911 # Inverse Boltzmann constant (K/eV)
h=4.14 # Planck constant [ev/fs]

# Duration
itt_max = 60000000 # number of Monte Carlo steps
t_max = 2000000 # total time for the simulation [fs]

# Ambient
temperature = 300 # [K] 
frequency = temperature/(h*inv_kb)*10**15 # [s^-1] equivalent to [10^-15 fs^-1] determined by transition state theory Arrhenius equation 

# Deposition
coverage = 0.25 # proportion of sites where O2 will be initially deposited on the surface 
maxdrop = 5000 # maximum number of (O2) particles that can be deposited

testno = 1

