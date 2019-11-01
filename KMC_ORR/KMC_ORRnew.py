import time
import sys
import numpy as np
# from reaclass import * # defines classes for reactions as well as move/event operations
from math import * 
import random
from uvalues2 import *
import matplotlib.pyplot as plt
#from gcn_bridge_onecluster import *
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
with open("inputgcn.dat", "r") as p:
	nl = int(p.readline())
	nm = int(p.readline())
	nh = int(p.readline())
	ccluster=int(p.readline())
        surfno=int(p.readline())
# Choice of architecture
# CONi@I == 0 , COSa == 1, TONi@I == 2 , TOSa == 3
i = 0

print(int(round(float(surfno)*nl/float(nl+nm+nh))), int(round(nm*float(surfno)/float(nl+nm+nh))), int(round(nh*float(surfno)*nl/float(nl+nm+nh))), surfno)
# collects various input conditions

if i==0:
    from initvalues8CONi import * 
if i==1:
    from initvalues8COSa import * 
if i==2:
    from initvalues8TONi import * 
if i==3:
    from initvalues8TOSa import * 


class Move(object):
    '''Data structure to store the kinematics of particle movement'''
    def __init__(self,rate,in_state,fi_state,react,position,gcn):
        self.rate = rate
        self.in_state = in_state
        self.fi_state = fi_state
        self.react = react
        self.position = position
        self.gcn = gcn


# Initialization
random.seed()
time1 = time.time()
t_total = 0
inv_kb=11604.51911 # Inverse Boltzmann constant (K/eV)
h=4.14 # Planck constant [ev/fs]
x = 0

# Sites on the nanoparticle

ntot = nl + nm + nh
ncover = ntot*coverage # number of particles needed on the given domain to satisfy the coverage


lsites = np.zeros((nl), dtype=np.int) # object containing the low gcn sites to evolve
msites = np.zeros((nm), dtype=np.int)
hsites = np.zeros((nh), dtype=np.int)

at_total = 0 # counter for each different domain
aritt = 0
als = [0,0,0,0,0,0,0,0] #state counter for each domain, O2 will be deposited first but not counted in the simulation
areact = np.zeros(48,dtype=np.int)
awater = 0 # water counter for different domains
at_total2 = 0
tonar = 0
tonar2 = 0
tonaracs = 0
tonarals = 0
tonaracs2 = 0
tonarals2 = 0

typco1 = 0
typco2 = 0
typco5 = 0
typco6 = 0
typco16 = 0

typcoin1 = 0
typcoin2 = 0
typcoin5 = 0
typcoin6 = 0
typcoin16 = 0

typcobindl = 0
typcobindm = 0
typcobindh = 0

typcodiffO2ml = 0
typcodiffO2lm = 0
typcodiffO2hm = 0
typcodiffO2mh = 0
typcodiffOOHhm = 0
typcodiffOOHmh = 0

dtypco1 = 0
dtypco2 = 0
dtypco5 = 0
dtypco6 = 0
dtypco16 = 0

dtypcoin1 = 0
dtypcoin2 = 0
dtypcoin5 = 0
dtypcoin6 = 0
dtypcoin16 = 0

dtypcobindl = 0
dtypcobindm = 0
dtypcobindh = 0

dtypcodiffO2ml = 0
dtypcodiffO2lm = 0
dtypcodiffO2hm = 0
dtypcodiffO2mh = 0
dtypcodiffOOHhm = 0
dtypcodiffOOHmh = 0


# structure for the energies: array of arrays of arrays, smallest arrays are by activation and reaction for low and high, next is by reaction number, finally by domain; [domains]->[[reactions]]->[[[energies]]]
E_all = [[E1lowact,E1lowreact,E1highact,E1highreact],[E2lowact,E2lowreact,E2highact,E2highreact],[E5lowact,E5lowreact,E5highact,E5highreact],[E6lowact,E6lowreact,E6highact,E6highreact],[E16lowact,E16lowreact,E16highact,E16highreact]]

#name = desg+"_collated.dat"
#name3 = "all_collated.dat"
#name2 = desg+"_test"+str(testno)+"_KMC_states"str(ccluster)++".dat"
name4 = desg+"_test"+str(testno)+"_KMC_SUMMARY"+str(ccluster)+".dat"

#col = open(name, 'a')
#allcol = open(name3, 'a')
#g = open(name2, 'w') # file for the states of each site (shows the reactants in each site in the array)
summ = open(str(dir_path)+ "/Output/"+name4, 'w')

lgco = open(str(dir_path)+"/Output/"+desg+"_test"+str(testno)+"_lgcnstatco"+str(ccluster)+".dat",'w')
mgco = open(str(dir_path)+"/Output/"+desg+"_test"+str(testno)+"_mgcnstatco"+str(ccluster)+".dat",'w')
hgco = open(str(dir_path)+"/Output/"+desg+"_test"+str(testno)+"_hgcnstatco"+str(ccluster)+".dat",'w')
agco = open(str(dir_path)+"/Output/"+desg+"_test"+str(testno)+"_agcnstatco"+str(ccluster)+".dat",'w')

lgco.write("#itt\tt\t\tO2g\tEmt\tO2a\t2O\tOH+O\tOOH\tH2O+O\t2OH\tH2O+OH\t2H2O \n")
mgco.write("#itt\tt\t\tO2g\tEmt\tO2a\t2O\tOH+O\tOOH\tH2O+O\t2OH\tH2O+OH\t2H2O \n")
hgco.write("#itt\tt\t\tO2g\tEmt\tO2a\t2O\tOH+O\tOOH\tH2O+O\t2OH\tH2O+OH\t2H2O \n")
agco.write("#itt\tt\t\tO2g\tEmt\tO2a\t2O\tOH+O\tOOH\tH2O+O\t2OH\tH2O+OH\t2H2O \n")

reactco = open(str(dir_path)+"/Output/"+desg+"_test"+str(testno)+"_reaction"+str(ccluster)+".dat",'w')
reactco.write("#itt\t1\t2\t5\t6\t16\t1in\t2in\t5in\t6in")
reactdiff = open(str(dir_path)+"/Output/"+desg+"_test"+str(testno)+"_binddiff"+str(ccluster)+".dat",'w')
reactdiff.write("#itt\tbinding(l,\tm,\th)\tO2 diff(lm,ml,mh,hm)\tOOH diff(mh,hm)")  

dtreactco = open(str(dir_path)+"/Output/"+desg+"_test"+str(testno)+"_dt_reactionc_test"+str(ccluster)+".dat",'w')
dtreactco.write("#itt   1  2  5  6  16  1in  2in  5in  6in")
dtreactdiff = open(str(dir_path)+"/Output/"+desg+"_test"+str(testno)+"_dt_reactiond_test"+str(ccluster)+".dat",'w')
dtreactdiff.write("#itt binding(l,m,h) O2 diff(lm,ml,mh,hm)   OOH diff(mh,hm)")  


E1 = E_all[0] # activation and reaction energies for first reaction for low and high GCN sites &c.
E2 = E_all[1]
E5 = E_all[2]
E6 = E_all[3]
E16 = E_all[4]
l=-1 # number corresponding to
ls = als # array for states
react = areact # array for reactions
t_total = 0 # initial 
water = 0
pt = 0
itt = 0
drop = 0 # total number of molecules deposited onto the nanoparticle for each domain
lsites = np.zeros((nl), dtype=np.int) # object containing the low gcn sites to evolve
msites = np.zeros((nm), dtype=np.int)
hsites = np.zeros((nh), dtype=np.int)
sitcover = 0 # counts how many of the sites are 'covered'

for no in range(nl): # initialising arrays to be empty when starting each domain
    lsites[no] = 0
for no in range(nm):
    msites[no] = 0
for no in range(nh):
    hsites[no] = 0
    
time_plot=[]
Emt=[]
O2=[]
TwoO=[]
OHplusO=[]
OOH=[]
H2OplusO=[]
TwoOH=[]
H2OplusOH=[]
water_plot=[]

for itt in range(itt_max + 1):
    '''Fill "moves_list", with Move obj. It is a list containing all the possible moves'''
    moves_list = []

    l=-1 # counter for number of LGCN free sites
    freels = []
    for osite in lsites:
        l+=1
        if osite == 0:
            freels.append(l)
    nfreels = len(freels)

    l=-1
    freems = []
    for osite in msites:
        l+=1
        if osite == 0:
            freems.append(l)
    nfreems = len(freems)

    l=-1
    freehs = []
    for osite in hsites:
        l+=1
        if osite == 0:
            freehs.append(l)
    nfreehs = len(freehs)

    l=-1 # we start with the position label as -1 so it can increase to 0 for the first site
    # 2.0) Reactions
    # Reactions for LGCN Sites:
    for osite in lsites:
        l+=1 # 'l' is made to increase for each subsequent site in the gcn space
        # 0) Reaction for Empty -> O2
        if osite == 0 and sitcover < ncover and drop < maxdrop and l0 == 'true': # requirements for step to take place (here that site is empty and particles remain to deposit)
            m = osite
            delta_total = en_bindl # energy barrier for this event
            rate = freq*exp(-delta_total*inv_kb/temperature) # rate for this event (product of prefactor with exp of ebarrier)
            n = 1
            rct = 0 # reaction label for given gcn (low)
            pos = l # position label for given gcn (low)
            gcn = 'low'
            moves_list.append(Move(rate,m,n,rct,pos,gcn)) # moves list updated with new event
        # 1) Reaction for O2 -> 2 O
        if osite == 1 and l1 == 'true':
            m = osite
            delta_total = E1[0] # the second element in the array E2 is the activation energy for a high GCN site for this reaction
            rate = freq1*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 2
            rct = 1
            pos = l
            gcn = 'low'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # 2) Reaction for 2 O -> OH + O
        if osite == 2 and l2 == 'true':
            m = osite
            delta_total = E2[0]
            rate = freq2*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 3
            rct = 2
            pos = l
            gcn = 'low'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # 4) Movement of O2 from LGCN to MGCN
        if osite == 1 and nfreems > 0 and l4 == 'true':
            m = osite
            delta_total = en_difflmO2
            rate = freq*exp(-delta_total*inv_kb/temperature)
            n = 0
            rct = 4
            pos = l
            gcn = 'low'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # 8) Reaction for 2 O -> O2
        if osite == 2 and l8 == 'true':
            m = osite
            Eact8 = E1[0] - E1[1] # the third element of array E3 is the reaction energy for a high GCN site for this reaction
            delta_total = Eact8
            rate = freqin1*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 1
            rct = 8
            pos = l
            gcn = 'low'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # 9) Reaction for OH + O -> 2 O
        if osite == 3 and l9 == 'true':
            m = osite
            Eact9 = E2[0] - E2[1]
            delta_total = Eact9
            rate = freqin2*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 2
            rct = 9
            pos = l
            gcn = 'low'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # 16) Reaction for OH + O -> H2O + O
        if osite == 3 and l16 == 'true':
            m = osite
            delta_total = E16[0]
            rate = freq16*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 5
            rct = 16
            pos = l
            gcn = 'low'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # 12) Reaction for OH + O -> 2OH
        if osite == 3 and l12 == 'true':
            m = osite
            delta_total = E2[0]
            rate = freq2*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 6
            rct = 12
            pos = l
            gcn = 'low'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # 13) Reaction for H2O + O -> H2O + OH
        if osite == 5 and l13 == 'true':
            m = osite
            delta_total = E2[0]
            rate = freq2*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 7
            rct = 13
            pos = l
            gcn = 'low'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # 23) Reaction for H2O + OH -> H2O + O
        if osite == 7 and l23 == 'true':
            m = osite
            delta_total = E2[0] - E2[1]
            rate = freq2*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 5
            rct = 23
            pos = l
            gcn = 'low'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # 14) Reaction for  2OH -> H2O + OH
        if osite == 6 and l14 == 'true':
            m = osite
            delta_total = E16[0]
            rate = freq16*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 7
            rct = 14
            pos = l
            gcn = 'low'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # 17) Reaction for H2O + OH -> 2H2O
        if osite == 7 and l17 == 'true':
            m = osite
            delta_total = E16[0]
            rate = freq16*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 0
            rct = 17
            pos = l
            gcn = 'low'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # 15) Reaction for 2OH -> OH + O
        if osite == 6 and l15 == 'true':
            m = osite
            delta_total = E2[0] - E2[1]
            rate = freqin2*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 3
            rct = 15
            pos = l
            gcn = 'low'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # 11) Reaction for OH + O -> OOH (LGCN to MGCN)
        if osite == 3 and nfreems > 0 and l11 == 'true':
            m = osite
            Eact11 = E6[2] - E6[3]
            delta_total = Eact11
            rate = freq6*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 0
            rct = 11
            pos = l
            gcn = 'low'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
    l=-1
    # Reactions for Medium GCN sites (near LGCN sites (one step)): 
    for osite in msites:
        l+=1
        # 7) Reaction for Empty -> O2
        if osite == 0  and sitcover < ncover and drop < maxdrop and l7 == 'true':
            m = osite
            delta_total = en_bindm
            rate = freq*exp(-delta_total*inv_kb/temperature)
            n = 1
            rct = 7
            pos = l
            gcn = 'med'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # 3) Movement of O2 from MGCN to LGCN
        if osite == 1 and nfreels > 0 and l3 == 'true':
            m = osite
            delta_total = en_diffmlO2
            rate = freq*exp(-delta_total*inv_kb/temperature)
            n = 0
            rct = 3
            pos = l
            gcn = 'med'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # 5) Reaction for O2 -> OOH
        if osite == 1 and l5 == 'true':
            m = osite
            delta_total = E5[2]
            rate = freq5*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 4
            rct = 5
            pos = l
            gcn = 'med'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # 6) Reaction for OOH -> OH + O (MGCN to LGCN)
        if osite == 4 and nfreels > 0 and l6 == 'true':
            m = osite
            delta_total = E6[2]
            rate = freq6*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 0
            rct = 6
            pos = l
            gcn = 'med'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # 10) Reaction for OOH -> O2
        if osite == 4 and l10 == 'true':
            m = osite
            Eact10 = E5[2] - E5[3]
            delta_total = Eact10
            rate = freq5*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 1
            rct = 10
            pos = l
            gcn = 'med'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # 24) Movement of O2 from MGCN to HGCN: 
        if osite == 1 and nfreehs > 0 and l24 == 'true':
            m = osite
            delta_total = en_diffmhO2
            rate = freq*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 0
            rct = 24
            pos = l
            gcn = 'med'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # 25) Movement of OOH form MGCN to HGCN: 
        if osite == 4 and nfreehs > 0 and l25 == 'true':
            m = osite
            delta_total = en_diff_OOH_temp
            rate = freq*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 0
            rct = 25
            pos = l
            gcn = 'med'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # NEW 26) Reaction for O2 -> 2O on med
        if osite == 1 and l26 == 'true':
            m = osite
            delta_total = E1[2]
            rate = freq1*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 2
            rct = 26
            pos = l
            gcn = 'med'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # NEW 27) Reaction for 2O -> OH + O on med
        if osite == 2 and l27 == 'true':
            m = osite
            delta_total = E2[2]
            rate = freq2*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 3
            rct = 27
            pos = l
            gcn = 'med'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # NEW 28) Reaction for 2O -> O2 on med
        if osite == 2 and l28 == 'true':
            m = osite
            delta_total = E1[2]-E1[3]
            rate = freqin1*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 1
            rct = 28
            pos = l
            gcn = 'med'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # NEW 29) Reaction for OH + O -> 2O on med
        if osite == 29:
            m = osite
            delta_total = E2[2]-E2[3]
            rate = freqin2*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 2
            rct = 29
            pos = l
            gcn = 'med'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # NEW 30) Reaction for OH + O -> H2O + O on med
        if osite == 3 and l30 == 'true':
            m = osite
            delta_total = E16[2]
            rate = freq16*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 5
            rct = 30
            pos = l
            gcn = 'med'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # NEW 31) Reaction for OH + O -> 2OH on med
        if osite == 3 and l31 == 'true':
            m = osite
            delta_total = E2[2]
            rate = freq2*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 6
            rct = 31
            pos = l
            gcn = 'med'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # NEW 32) Reaction for 2OH -> OH + O on med
        if osite == 6 and l32 == 'true':
            m = osite
            delta_total = E2[2]-E2[3]
            rate = freqin2*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 3
            rct = 32
            pos = l
            gcn = 'med'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # NEW 33) Reaction for H2O + O -> H2O + OH on med
        if osite == 5 and l33 == 'true':
            m = osite
            delta_total = E2[2]
            rate = freq2*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 7
            rct = 33
            pos = l
            gcn = 'med'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # NEW 34) Reaction for H2O + OH -> H2O + O on med
        if osite == 7 and l34 == 'true':
            m = osite
            delta_total = E2[2]-E2[3]
            rate = freqin2*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 5
            rct = 34
            pos = l
            gcn = 'med'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # NEW 35) Reaction for 2OH -> H2O + OH on med
        if osite == 6 and l35 == 'true':
            m = osite
            delta_total = E16[2]
            rate = freq16*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 7
            rct = 35
            pos = l
            gcn = 'med'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # NEW 36) Reaction for H2O + OH -> 2H2O on med
        if osite == 7 and l36 == 'true':
            m = osite
            delta_total = E16[2]
            rate = freq16*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 0
            rct = 36
            pos = l
            gcn = 'med'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
    l=-1
    # Reactions for High GCN (far from LGCN (more than one step)):
    for osite in hsites:
        l+=1
        # 18) Reaction for Empty -> O2
        if osite == 0  and sitcover < ncover and drop < maxdrop and l18 == 'true':
            m = osite
            delta_total = en_bindh
            rate = freq*exp(-delta_total*inv_kb/temperature)
            n = 1
            rct = 18
            pos = l
            gcn = 'high'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # 19) Reaction for O2 -> OOH
        if osite == 1 and l19 == 'true':
            m = osite
            delta_total = E5[2]
            rate = freq5*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 4
            rct = 19
            pos = l
            gcn = 'high'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # 20) Movement of O2 from HGCN to MGCN: (when counter 3 is too great it moves to mgcn)
        if osite == 1 and nfreems > 0 and l20 == 'true':
            m = osite
            delta_total = en_diffhmO2
            rate = freq*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 0
            rct = 20
            pos = l
            gcn = 'high'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # 21) Movement of OOH form HGCN to MGCN: (when counter 3 is too great it moves to mgcn)
        if osite == 4 and nfreems > 0 and l21 == 'true':
            m = osite
            delta_total = en_diff_OOH_temp - 0.7 # diffusion barrier is temporary value and -.7 is there to make the movement to mgcn more favourable
            rate = freq*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 0
            rct = 21
            pos = l
            gcn = 'high'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # 22) Reaction of OOH -> O2
        if osite == 4 and l22 == 'true':
            m = osite
            delta_total = E5[2] - E5[3]
            rate = freqin5*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 1
            rct = 22
            pos = l
            gcn = 'high'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # NEW 37) Reaction for O2 -> 2O on high
        if osite == 1 and l37 == 'true':
            m = osite
            delta_total = E1[2]
            rate = freq1*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 2
            rct = 37
            pos = l
            gcn = 'high'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # NEW 38) Reaction for 2O -> OH + O on high
        if osite == 2 and l38 == 'true':
            m = osite
            delta_total = E2[2]
            rate = freq2*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 3
            rct = 38
            pos = l
            gcn = 'high'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # NEW 39) Reaction for 2O -> O2 on high
        if osite == 2 and l39 == 'true':
            m = osite
            delta_total = E1[2]-E1[3]
            rate = freqin1*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 1
            rct = 39
            pos = l
            gcn = 'high'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # NEW 40) Reaction for OH + O -> 2O on high
        if osite == 3 and l40 == 'true':
            m = osite
            delta_total = E2[2]-E2[3]
            rate = freqin2*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 2
            rct = 40
            pos = l
            gcn = 'high'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # NEW 41) Reaction for OH + O -> H2O + O on high
        if osite == 3 and l41 == 'true':
            m = osite
            delta_total = E16[2]
            rate = freq16*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 5
            rct = 41
            pos = l
            gcn = 'high'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # NEW 42) Reaction for OH + O -> 2OH on high
        if osite == 3 and l42 == 'true':
            m = osite
            delta_total = E2[2]
            rate = freq2*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 6
            rct = 42
            pos = l
            gcn = 'high'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # NEW 43) Reaction for 2OH -> OH + O on high
        if osite == 6 and l43 == 'true':
            m = osite
            delta_total = E2[2]-E2[3]
            rate = freqin2*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 3
            rct = 43
            pos = l
            gcn = 'high'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # NEW 44) Reaction for H2O + O -> H2O + OH on high
        if osite == 5 and l44 == 'true':
            m = osite
            delta_total = E2[2]
            rate = freq2*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 7
            rct = 44
            pos = l
            gcn = 'high'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # NEW 45) Reaction for H2O + OH -> H2O + O on high
        if osite == 7 and l45 == 'true':
            m = osite
            delta_total = E2[2]-E2[3]
            rate = freqin2*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 5
            rct = 45
            pos = l
            gcn = 'high'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # NEW 46) Reaction for 2OH -> H2O + OH on high
        if osite == 6 and l46 == 'true':
            m = osite
            delta_total = E16[2]
            rate = freq16*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 7
            rct = 46
            pos = l
            gcn = 'high'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
        # NEW 47) Reaction for H2O + OH -> 2H2O on high
        if osite == 7 and l47 == 'true':
            m = osite
            delta_total = E16[2]
            rate = freq16*exp(-delta_total*inv_kb/temperature) # standard rate of reaction
            n = 0
            rct = 47
            pos = l
            gcn = 'high'
            moves_list.append(Move(rate,m,n,rct,pos,gcn))
           

    '''Pick one of the Move obj using a "Monte Carlo" procedure'''
    # Total rate: sum each move's rate
    R=0                             
    for move in moves_list:
        R += move.rate
    # Pick a random number between (0,1)    
    chi = 0.0
    while chi == 0.0:
        chi = random.random()
    t = log(1/chi)/R     # time interval
    t_total+=t           # total time taken
    # Realize the corresponding move
    chi = random.random()
    norm_rsum=0                        
    for move in moves_list:
        norm_rsum += move.rate/R
        if norm_rsum >= chi:
            if move.react == 0: # in the case of deposition we reduce the number of particles to be dropped by 1
                lsites[move.position] = move.fi_state
                drop +=1
                sitcover +=1

            if move.react == 7:
                msites[move.position] = move.fi_state
                drop +=1
                sitcover+=1

            if move.gcn == 'med' and move.react == 3: # here we move O2 from med to low GCN
                lsites[random.choice(freels)] = 1
                
            if move.gcn == 'low' and move.react == 4: # here we move O2 from low to med GCN
                msites[random.choice(freems)] = 1

            if move.gcn == 'med' and move.react == 6: # OOH (med) -> OH + O (low)
                lsites[random.choice(freels)] = 3

            if move.gcn == 'low' and move.react == 11: # OH + O (low) -> OOH (med)
                msites[random.choice(freems)] = 4

            if move.gcn == 'high' and move.react == 20: # movement of reactant O2 on the surface
                msites[random.choice(freems)] = 1

            if move.gcn == 'high' and move.react == 21: # movement of reactant OOH on the surface
                msites[random.choice(freems)] = 4

            if move.gcn == 'med' and move.react == 24: # movement of reactant O2 on the surface
                hsites[random.choice(freehs)] = 1

            if move.gcn == 'med' and move.react == 25: # movement of reactant OOH on the surface
                hsites[random.choice(freehs)] = 4

            if move.gcn == 'low': # evolution of site if LGCN
                lsites[move.position] = move.fi_state

            if move.gcn == 'med': # evolution of site if MGCN
                msites[move.position] = move.fi_state

            if move.react == 17: # additional action for the event of water production
                lsites[move.position] = move.fi_state
                water +=1
                sitcover -=1

            if move.gcn == 'high':
                hsites[move.position] = move.fi_state

            if move.react == 18:
                hsites[move.position] = move.fi_state
                drop +=1
                sitcover+=1
                                        
            if (move.react == 0):
                typcobindl +=1
            if (move.react == 7):
                typcobindm +=1
            if (move.react == 18):
                typcobindh +=1

            if (move.react == 1):
                typco1 +=1
            if (move.react == 2) or (move.react == 12) or (move.react == 13):
                typco2+=1

            if (move.react == 3):
                typcodiffO2ml += 1
            if (move.react == 4):
                typcodiffO2lm += 1

            if (move.react == 5) or (move.react == 19):
                typco5 += 1
            if (move.react == 6):
                typco6 += 1
            if (move.react == 8):
                typcoin1 += 1 
            if (move.react == 9) or (move.react == 15) or (move.react == 23):
                typcoin2 += 1 
            if (move.react == 10) or (move.react == 22):
                typcoin5 += 1 
            if (move.react == 11):
                typcoin6 += 1 
            if (move.react == 14) or (move.react == 16) or (move.react == 17):
                typco16 += 1 
            if (move.react == 20):
                typcodiffO2hm += 1
            if (move.react == 21):
                typcodiffOOHhm += 1
            if (move.react == 24):
                typcodiffO2mh += 1
            if (move.react == 25):
                typcodiffOOHmh += 1

## Time step register
            if (move.react == 0):
                dtypcobindl += t
            if (move.react == 7):
                dtypcobindm += t
            if (move.react == 18):
                dtypcobindh += t

            if (move.react == 1):
                dtypco1 += t
            if (move.react == 2) or (move.react == 12) or (move.react == 13):
                dtypco2+= t

            if (move.react == 3):
                dtypcodiffO2ml += t
            if (move.react == 4):
                dtypcodiffO2lm += t

            if (move.react == 5) or (move.react == 19):
                dtypco5 += t
            if (move.react == 6):
                dtypco6 += t
            if (move.react == 8):
                dtypcoin1 += t 
            if (move.react == 9) or (move.react == 15) or (move.react == 23):
                dtypcoin2 += t 
            if (move.react == 10) or (move.react == 22):
                dtypcoin5 += t 
            if (move.react == 11):
                dtypcoin6 += t 
            if (move.react == 14) or (move.react == 16) or (move.react == 17):
                dtypco16 += t 
            if (move.react == 20):
                dtypcodiffO2hm += t
            if (move.react == 21):
                dtypcodiffOOHhm += t
            if (move.react == 24):
                dtypcodiffO2mh += t
            if (move.react == 25):
                dtypcodiffOOHmh += t
            reactco.write("\n"+str(itt)+'\t'+str(typco1)+'\t'+str(typco2)+'\t'+str(typco5)+'\t'+str(typco6)+'\t'+str(typco16)+'\t'+str(typcoin1)+'\t'+str(typcoin2)+'\t'+str(typcoin5)+'\t'+str(typcoin6))
            reactdiff.write("\n"+str(itt)+'\t'+str(typcobindl)+'\t'+str(typcobindm)+'\t'+str(typcobindh)+'\t'+str(typcodiffO2lm)+'\t'+str(typcodiffO2ml)+'\t'+str(typcodiffO2mh)+'\t'+str(typcodiffO2hm)+'\t'+str(typcodiffOOHmh)+'\t'+str(typcodiffOOHhm))


            dtreactco.write("\n"+str(itt)+'\t'+str(dtypco1)+'\t'+str(dtypco2)+'\t'+str(dtypco5)+'\t'+str(dtypco6)+'\t'+str(dtypco16)+'\t'+str(dtypcoin1)+'\t'+str(dtypcoin2)+'\t'+str(dtypcoin5)+'\t'+str(dtypcoin6))
            dtreactdiff.write("\n"+str(itt)+'\t'+str(dtypcobindl)+'\t'+str(dtypcobindm)+'\t'+str(dtypcobindh)+'\t'+str(dtypcodiffO2lm)+'\t'+str(dtypcodiffO2ml)+'\t'+str(dtypcodiffO2mh)+'\t'+str(dtypcodiffO2hm)+'\t'+str(dtypcodiffOOHmh)+'\t'+str(dtypcodiffOOHhm))


            al = move.fi_state
            ls[al] +=1
            arct = move.react
            react[arct] +=1
            tot_water = 2*water
            nfrees = nfreels + nfreems

            if itt%100 == 0 and itt > 0:
                lgco.write("#itt\tt\t\tO2g\tEmt\tO2a\t2O\tOH+O\tOOH\tH2O+O\t2OH\tH2O+OH\t2H2O \n")
                mgco.write("#itt\tt\t\tO2g\tEmt\tO2a\t2O\tOH+O\tOOH\tH2O+O\t2OH\tH2O+OH\t2H2O \n")
                hgco.write("#itt\tt\t\tO2g\tEmt\tO2a\t2O\tOH+O\tOOH\tH2O+O\t2OH\tH2O+OH\t2H2O \n")
                agco.write("#itt\tt\t\tO2g\tEmt\tO2a\t2O\tOH+O\tOOH\tH2O+O\t2OH\tH2O+OH\t2H2O \n")

            lgco.write(str(itt)+'\t'+str(format(t_total,'.5e'))+'\t'+str(maxdrop-drop)+'\t')
            mgco.write(str(itt)+'\t'+str(format(t_total,'.5e'))+'\t'+str(maxdrop-drop)+'\t')
            hgco.write(str(itt)+'\t'+str(format(t_total,'.5e'))+'\t'+str(maxdrop-drop)+'\t')
            agco.write(str(itt)+'\t'+str(format(t_total,'.5e'))+'\t'+str(maxdrop-drop)+'\t')
            massb = drop - water
            reactants_toplot=[]
            for st in range(8):
                statcol=0
                statcom=0
                statcoh=0
                statcoa=0
                for pl in range(nl):
                    if lsites[pl] == st:
                        statcol+=1
                        statcoa+=1
                        if st != 0:
                            massb-=1
                lgco.write(str(statcol)+'\t')
                for pm in range(nm):
                    if msites[pm] == st:
                        statcom+=1
                        statcoa+=1
                        if st != 0:
                            massb-=1
                mgco.write(str(statcom)+'\t')
                for ph in range(nh):
                    if hsites[ph] == st:
                        statcoh+=1
                        statcoa+=1
                        if st != 0:
                            massb-=1
                hgco.write(str(statcoh)+'\t')
                agco.write(str(statcoa)+'\t')
                reactants_toplot.append(statcoa)
                
            lgco.write(str(water)+'\t')
            mgco.write(str(water)+'\t')
            hgco.write(str(water)+'\t')
            agco.write(str(water)+'\t')
            
            Emt.append(reactants_toplot[0])
            O2.append(reactants_toplot[1])
            TwoO.append(reactants_toplot[2])
            OHplusO.append(reactants_toplot[3])
            OOH.append(reactants_toplot[4])
            H2OplusO.append(reactants_toplot[5])
            TwoOH.append(reactants_toplot[6])
            H2OplusOH.append(reactants_toplot[7])
            water_plot.append(water)
            

            
            time_plot.append(t_total)
            
			
            if massb != 0:
                print massb

            lgco.write('\n')
            mgco.write('\n')
            hgco.write('\n')
            agco.write('\n')

            if water == x*maxdrop/10:
                print 10*x,"%"
                x += 1
            if itt == pt*1:
                lgc = str(lsites)
                mgc = str(msites)
                hgc = str(hsites)
                #g.write(str(itt))
                #g.write(lgc)
                #g.write(mgc)
                #g.write(hgc)
                #g.write("\n")
                sareact = str(areact)
                sals = str(als)
                pt+=1
            elif maxdrop-drop < 15:
                lgc = str(lsites)
                mgc = str(msites)
                hgc = str(hsites)
                g.write("\n"+lgc)
                g.write(mgc)
                g.write(hgc)
                sareact = str(areact)
                sals = str(als)
            break
    
    # events causing the loop to end (time reaching maximum or all particles being used up)
    if water == maxdrop:
        lgco.write('\n')
        mgco.write('\n')
        hgco.write('\n')            
        break


plt.plot(time_plot, Emt, color="r")
plt.savefig("Emt_cluster"+str(ccluster))
plt.close()
plt.plot(time_plot, O2, color="g")
plt.savefig("O2_cluster"+str(ccluster))
plt.close()
plt.plot(time_plot, TwoO, color="y")
plt.savefig("2O_cluster"+str(ccluster))
plt.close()
plt.plot(time_plot, OHplusO, color="k")
plt.savefig("OH+O_cluster"+str(ccluster))
plt.close()
plt.plot(time_plot, OOH, color="b")
plt.savefig("OOH_cluster"+str(ccluster))
plt.close()
plt.plot(time_plot, H2OplusO, color="m")
plt.savefig("H20+O_cluster"+str(ccluster))
plt.close()
plt.plot(time_plot, TwoOH, color="c")
plt.savefig("2OH_cluster"+str(ccluster))
plt.close()
plt.plot(time_plot, H2OplusOH, color="burlywood")
plt.savefig("H20+OH_cluster"+str(ccluster))
plt.close()
plt.plot(time_plot, water_plot, color="b")
plt.savefig("Water_cluster"+str(ccluster))
plt.close()





at_total += t_total
aritt += itt
awater += water*2
tonar += awater/at_total
tonaracs += awater/(at_total*nl)
tonarals += awater/(at_total*(nl+nm+nh))

#g.write(" END ")

if i != 1:
    at_total2 = format(at_total,'.4f')
    tonar2 = format(tonar,'.2f') + " s^-1"
    tonaracs2 = format(tonaracs,'.2f') + " s^-1"
    tonarals2 = format(tonarals,'.2f') + " s^-1"
elif i == 1:
    at_total2 = format(at_total,'.4f')
    tonar2 = format(tonar,'.6f') + " s^-1"
    tonaracs2 = format(tonaracs,'.6f') + " s^-1"
    tonarals2 = format(tonarals,'.6f') + " s^-1"

# printing of values at the end of the simulation

print("\n")
print(desg)
print("\n")
print("Times:", at_total2)
print("\n")
print("Number of MC steps:", aritt)
print("Water produced:", awater)

summ.write("\nSummary of simulation for " + desg)
#summ.write("\n\nNanoparticle with:\n"+str(nl1)+" Low sites\n"+str(nm1)+" Medium sites\n"+str(nh1)+" High sites\n")
summ.write("Simulated time interval: ") 
summ.write(str(at_total2)+"s")
summ.write("\n\n\n")
summ.write("Water produced: ")
summ.write(str(awater))

summ.write("\n\nTurnover number**:\n")
summ.write(str(tonar2))
summ.write("\n\nTurnover number per active site**:\n")
summ.write(str(tonaracs2))
summ.write("\n\nTurnover number per all sites**:\n")
summ.write(str(tonarals2))

summ.write("\n\nTemperature:\n")
summ.write(str(temperature))
summ.write("K")

summ.write("\n")

time2 = time.time()
summ.write("Time taken to run the simulation "+ str(time2-time1))

#col.write(str(testno)+'\t'+str(nl1)+'\t'+str(nm1)+'\t'+str(nh1)+'\t'+str(nl1+nm1+nh1)+'\t'+str(coverage)+'\t'+str(awater)+'\t'+str(at_total2)+'\t'+str(time2-time1)+'\t'+str(temperature)+'\n')
#allcol.write(desg+'\t'+str(testno)+'\t'+str(nl1)+'\t'+str(nm1)+'\t'+str(nh1)+'\t'+str(nl1+nm1+nh1)+'\t'+str(coverage)+'\t'+str(awater)+'\t'+str(at_total2)+'\t'+str(time2-time1)+'\t'+str(temperature)+'\n')

#CHECK: time the simulation runs for
print("Time for simulation to run:", time2 - time1)
