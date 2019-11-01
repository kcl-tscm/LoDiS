# The values below are those that are specific to the TONi structure, only the l,m,h numbers should change

desg = 'TONi'

ncat = 2 # parameter for the size of the nanoparticle

# Droplets type
# Energies sorted by 
# reactions: 
#1) O2 -> 2O
#2) 2O -> OH + O
#5) O2 -> OOH
#6) OOH -> OH + O
#16) OH + O -> H2O

freq = 6.245*10**12
freq1 = 1.07*10**14
freq2 = 4.2*10**13
freq5 = 3*10**13
freq6 = 1.31*10**17
freq16 = 1.07*10**12

freqin1 = 1.85*10**14
freqin2 = 7.73*10**14
freqin5 = 2.38*10**15
freqin6 = 2.29*10**15
freqin16 = 1.07*10**12

## First Binding energies, then diffusion barriers, then the activation energies (for the different nanoparticles)

# binding energy TO NI@I low: -1.47 eV  high: -0.09 eV
#          [TONi@I]
en_bindl = -1.42
en_bindm = -0.75
en_bindh = -0.09

# en_bind = 0.4 # [eV] binding energy for the reactants on the domains 
# en_diffC = 0.37 O2 Pt(111) Bray 2011 Langmuir

en_diff_OOH_temp = 0.7
en_diffconstO2 = 0.37 # [eV] diffusion energy for the reactants on the domains

en_diffmlO2 = en_diffconstO2
en_difflmO2 = en_bindm-en_bindl + en_diffconstO2
en_diffhmO2 = en_diffconstO2
en_diffmhO2 = en_bindh-en_bindm + en_diffconstO2


####!!!! The energies below are split by the type of nanoparticle first, then by the GCN type, then into the various reactions, and finally activation and reaction energy
stno = "no of states: 8 (7 actual + 1 empty)"
evno = "no of different events possible: 48"  

## Energies for TONi@I
E1lowact = 0.65
E1lowreact = -0.75
E2lowact = 0.75
E2lowreact = -0.52
E5lowact = 0.4   #### not used
E5lowreact = -0.1   #### not used
E6lowact = 0.4   #### not used
E6lowreact = -0.1   #### not used
E16lowact = 0.70
E16lowreact = -0.10


E1highact = 1.11   #### not used
E1highreact = -0.22   #### not used
E2highact = 0.70   #### not used
E2highreact = -0.56   #### not used
E5highact = 0.19  # very low so may cause problems
E5highreact = -0.36
E6highact = 0.84
E6highreact = -0.29
E16highact  = 1   #### not used
E16highreact  = -1    #### not used

l0 = 'true' # L: O2(g) + * -> O2(a)
l1 = 'true' # L: O2* -> 2O*
l2 = 'true' # L: 2O -> OH + O
l3 = 'true' # ML: O2 diff
l4 = 'true' # LM: O2 diff
l5 = 'false' # M: O2 -> OOH
l6 = 'false' # ML: OOH -> OH + O 
l7 = 'true' # M: O2(g) + * -> O2(a)
l8 = 'true' # L: 2O -> O2
l9 = 'true' # L: OH + O -> 2O
l10 = 'false' # M: OOH -> O2
l11 = 'false' # LM: OH + O -> OOH
l12 = 'true' # L: OH + O -> 2OH
l13 = 'true' # L: H2O + O -> H2O + OH
l14 = 'true' # L: 2OH -> H2O + OH
l15 = 'true' # L: 2OH -> OH + O
l16 = 'true' # L: OH + O -> H2O + O
l17 = 'true' # L: H2O + OH -> 2H2O + *
l18 = 'true' # H: O2(g) + * -> O2(a)
l19 = 'false' # H: O2 -> OOH
l20 = 'true' # HM: O2 diff
l21 = 'false' # HM: OOH diff
l22 = 'false' # H: OOH -> O2
l23 = 'true' # L: H2O + OH -> H2O + O
l24 = 'false' # MH: O2 diff
l25 = 'false' # MH: OOH diff
l26 = 'false' # M: O2 -> 2O
l27 = 'false' # M: 2O -> OH + O
l28 = 'false' # M: 2O -> O2
l29 = 'false' # M: OH + O -> 2O
l30 = 'false' # M: OH +O -> H2O + O
l31 = 'false' # M: OH +O -> 2OH
l32 = 'false' # M: 2OH -> OH + O
l33 = 'false' # M: H2O + O -> H2O + OH
l34 = 'false' # M: H2O + OH -> H2O + O
l35 = 'false' # M: 2OH -> H2O + OH
l36 = 'false' # M: H2O + OH -> 2H2O
l37 = 'false' # H: O2 -> 2O
l38 = 'false' # H: 2O -> OH + O
l39 = 'false' # H: 2O -> O2
l40 = 'false' # H: OH +O -> 2O
l41 = 'false' # H: OH +O -> H2O + O
l42 = 'false' # H: OH +O -> 2OH
l43 = 'false' # H: 2OH -> OH + O
l44 = 'false' # H: H2O + O -> H2O + OH
l45 = 'false' # H: H2O + OH -> H2O + O
l46 = 'false' # H: 2OH -> H2O + OH
l47 = 'false' # H: H2O + OH -> 2H2O

