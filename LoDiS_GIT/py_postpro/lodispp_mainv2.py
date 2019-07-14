import os
import sys
import math
import pp_init2 as ppi
import numpy as np
import bin_postplot2 as bpp
import gcn_cal2 as gcn
#Necessary inputs to analyise the data
"""
Atomic radius (Angstrom):                  Lattice constant (Angstrom):                 Mass:
Au = 1.4600000000000000                    Au = 4.0782                                  Au = 197.0
Cu = 1.2800000000000000                    Cu = 3.6149                                  Cu = 64.0
Ag = 1.4450000999999999                    Ag = 4.0853                                  Ag = 108.0
Pt = 1.3850000000000000                    Pt = 3.9242                                  Pt = 195.0
Pd = 1.3753226894078000                    Pd =                                         Pd = 106.0
"""
radcons = {'Au':1.4600, 'Cu':1.2800, 'Ag':1.4450, 'Pt':1.3850}
#latcons = {'Au': 4.07820, 'Cu': 3.61490, 'Ag': 4.08530, 'Pt': 3.92420}
masscons = {'Au': 197.0, 'Cu':64.0, 'Ag':108, 'Pt':195.0}

initial = ppi.postprocess_init()
ev_dat = ppi.retrieve_energy(initial)

itmd_check = str(initial[0])                        #Melting or Freezing?
min_ad = float(initial[2])                         #Minimum analysed pair distance
cutoff_mult = float(initial[3])                    #the factor to multiply with the atomic radius to determine the cutoff distance (usually 4 or 5)
bin_mult = float(initial[4])                       #The factor to mulitply with the radius to determine the bin width




if itmd_check.lower() in ['y','yes']:                       #produce pdf(T)/rdf(T)
	#Read Snapshots of movie.xyz close to selected temperature values
	simtlist = [float(i[6]) for i in ev_dat]                    #list of all temperatures
	plot_unit = 'K'
	print('Using temperature dependence')

elif itmd_check.lower() in ['n', 'no', 'auto']:
	simtlist = [float(i[0]) for i in ev_dat]
	plot_unit = 'ps'
	print('Using time dependence')
else:
	print('Error: Incorrect input for itMD check')
	sys.exit()
tc_list = []                               #list of pdf/rdf temperatures or times
tc_index = []                                #list of temperature/time indicies close to those chosen by the user

if len(initial[1]) > 1:                      #If a set of t/T are inputed:
	print('List of snapshots inputted for calculation:')
	t_choice = [float(i) for i in initial[1]]       #list of t/T
	t_choice.sort()       #sort into ascending order
	for tc in t_choice:
		tsnapi= min(simtlist, key = lambda x: abs(x - tc))  #obtain the closest t/T
		tc_list.append(tsnapi)
		tc_index.append(simtlist.index(tsnapi))
		print(' %s '% tsnapi)
elif initial[1][0].lower() == 'auto':              #If set to 'auto' (no t/T selected):
	print('Automatic run selected, generating 4 snapshots:')
	tmin = min(simtlist)
	tmax = max(simtlist)
	tcrun = np.linspace(tmin,tmax,4)
	for tc in tcrun:
		tc_i = min(simtlist, key = lambda x: abs(x - tc))   #find closest t/T
		tc_list.append(tc_i)
		tc_index.append(simtlist.index(tc_i))
	print(tc_list)

else:                                              #If only one t/T given:
	print('One snapshot selected for calculation:')
	t_choice = float(initial[1][0])
	tsnap = min(simtlist, key = lambda x: abs(x - t_choice))  #obtain the closest temperature
	tc_list = [tsnap]
	tc_index = [simtlist.index(tsnap)]
	print(tsnap)

print()
#Read movie.xyz
pos_dat = ppi.retrieve_position(initial, tc_index)


snapshots = pos_dat[0]                              #list of snapshots
snap_no = pos_dat[1]                                #Number of snapshot
atoms = pos_dat[2]                                 #list of chemical species e.g. ['Au', 'Ag']
atom_no = pos_dat[3]                              #list of number of atoms of each species e.g [111, 36]
atom_tot = atom_no[0] + atom_no[1]                  #Total number of atoms
duo_rad = {atoms[0]:radcons[atoms[0]], atoms[1]:radcons[atoms[1]]}    #dictionary with both species and their radius
atom_rad = max([radcons[atoms[0]], radcons[atoms[1]]])           #The larger atomic radius between the 2 species
bin_width = atom_rad*bin_mult                          #The bin width
rc = atom_rad*cutoff_mult                              #The cutoff distance

#Check nature of cluster
if atoms[0] != atoms[1]:
	print('Bimetallic cluster detected')
	bim_check = True

else:
	print('Monometallic cluster detected')
	bim_check = False


#Need to find the centre of mass for each snapshot

comx_list = []
comy_list = []
comz_list = []
for snap in range(snap_no):
	com_x = 0                                               #COMx = (m1x1 + m2x2+ ...+ mnxn)/(m1 + m2 +....+ mn)
	com_y = 0
	com_z = 0
	sum_mass = 0
	for i in range(atom_tot):
		atom_i = str(snapshots[snap][i][0])
		if atom_i in atoms:
			m_i = masscons[atom_i]
		else:
			print('Unexpected element within the cluster, system exiting')
			print(atom_i)
			sys.exit()

		com_x += float(snapshots[snap][i][1])*m_i
		com_y += float(snapshots[snap][i][2])*m_i
		com_z += float(snapshots[snap][i][3])*m_i
		sum_mass += m_i
	comx_list.append(com_x/sum_mass)
	comy_list.append(com_y/sum_mass)
	comz_list.append(com_z/sum_mass)


#Initialise the bins
bin_no = int((rc - min_ad)/bin_width) +1                       #Number of bins
bin_minx = min_ad + 0.5*bin_width                                     #Minimum central bin value
bin_max = bin_minx + bin_width*bin_no
bin_x = np.linspace(bin_minx,bin_max, bin_no)


col = ['r', 'y', 'g', 'c', 'b', 'black']
#Generate pdfs ,rdfs and CNi's
if bim_check:
	cn_dict = bpp.bim_pddf(snapshots, snap_no, tc_list,
						   atoms, atom_tot, bin_x, bin_width,
						   plot_unit, duo_rad, col, min_ad, bin_max)
else:
	cn_dict = bpp.mon_pddf(snapshots, snap_no, tc_list,
						   atoms, atom_tot, bin_x, bin_width,
						   plot_unit, duo_rad, col, min_ad, bin_max)

rd_d = bpp.rddf(snapshots, snap_no, tc_list,
	comx_list, comy_list,
 	comz_list, atom_tot, bin_x,
	bin_width,plot_unit, min_ad, bin_max)              #Plots rdfs and returns radial distances

#Generate plot of radial distance vs CNi
gcn.r_vs_cni_plot(tc_list, snap_no, cn_dict,
				  atom_tot, rd_d, min_ad,
				  bin_max, bin_width, plot_unit)

#Generate 3D plot of radial distance vs CNi vs GCNi for surface atoms
gcn.gcn_atop_plot(cn_dict, tc_list, snap_no, rd_d,
				  plot_unit, atom_tot,
				  min_ad, bin_max, bin_width, col)

#Write to file the GCNi occurrences
gcn.gcn_genome(cn_dict, tc_list, snap_no, atoms, atom_no, plot_unit)

print('Process completed')