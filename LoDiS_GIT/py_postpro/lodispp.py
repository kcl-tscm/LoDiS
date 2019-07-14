import os
import sys
import math
import pp_init2 as ppi
import numpy as np
import bin_postplot2 as bpp
import gcn_cal2 as gcn
import vacf_pp as vp

def md_postpro(run_process, tstep, snap_id, temp_id, min_ad,
			   cutoff, bin_mult, natom, ndepmax,
			   radcons, masscons):

	#Necessary inputs to analyise the data
	"""
	radcons = atom_rad = distionary of elements and their radius
	masscons = atom_mass = dictionary of elements and their mass

	Atomic radius (Angstrom):                  Lattice constant (Angstrom):                 Mass:
	Au = 1.4600000000000000                    Au = 4.0782                                  Au = 197.0
	Cu = 1.2800000000000000                    Cu = 3.6149                                  Cu = 64.0
	Ag = 1.4450000999999999                    Ag = 4.0853                                  Ag = 108.0
	Pt = 1.3850000000000000                    Pt = 3.9242                                  Pt = 195.0
	Pd = 1.3753226894078000                    Pd =                                         Pd = 106.0
	"""


	ev_dat = ppi.retrieve_energy()
	if ev_dat != 'none':

		if temp_id in [1, True, 1.0]:                       #produce pdf(T)/rdf(T)
			#Read Snapshots of movie.xyz close to selected temperature values
			simtlist = [float(i[6]) for i in ev_dat]                    #list of all temperatures
			plot_unit = 'K'
			print('Using temperature dependence')

		elif temp_id in [0, False]:
			simtlist = [float(i[0]) for i in ev_dat]                    #list of all times
			plot_unit = 'ps'
			print('Using time dependence')
		else:
			print('Error: Incorrect input for itMD check')
			sys.exit()
		tc_list = []                               #list of pdf/rdf temperatures or times
		tc_index = []                                #list of temperature/time indicies close to those chosen by the user

		if len(snap_id) > 1:                      #If a set of t/T are inputed:
			print('List of snapshots inputted for calculation:')
			snap_id.sort()       #sort into ascending order
			for tc in snap_id:
				tsnapi= min(simtlist, key = lambda x: abs(x - tc))  #obtain the closest t/T
				tc_list.append(tsnapi)
				tc_index.append(simtlist.index(tsnapi))
				print(' %s '% tsnapi)
		elif snap_id == ['auto']:              #If set to 'auto' (no t/T selected):
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
			t_choice = snap_id[0]
			tsnap = min(simtlist, key = lambda x: abs(x - t_choice))  #obtain the closest temperature
			tc_list = [tsnap]
			tc_index = [simtlist.index(tsnap)]
			print(tsnap)

		print()

		#Read movie.xyz
		pos_dat = ppi.retrieve_position(tc_index)

		if pos_dat != 'none':
			snapshots = pos_dat[0]                              #list of snapshots
			snap_no = pos_dat[1]                                #Number of snapshot
			atoms = pos_dat[2]                                 #list of chemical species e.g. ['Au', 'Ag']
			atom_no = pos_dat[3]                              #list of number of atoms of each species e.g [111, 36]
			if run_process == 'Growth':
				atom_tot = sum(atom_no)
			else:
				atom_tot = natom

			duo_rad = {atoms[0]:radcons[atoms[0]], atoms[1]:radcons[atoms[1]]}    #dictionary with both species and their radius
			atom_rad = max([radcons[atoms[0]], radcons[atoms[1]]])           #The larger atomic radius between the 2 species
			bin_width = atom_rad*bin_mult                          #The bin width

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
			bin_no = int((cutoff - min_ad)/bin_width) +1                       #Number of bins
			bin_minx = min_ad + 0.5*bin_width                                     #Minimum central bin value
			bin_max = bin_minx + bin_width*bin_no
			bin_x = np.linspace(bin_minx,bin_max, bin_no)
			print()
			print('Min distance: %s, Max distance: %s, Bin width: %s'%(min_ad, bin_max+0.5*bin_width, bin_width))

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

			print('Initiating VACF calculations')
		else:
			if run_process == 'Growth':
				atom_tot = natom + ndepmax
			else:
				atom_tot = natom
	else:
		print('Checking viability to run VACF calculation')
		print('Automatic population counting not possible')
		if run_process == 'Growth':
			atom_tot = natom + ndepmax
		else:
			atom_tot = natom

	#Run VACF function to produce meV vs VDOS
	vacf_cal = vp.acf_ff(tstep, atom_tot)

	if ev_dat == 'none' and vacf_cal == 'none':
		print()
		print('Post-process analysis is impossible, please provide the necessary ')
		print('*energy.out*, *movie.xyz* and *fort.56* files before running this program')
	else:
		print('Process completed')