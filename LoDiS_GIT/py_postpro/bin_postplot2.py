import os
import math
import matplotlib.pyplot as plt
import gcn_cal2 as gcn

def bim_pddf(snapshots, snap_no, tc_list,
		 atoms, atom_tot, bin_x, bin_width,
		 plot_unit, duo_rad, col, min_ad, bin_max):                      #Bimetallic cluster CN calculations

	print('Calculating pair distance distributions for the snapshots')
	fig = plt.figure()
	pdf_dict = {}                    #dictionary of pdfs
	pd_d = {}
	m_pwr = 12
	n_pwr = 6
	cn_dict = {}
	for snap in range(snap_no):
		snap_id = tc_list[snap]       #Time or Temperature of pdf
		pdf_s = []
		cn = [0]*atom_tot
		nn_list = {i:[] for i in range(atom_tot)}
		for i in range(atom_tot-1):
			atom_i = snapshots[snap][i][0]             #element of atom i
			for j in range(i+1,atom_tot):
				atom_j = snapshots[snap][j][0]           #element of atom j
				dx = float(snapshots[snap][j][1]) - float(snapshots[snap][i][1])
				dy = float(snapshots[snap][j][2]) - float(snapshots[snap][i][2])
				dz= float(snapshots[snap][j][3]) - float(snapshots[snap][i][3])
				drij = math.sqrt(dx**2 + dy**2 + dz**2)
				pdf_s.append(drij)
				if atom_i == atom_j:                        #d0 = 2r or d0 = a/sqrt(2) for fcc
					d0 = 2*duo_rad[atom_i]                  #d0*1.2 for tolerance factor for discreet CNi
				else:
					d0 = (duo_rad[atoms[0]] + duo_rad[atoms[1]])
				if drij <= d0:
					cn[i] += 1
					cn[j] += 1
					nn_list[i].append(j)
					nn_list[j].append(i)
				else:
					r0 = 0.147*d0*math.sqrt(2)
					frac = (drij - d0)/r0
					cn[i] += (1 - frac**n_pwr)/(1 - frac**m_pwr)
					cn[j] += (1 - frac**n_pwr)/(1 - frac**m_pwr)

		#Run GCN calculations
		gcn_return = gcn.gcn_atop(atom_tot, cn, nn_list, col)

		cn_dict[snap_id] = [cn]
		for g in gcn_return:
			cn_dict[snap_id].append(g)

		pdf_y = []
		for x in bin_x:                         #Loop through the bins calculated for the initial temperature
			pop = 0
			for r in pdf_s:
				if r >= x-0.5*bin_width and r<= x+0.5*bin_width:
					pop += 1
			pdf_y.append(pop)

		pdf_dict[snap_id] = pdf_y

	min_snap = min(tc_list)
	max_freq = max(pdf_dict[min_snap])
	for i in range(snap_no):
		snap_id = tc_list[i]
		ax=plt.subplot(snap_no, 1, i+1)
		ax.bar(bin_x,pdf_dict[snap_id], bin_width, align = 'center')
		plt.annotate(('%.2f%s' %(snap_id, plot_unit)), xy = (0.01, 0.80), xycoords = 'axes fraction', size = 8)
		ax.set_ylim(0,max_freq)
		ax.set_xlim(min_ad, bin_max+0.5*bin_width)
		plt.xticks(fontsize = 10)
		plt.yticks(fontsize = 10)
		if i +1< snap_no:
			plt.setp(ax.get_xticklabels(), visible=False)
			plt.setp(ax.get_yticklabels(), visible=False)
		else:
			plt.xlabel('Pair distance distribution (Angstrom)', fontsize = 20)
			plt.ylabel('Frequency', fontsize = 20)

	plt.tight_layout()
	plt.savefig('pdf_plot.png')
	plt.close()
	return cn_dict


def mon_pddf(snapshots, snap_no, tc_list,
		 atoms, atom_tot, bin_x, bin_width,
		 plot_unit, duo_rad, col, min_ad, bin_max):                      #Bimetallic cluster CN calculations

	print()
	print('Calculating pdfs, CNis and GCNis for the snapshots')
	fig = plt.figure()
	pdf_dict = {}                    #dictionary of pdfs
	pd_d = {}
	m_pwr = 12
	n_pwr = 6
	cn_dict = {}
	d0 = 2*duo_rad[atoms[0]]         #d0*1.2 for discreet CNi
	for snap in range(snap_no):
		snap_id = tc_list[snap]       #Time or Temperature of pdf
		pdf_s = []
		cn = [0]*atom_tot
		nn_list = {i:[] for i in range(atom_tot)}
		for i in range(atom_tot-1):
			atom_i = snapshots[snap][i][0]             #element of atom i
			for j in range(i+1,atom_tot):
				atom_j = snapshots[snap][j][0]           #element of atom j
				dx = float(snapshots[snap][j][1]) - float(snapshots[snap][i][1])
				dy = float(snapshots[snap][j][2]) - float(snapshots[snap][i][2])
				dz= float(snapshots[snap][j][3]) - float(snapshots[snap][i][3])
				drij = math.sqrt(dx**2 + dy**2 + dz**2)
				pdf_s.append(drij)
				if drij <= d0:
					cn[i] += 1
					cn[j] += 1
					nn_list[i].append(j)
					nn_list[j].append(i)
				else:
					r0 = 0.147*d0*math.sqrt(2)
					frac = (drij - d0)/r0
					cn[i] += (1 - frac**n_pwr)/(1 - frac**m_pwr)
					cn[j] += (1 - frac**n_pwr)/(1 - frac**m_pwr)

		#Run GCN calculations
		gcn_return = gcn.gcn_atop(atom_tot, cn, nn_list, col)

		cn_dict[snap_id] = [cn]
		for g in gcn_return:
			cn_dict[snap_id].append(g)        #atop - [cn, gcn, gcn_col, gcn_occur]

		pdf_y = []
		for x in bin_x:                         #Loop through the bins calculated for the initial temperature
			pop = 0
			for r in pdf_s:
				if r >= x-0.5*bin_width and r<= x+0.5*bin_width:
					pop += 1
			pdf_y.append(pop)

		pdf_dict[snap_id] = pdf_y

	min_snap = min(tc_list)
	max_freq = max(pdf_dict[min_snap])
	for i in range(snap_no):
		snap_id = tc_list[i]
		ax=plt.subplot(snap_no, 1, i+1)
		ax.bar(bin_x,pdf_dict[snap_id], bin_width, align = 'center')
		plt.annotate(('%.2f%s' %(snap_id, plot_unit)), xy = (0.01, 0.80), xycoords = 'axes fraction', size = 8)
		ax.set_ylim(0,max_freq)
		ax.set_xlim(min_ad, bin_max+0.5*bin_width)
		plt.xticks(fontsize = 10)
		plt.yticks(fontsize = 10)
		if i +1< snap_no:
			plt.setp(ax.get_xticklabels(), visible=False)
		else:
			plt.xlabel('Pair distance distribution (Angstrom)', fontsize = 20)
			plt.ylabel('Frequency', fontsize = 20)

	plt.tight_layout()
	plt.savefig('pdf_plot.png')
	plt.close()
	return cn_dict

def rddf(snapshots, snap_no, tc_list, comx_list, comy_list,
		 comz_list, atom_tot, bin_x,
		 bin_width,plot_unit, min_ad, bin_max):

	print()
	print('Calculating the radial distribution for the snapshots')
	fig = plt.figure()
	rdf_dict = {}                          #dictionary of rdfs
	rd_d = {}                              #dictionary of radial distances
	for snap in range(snap_no):                         #time = (s+1)*100.0000, range(snap_no)
		snap_id = tc_list[snap]
		rdf_s = []
		for i in range(atom_tot):
			dx = comx_list[snap] - float(snapshots[snap][i][1])
			dy = comy_list[snap] - float(snapshots[snap][i][2])
			dz = comz_list[snap] - float(snapshots[snap][i][3])
			dr_com = math.sqrt(dx**2 + dy**2 + dz**2)
			rdf_s.append(dr_com)

		rdf_y = []
		for x in bin_x:                         #Loop through the bins calculated for the initial temperature
			pop = 0
			for r in rdf_s:
				if r >= x-0.5*bin_width and r<= x+0.5*bin_width:
					pop += 1
			rdf_y.append(pop)

		rdf_dict[snap_id] = rdf_y
		rd_d[snap_id] = rdf_s
	min_snap = min(tc_list)                     #lowest temperature has highest frquency
	max_freq = max(rdf_dict[min_snap])
	for i in range(snap_no):
		snap_id = tc_list[i]
		ax=plt.subplot(snap_no, 1, i+1)
		ax.bar(bin_x,rdf_dict[snap_id], bin_width, align = 'center')
		plt.annotate(('%.2f%s' %(snap_id, plot_unit)), xy = (0.01, 0.80), xycoords = 'axes fraction', size = 8)
		ax.set_ylim(0,max_freq)
		ax.set_xlim(min_ad, bin_max+0.5*bin_width)
		plt.xticks(fontsize = 10)
		plt.yticks(fontsize = 10)

		if i +1< snap_no:
			plt.setp(ax.get_xticklabels(), visible=False)
			plt.setp(ax.get_yticklabels(), visible=False)
		else:
			plt.xlabel('Radial distance distribution (Angstrom)', fontsize = 20)
			plt.ylabel('Frequency', fontsize = 20)

	plt.tight_layout()
	plt.savefig('rdf_plot.png')
	plt.close()
	return rd_d