import os
import sys


def retrieve_energy():
	if os.path.isfile('energy.out'):
		f_name = open('energy.out', 'r')
		f_read = f_name.readlines()  # ['line1', line2',...]
		f_read = f_read[1:]  # remove column titles
		ev_dat = [y.split() for y in f_read]  # [[line1], [line2] ...], line1 = 'time', 'etot', ....
		f_name.close()
		return ev_dat
	else:
		print('No energy.out files were detected within this directory')
		print('Forgoing PDF,RDF, CNi and GCNi calculations')
		return 'none'


def retrieve_position(tc_index):
	if os.path.isfile('movie.xyz'):  # Only loop through folders with movie.xyz files
		print('Opening movie.xyz for reading')
		f_name = open('movie.xyz', 'r')
		charrep = [' ','\t', '\n', ',']            #Clean file of any characters
		f_read = f_name.readlines()
		num_lines = len(f_read)  # Total number of lines
		a_t = f_read[0].split()
		atom_tot = int(a_t[0])
		e_list = []
		atoms = []
		atom_no = []
		for f in range(num_lines):
			for c in charrep:
				f_read[f] = f_read[f].replace(c, ' ')

		f_list = [y.split() for y in f_read]  # Each line is converted into a list within i.e a list of line lists
		for f in range(atom_tot + 2):
			if f >= 2:
				e_list.append(f_list[f][0])
		for e in e_list:
			if e not in atoms:
				atoms.append(e)
				atom_no.append(1)
			else:
				atom_no[atoms.index(e)] += 1
		atoms.sort()
		if len(atoms) == 1:
			atoms = atoms*2
			atom_no = [atom_no[0], 0]

		print('Detected cluster: %s %s, %s %s'% (atoms[0], atom_no[0], atoms[1], atom_no[1]))
		snap_no = int(num_lines / (atom_tot + 2))  # number of snapshots
		f_cut = [z[:4] for z in f_list]  # Removes last column
		snaps = [f_cut[(atom_tot + 2)*n + 2:(atom_tot + 2)*(n + 1)] for n in range(snap_no)]  # Removes lines with text and seperate different snapshots in lists
		snapshots = [snaps[tc] for tc in tc_index]  # Obtain snapshots for selected temperatures or calculated times
		snap_no = len(snapshots)
		pos_dat = [snapshots, snap_no, atoms, atom_no]  # pos_dat contains snapshots and the number of snapshots.
		f_name.close()  # [*snaps,snp][*snap0,*snap1,...][*line0,*line1,...][x,y,z], * is a list variable
		return pos_dat
	else:
		print('No movie.xyz files were detected within this directory')
		print('Forgoing PDF,RDF, CNi and GCNi calculations')
		return 'none'
