from pdf import *
from cn import *
from ptm import *
from cna import *
filename = "Au2223-FCC-Dh-Ih_test.xyz"

if __name__ == '__main__':
	tic_total=time.time()
	tic = time.time()
	distances = read_trajectory(filename, r_cut) #list of all neighbor distances
	toc = time.time()
	print("Time to calculate distances: %.3f [s]" %(toc-tic))

	tic = time.time()
	spline, values = create_function(distances, r_cut)
	r_cut_cn = get_cutoff_distance(distances, r_cut)
	toc = time.time()
	print("Time to create spline and cut-off distance: %.3f [s]" %(toc-tic))

	tic = time.time()
	cn, gcn = get_coord_number(filename, r_cut_cn)
	toc = time.time()
	print("Time to create coordination number (using student method): %.3f [s]" %(toc-tic))

	tic = time.time()
	cn2 = get_coord_number_asap(filename, r_cut_cn)
	toc = time.time()
	print("Time to create coordination number (using asap3): %.3f [s]" %(toc-tic))

	tic = time.time()
	ptm_data=get_PTM_array(filename,20)
	condition=np.less(ptm_data['rmsd'],0.12)
	ptm_rmsd_for_hist=np.compress(condition,ptm_data['rmsd'])
	toc = time.time()
	print("Time to run PTM analysis (using asap3): %.3f [s]" %(toc-tic))

	traj = read(filename, index = ':')
	tic = time.time()
	atomic_cnas, cna_classes = extract_cnas(traj, r_cut_cn)
	toc = time.time()
	print("Time to run CNA analysis (using asap3): %.3f [s]" %(toc-tic))

	print("CUT-OFF DISTANCE FOR THIS STRUCTURE", r_cut_cn)
	if(cn==list(cn2)):print('BOTH CN METHODS ARE EQUAL')
	toc_total=time.time()
	print("Total time taken for the calculations: %.3f [s]" %(toc_total-tic_total))


	""" (Armand)
	All of the code further below are graphs.
	I'm still implementing the Steinhardt Parameters, at the moment they wont compile. working on it.
	I'm currently attempting to implement a scatter matrix plot of all descriptors found
	also, Claudio if you read this, I hope one day I can code as well as your cna code
	god it just works so well... Thank you for that code
	"""

	"""
	for i in range(0,6):
		condition2=np.equal(ptm_data['structure'],i)
		ptm_cna_contributions=np.compress(condition2,atomic_cnas,axis=0)
		number_of_PTM=np.shape(ptm_cna_contributions)
		ptm_CNA=plt.figure(i)
		ptm_CNA=plt.title('CNA signatures in PTM=%1d with N_PTM=%1d'%(i,number_of_PTM[0]))
		ptm_CNA=plt.xlabel('CNA classes')
		ptm_CNA=plt.xticks(np.arange(19),labels=cna_classes.keys(),rotation=45)
		ptm_CNA=plt.ylabel('N')
		ptm_CNA=plt.bar(np.arange(19), np.sum(ptm_cna_contributions,axis=0), color='g')


	histogram_pdf=plt.figure(6)
	histogram_pdf=plt.title('Pair Distribution Function')
	histogram_pdf=plt.xlabel('Distance (A)')
	histogram_pdf=plt.ylabel('N')
	histogram_pdf=plt.plot(values,spline(values))

	histogram_cn2=plt.figure(7)
	histogram_cn2=plt.title('CN Histogram')
	histogram_cn2=plt.xlabel('CN values')
	histogram_cn2=plt.ylabel('N')
	histogram_cn2=plt.hist(cn2)

	Ptm_rms_histogram=plt.figure(8)
	Ptm_rms_histogram=plt.title('PTM RMSD accepted histogram')
	Ptm_rms_histogram=plt.xlabel('RMSD accepted')
	Ptm_rms_histogram=plt.ylabel('N')
	Ptm_rms_histogram=plt.hist(ptm_rmsd_for_hist)

	Ptm_histogram=plt.figure(9)
	Ptm_histogram=plt.title('PTM accepted histogram')
	Ptm_histogram=plt.xlabel('PTM values')
	Ptm_histogram=plt.ylabel('N')
	Ptm_histogram=plt.hist(ptm_data['structure'])

	CNA_histogram=plt.figure(10)
	CNA_histogram=plt.title('CNA histogram')
	CNA_histogram=plt.xlabel('CNA classes')
	CNA_histogram=plt.xticks(np.arange(19),labels=cna_classes.keys(),rotation=45)
	CNA_histogram=plt.ylabel('N')
	CNA_histogram=plt.bar(np.arange(19), cna_classes.values(), color='g')


	plt.show()
	"""
#print('CN PER ATOM IS', cn)
#print(max(cn))
#print('CN PER ATOM IS', cn2)
#print(max(cn2))
#print(list(cn2))
#print('GCN PER ATOM IS', gcn)
