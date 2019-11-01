#!/usr/bin/env python
# coding: utf-8

# In[ ]:

##### Code written by Russell Brooks @Russellzoom #########
# A Python code to calculate Spider algorithm
#### !!!!!WARNING: this code is still not fully functional, anyone is welcome to try to fix and improve it ####

'''retrieve XYZ'''
def retrieve_position(tc_index):

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


# In[2]:


'''#Monometallic cluster CN calculations'''

import math
import matplotlib.pyplot as plt

#eta sets the size of the molecule which can land on site
eta=1.0  

r=1.4450000999999999
d_min=r+r*eta*1.2                           # <-----m ay need to *1.2 for discrete CN
pos_dat = retrieve_position(range(0,100))   # <-- list of snapshots with the range chosen by the number of snapshots in xyz file

snapshots=pos_dat[0]
snap_no = pos_dat[1]                                #Number of snapshot
atoms = pos_dat[2]                                 #list of chemical species e.g. ['Au', 'Ag']
atom_no = pos_dat[3] 
atom_tot = sum(atom_no)
snap=0                        #only consider the first snap/frame here
m_pwr = 12
n_pwr = 6
pdf_s = []
cn_dict={}
pdf_dict = {}                    #dictionary of pdfs
cn = [0]*atom_tot
nn_list = {i:[] for i in range(0,atom_tot)}
#nn_list=[]
for i in range(atom_tot-1):
    pdf_i=[]
    atom_i = snapshots[snap][i][0] #element of atom i
    for j in range(i+1,atom_tot):
        atom_j = snapshots[snap][j][0]           #element of atom j
        dx = float(snapshots[snap][j][1]) - float(snapshots[snap][i][1])
        dy = float(snapshots[snap][j][2]) - float(snapshots[snap][i][2])
        dz= float(snapshots[snap][j][3]) - float(snapshots[snap][i][3])
        drij = math.sqrt(dx**2 + dy**2 + dz**2)
        pdf_i.append(drij)
        if drij <=d_min:        
            cn[i] += 1
            cn[j] += 1
            nn_list[i].append(j)
            nn_list[j].append(i)
            
        else:
            r0 = 0.147*d_min*math.sqrt(2)
            frac = (drij - d_min)/r0
            cn[i] += (1 - frac**n_pwr)/(1 - frac**m_pwr)
            if (1 - frac**n_pwr)/(1 - frac**m_pwr)>0.9:
                nn_list[i].append(j)
            cn[j] += (1 - frac**n_pwr)/(1 - frac**m_pwr)
            if (1 - frac**n_pwr)/(1 - frac**m_pwr)>0.9:
                nn_list[j].append(i)

    #print(len(pdf_i))
    cn_dict[i] = cn[i]
    pdf_dict[i]= pdf_i
if nn_list[i]==[]: #these are atoms with zero neighbours 
        print("no neighbours", i)
else:print("no anomolies")
print(snapshots[snap][20])
print(len(nn_list[20]))


# In[3]:


'''Radial distribtion calculation '''

                           
rdf_s = []  #List of radial distances in order of atom number

#Need to find the centre of mass for the cluster

comx_list = []
comy_list = []
comz_list = []
com_x = 0                                               #COMx = (m1x1 + m2x2+ ...+ mnxn)/(m1 + m2 +....+ mn)
com_y = 0
com_z = 0
sum_mass = 0
m_i = 108.0
for i in range(atom_tot):
    com_x += float(snapshots[snap][i][1])*m_i
    com_y += float(snapshots[snap][i][2])*m_i
    com_z += float(snapshots[snap][i][3])*m_i
    sum_mass += m_i
    comx_list.append(com_x/sum_mass)
    comy_list.append(com_y/sum_mass)
    comz_list.append(com_z/sum_mass)
#Radial Distribution calculation

for i in range(atom_tot):
    dx = comx_list[i] - float(snapshots[snap][i][1])
    dy = comy_list[i] - float(snapshots[snap][i][2])        
    dz = comz_list[i] - float(snapshots[snap][i][3])
    dr_com = math.sqrt(dx**2 + dy**2 + dz**2)
    rdf_s.append(dr_com)
print("calculated RD's")


# In[4]:


'''atop GCN calculation'''

cn_bulk = 12
gcn = [0]*atom_tot
non_bulk=[]
for i in range(atom_tot):
    nni = nn_list[i] #list of NN atoms for i
    if cn[i] < 13:               #all atoms
        for n in nni:  #discreet approach
            gcn[i] += cn[n]                #GCNi = sum(CNj)/12
        gcn[i] = gcn[i]/cn_bulk
        gcn[i] = float("{0:.0f}".format(gcn[i]))
        if gcn[i]<=11:
            non_bulk.append(gcn[i])
print("Total atoms", atom_tot)
print("number of non bulk",len(non_bulk))
print("number of bulk atoms", atom_tot-len(non_bulk))


# In[5]:


'''plot GCN & Cn vs RD'''

plt.figure(figsize=(20, 5))
plt.subplot(131)
plt.title("GCN and CN vs radial distribution")
plt.scatter(cn, rdf_s)
plt.scatter(gcn, rdf_s)
plt.xlabel("GCN/CN")
plt.ylabel("Radial distribution")
plt.legend(["CN", "GCN"])

plt.subplot(132)
plt.title("Number of atoms at each G/CN")
num_bins =20
n, bins, patches = plt.hist(cn, num_bins, facecolor='blue')
n, bins, patches = plt.hist(gcn, num_bins, facecolor='orange')

plt.xlabel("GCN/CN")
plt.ylabel("bin")
plt.show()
plt.show()


# In[ ]:


'''Gather coordinates'''
 for i in range(0,atom_tot):
    cluster.append(snapshots[snap][i])
    #Store coords in lists of every atom for graphing
    x_i=float(cluster[i][1])
    y_i=float(cluster[i][2])
    z_i=float(cluster[i][3])
    X_array.append(x_i)
    Y_array.append(y_i)
    Z_array.append(z_i)


# In[15]:


'''SPIDER ALGORITHM'''

import random 
import math

'''
ideally such parameters shoudl be read by the .pot file
epsilon start at e.g. helium and leave as parameter

	Atomic radius (Angstrom):                  Lattice constant (Angstrom):                 Mass (amu):
	Au = 1.4600000000000000                    Au = 4.0782                                  Au = 197.0
	Cu = 1.2800000000000000                    Cu = 3.6149                                  Cu = 64.0
	Ag = 1.4450000999999999                    Ag = 4.0853                                  Ag = 108.0
	Pt = 1.3850000000000000                    Pt = 3.9242                                  Pt = 195.0
	Pd = 1.3753226894078000                    Pd = 3.8899                                  Pd = 106.0
        He = 0.30
	O  = ...
	OH = ...
	O2 = ...
	NO = ...
'''
eta_val=[]  #eta sets the size of the molecule which can land on site
ss_list=[]  #Surface Sites at each eta 0=bulk 1=surface site
ss_no=[]    #sum of surface sites at each et
r=1.4450000999999999 #this works only for one type of atoms
count=0

#loop through eta values
for eta in [float(j)/10 for j in range(0,30)]: # <----choose range of eta
    count+=1
    eta_val.append(eta)
    d_min=r+r*eta*1.2    #may need to include *1.2 discrete
    surface_site=[]
    non_bulk=[]
    ss_list_eta=[]
    cluster=[]

    #Picks each atom to evaluate
    for i in range(0,atom_tot):
        cluster.append(snapshots[snap][i])
        x_i=float(cluster[i][1])
        y_i=float(cluster[i][2])
        z_i=float(cluster[i][3])
        #Select only atoms with cn<9
        #if cn[i]<9:                           
        if gcn[i]<11:
            atom_i=[x_i,y_i,z_i]
            rand_points=[]
            non_bulk.append(i)
            
            #generate random points around each surface atom at a distance d_min
            while True:
                #generates spherical angles theta and phi
                theta=random.uniform(0,2*math.pi)
                phi=random.uniform(0,math.pi)
                rand_point=[]
                #cartesian coord transform of random point and translation to atom_i
                x=(r*eta)*math.cos(theta)*math.sin(phi) + atom_i[0]
                y=(r*eta)*math.sin(theta)*math.sin(phi) + atom_i[1]
                z=(r*eta)*math.cos(phi) + atom_i[2]
                rand_point=[x,y,z]
                rand_points.append(rand_point)
                gap=[] #store the gaps that pass test between nn and rand point
                
                #measures d for all Nearest Neighbour to atom i (nn atoms to i are nn_list[i] with coordinates cluster[nn_list[i]])
                for k in nn_list[i]:            
                    dx=x-float(snapshots[snap][k][1])
                    dy=y-float(snapshots[snap][k][2])
                    dz=z-float(snapshots[snap][k][3])
                    s=math.sqrt(dx**2 + dy**2 + dz**2)
                    if s>=d_min:
                        gap.append(k)
                    else: 
                        break
                if len(gap)==len(nn_list[i]):
                    surface_site.append(i)      # contains the coordinate of surface site for graphing
                    ss_list_eta.append(1)       # append 1 for surface site

                    break 
                if len(rand_points)==5000:
                    ss_list_eta.append(0) 
                    break
        else: 
            ss_list_eta.append(0) # append 0 for bulk
            continue 
    ss_list.append(ss_list_eta)
    ss_no.append(len(surface_site))
    
    print('Surface sites=',len(surface_site),"   ", 'Non bulk=',len(non_bulk))
#print("Surface sites=", ss_list)
print("surface sites", ss_no)


# In[16]:


'''Plot for Surface atoms vs Eta'''

plt.figure(figsize=(20, 5))
print("How surface atoms compares with eta")
plt.subplot(131)
plt.scatter(eta_val, ss_no)
plt.title("scatter")
plt.xlabel("eta")
plt.ylabel("Surface atoms")
plt.subplot(132)
plt.plot(eta_val, ss_no)
plt.title('line')
plt.xlabel("eta")
plt.ylabel("Surface atoms")
plt.show()


# In[11]:


'''3D graph of coloured surface atoms with varying eta'''

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

for snap in range(0,10): # <---can produce graph for any snap frame but need to re-run surface atom calc for that snap number
    pos_dat = retrieve_position(range(0,100))   #list of snapshots
    snapshots=pos_dat[0]
    atom_no = pos_dat[3] 
    atom_tot = sum(atom_no)
    cluster=[] #[atom number][Ag, x, y, z]
    X_array=[]
    Y_array=[]
    Z_array=[]
    
    for i in range(0,atom_tot):
        cluster.append(snapshots[snap][i])
        #Store coords in lists of every atom for graphing
        x_i=float(cluster[i][1])
        y_i=float(cluster[i][2])
        z_i=float(cluster[i][3])
        X_array.append(x_i)
        Y_array.append(y_i)
        Z_array.append(z_i)
        
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    colours=[]

    for i in ss_list[snap]:
        if i==1:
            colours.append("blue")
        else: colours.append("orange")
    ax.scatter(X_array, Y_array, Z_array, c=colours ,s=100, marker='o')
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')

    plt.show()


# In[ ]:



