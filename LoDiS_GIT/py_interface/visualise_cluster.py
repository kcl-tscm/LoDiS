import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pt
import matplotlib.lines as mpll
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d


def cluster_visual(posfile, atom_rad, mgo_pres, bimon):
    #Read .xyz file and array
    f_init = open(posfile,'r')      # Original Structure


    f_read = f_init.readlines()
    charrep = [' ','\t', '\n', ',']            #Clean file of any characters
    N_r = f_read[0].split()            # Number of atoms
    N = int(N_r[0])
    elem_i = []                               # Element of atom i
    x_i, y_i, z_i = [], [], []
    for f in range(N+2):
        for c in charrep:
            f_read[f] = f_read[f].replace(c, ' ')
        f_read[f] = f_read[f].split()
        if f >= 2:
            elem_i.append(f_read[f][0])
            x_i.append(float(f_read[f][1]))
            y_i.append(float(f_read[f][2]))
            z_i.append(float(f_read[f][3]))


    f_read = f_read[2:N+2]            #Remove the headers on the first two lines
    f_elem = []
    for e in elem_i:
        if e not in f_elem:
            f_elem.append(e)          #list of unique elements found in the .xyz file
    if len(f_elem) == 1:
        f_elem = f_elem*2            # converts ['Ag'] ==> ['Ag', 'Ag'] for instance

    eldet = []
    for e in f_elem:
        if e not in atom_rad:
            print('visualise_cluster> Unknown atom detected')
            sys.exit()
        else:
            eldet.append(e)


    rad_i = [atom_rad[elem_i[i]] for i in range(N)]
    area_i = [math.pi*r**2 for r in rad_i]
    col = ['c']*N

    if bimon:
        for i in range(N):
            if elem_i[i] not in eldet:
                print('visualise_cluster> Anomalous atom detected: atom %s element %s'%(i, elem_i[i]))
                print('visualise_cluster> Halting process')
                sys.exit()
            if elem_i[i] == eldet[1]:
                col[i] = 'yellow'
    else:
        for i in range(N):
            if elem_i[i] not in eldet:
                print('visualise_cluster> Anomalous atom detected: atom %s element %s'%(i, elem_i[i]))
                print('visualise_cluster> Halting process')
                sys.exit()

    #plot a 3D representation of the cluster

    fig = plt.figure()
    ax=plt.subplot(111,projection='3d')

    sca = ax.scatter(x_i, y_i, z_i, zdir = 'z', color =col)

    minx = min([x-r for x,r in zip(x_i,rad_i)])
    maxx = max([x+r for x,r in zip(x_i,rad_i)])
    miny = min([y-r for y,r in zip(y_i,rad_i)])
    maxy = max([y+r for y,r in zip(y_i,rad_i)])
    minz = min([z-r for z,r in zip(z_i,rad_i)])
    maxz = max([z+r for z,r in zip(z_i,rad_i)])


    ax.set_xlim(minx, maxx)   #reset range for cluster
    ax.set_ylim(miny, maxy)
    ax.set_zlim(minz, maxz)


    #Resize scatter point to atomic size
    ax.set_aspect(1)
    fig.canvas.draw()

    bbox = fig.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width = bbox.width*fig.dpi               #width in pixels
    xr_l = ax.get_xlim()                     #[min x, max x]
    xrange = xr_l[1] - xr_l[0]                #width in Angstrom
    ratio = width/xrange                      #units = pixels/Angstrom
    rad_s = []
    #Note marker size is in pixels**2 and pixels = diameter for circle
    for i in rad_i:
        rad_s.append((ratio * np.sqrt(2)*i/2) ** 2)

    sca.remove()
    ax.scatter(x_i, y_i, z_i, zdir = 'z', s = rad_s, color = col, alpha = 1, edgecolor = 'black')
    plt.axis('off')
    ax.set_xlim(minx, maxx)   #reset range for cluster
    ax.set_ylim(miny, maxy)
    ax.set_zlim(minz, maxz)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    #Legend
    if bimon:
        scatter1_proxy = mpll.Line2D([0],[0], linestyle="none", c='c', marker = 'o')
        scatter2_proxy = mpll.Line2D([0],[0], linestyle="none", c='yellow', marker = 'o')
        ax.legend([scatter1_proxy, scatter2_proxy], [eldet[0], eldet[1]], numpoints = 1)
    else:
        scatter1_proxy = mpll.Line2D([0],[0], linestyle="none", c='c', marker = 'o')
        ax.legend([scatter1_proxy], [eldet[0]], numpoints = 1)

    if mgo_pres in [1, True]:
        if minz <= 0:
            rect = pt.Rectangle((minx, miny), maxx-minx, maxy-miny, facecolor= 'none', edgecolor= 'r')
        else:
            rect = pt.Rectangle((minx, miny), maxx-minx, maxy-miny, facecolor= 'none', edgecolor= 'g')
        ax.add_patch(rect)
        art3d.pathpatch_2d_to_3d(rect, z=0, zdir = 'z')


    plt.show()

    f_init.close()
