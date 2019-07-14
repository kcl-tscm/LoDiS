import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pt
import matplotlib.lines as mpll
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d


def pos2_visual(posfile2, atom_rad, elem1, elem2):

    err_out = 0
    #Read .xyz file and array
    f_init = open(posfile2,'r')      # Original Structure

    f_read = f_init.readlines()
    charrep = [' ','\t', '\n', ',']            #Clean file of any characters
    N2_r = f_read[0].split()            # Number of atoms
    N2 = int(N2_r[0])
    elem_i = []                               # Element of atom i
    x_i, y_i, z_i = [], [], []
    for f in range(N2+2):
        for c in charrep:
            f_read[f] = f_read[f].replace(c, ' ')
        f_read[f] = f_read[f].split()
        if f >= 2:
            elem_i.append(f_read[f][0])
            x_i.append(float(f_read[f][1]))
            y_i.append(float(f_read[f][2]))
            z_i.append(float(f_read[f][3]))

    f_read = f_read[2:N2+2]
    c_elem = []
    for e in elem_i:
        if e not in c_elem:
            c_elem.append(e)
    if len(c_elem) == 1:
        c_elem = c_elem*2
    c_elem.sort()
    eldet = []
    for e in c_elem:
        if e not in atom_rad:
            err_out = 'Elements in the .xyz are not contained in the .pot file'
            return c_elem[0], c_elem[1], N2, err_out
        else:
            eldet.append(e)

    if eldet[0] == eldet[1]:
        bimon2 = False
    else:
        bimon2 = True

    rad_i = [atom_rad[elem_i[i]] for i in range(N2)]
    area_i = [math.pi*r**2 for r in rad_i]
    col = ['c']*N2


    for i in range(N2):
        if elem_i[i] not in eldet:
            err_out = 'Anomalous atom detected: atom %s element %s'%(i, elem_i[i])
            return c_elem[0], c_elem[1], N2, err_out
        if elem_i[i] != elem1:
            col[i] = 'yellow'

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
    if bimon2:
        scatter1_proxy = mpll.Line2D([0],[0], linestyle="none", c='c', marker = 'o')
        scatter2_proxy = mpll.Line2D([0],[0], linestyle="none", c='yellow', marker = 'o')
        if eldet[0] == elem1:
            ax.legend([scatter1_proxy, scatter2_proxy], [eldet[0], eldet[1]], numpoints = 1)
        else:
            ax.legend([scatter1_proxy, scatter2_proxy], [eldet[1], eldet[0]], numpoints = 1)
    else:
        if eldet[0] == elem1:
            scatter1_proxy = mpll.Line2D([0],[0], linestyle="none", c='c', marker = 'o')
            ax.legend([scatter1_proxy], [eldet[0]], numpoints = 1)
        else:
            scatter1_proxy = mpll.Line2D([0],[0], linestyle="none", c='yellow', marker = 'o')
            ax.legend([scatter1_proxy], [eldet[0]], numpoints = 1)


    plt.show()

    f_init.close()

    return c_elem[0], c_elem[1], N2, err_out
