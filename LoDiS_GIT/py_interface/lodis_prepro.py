import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pt
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d

#Read .xyz file and array
in_file = raw_input('Give the .xyz file name to be read: ')
f_init = open(in_file,'r')      # Original Structure

atom_rad = {'Pt': 1.385, 'Ag':1.44, 'Pd':1.375, 'Au':1.46}
mgo_pres = True

f_read = f_init.readlines()
charrep = [' ','\t', '\n', ',']
for c in charrep:
    if c in f_read[0]:
        f_read[0].replace(c, ' ')
    if c in f_read[1]:
        f_read[1].replace(c, ' ')

N= int(f_read[0])          # Number of atoms
f_elem = f_read[1].split()    #[elem1, elem2]
eldet = []
for e in f_elem:
    if e not in atom_rad:
        print('Unknown atom detected')
        sys.exit()
    else:
        eldet.append(e)

f_read = f_read[2:]                   #Remove headers
elem_i = []
x_i, y_i, z_i = [], [], []
for l in f_read:
    for c in charrep:
        if c in l:
            l.replace(c, ' ')
    l = l.split()
    elem_i.append(l[0])                     #list of atom elements
    x_i.append(float(l[1]))
    y_i.append(float(l[2]))
    z_i.append(float(l[3]))
rad_i = [atom_rad[elem_i[i]] for i in range(N)]
area_i = [math.pi*r**2 for r in rad_i]
col = ['b']*N
print('%s\t%s\t%s'%(x_i[:3], y_i[:3], z_i[:3]))
print('%s\t%s'%(elem_i[:3], rad_i[:3]))
print()
for i in range(N):
    if elem_i[i] not in eldet:
        print('Anomalous atom detected: atom %s element %s'%(i, elem_i[i]))
        print('Halting process')
        sys.exit()
    if elem_i[i] == eldet[1]:
        col[i] = 'r'
#plot a 3D representation of the cluster

fig = plt.figure()
ax=plt.subplot(111,projection='3d')
#plt.ion()
"""
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
    rad_s.append((ratio * np.sqrt(2)*i) ** 2)

sca.remove()
ax.scatter(x_i, y_i, z_i, zdir = 'z', s = rad_s, color = 'b', alpha = 1, edgecolor = 'black')
ax.set_xlim(minx, maxx)   #reset range for cluster
ax.set_ylim(miny, maxy)
ax.set_zlim(minz, maxz)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

if mgo_pres:
    rect = pt.Rectangle((minx, miny), maxx-minx, maxy-miny, facecolor= 'none', edgecolor= 'r')
    ax.add_patch(rect)
    art3d.pathpatch_2d_to_3d(rect, z=0, zdir = 'z')
"""
sca = ax.scatter(x_i[0], y_i[0], z_i[0], zdir = 'z', s = 10, color =col)
ax.set_xlim(-3, 3)
ax.set_ylim(-3, 3)
ax.set_zlim(-3, 3)


bbox = fig.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
width = bbox.width*fig.dpi               #width in pixels
xr_l = ax.get_xlim()                     #[min x, max x]
xrange = xr_l[1] - xr_l[0]                #width in Angstrom
ratio = width/xrange                      #units = pixels/Angstrom
rad_s = (ratio*2*rad_i[0])**2
sca.remove()
ax.scatter(x_i[0], y_i[0], z_i[0], zdir = 'z', s = rad_s, color = 'b', alpha = 1, edgecolor = 'black')
ax.set_xlim(-3, 3)
ax.set_ylim(-3, 3)
ax.set_zlim(-3, 3)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')



circle = pt.Circle((x_i[0], y_i[0]), rad_i[0], facecolor= 'none', edgecolor= 'g')
ax.add_patch(circle)
art3d.pathpatch_2d_to_3d(circle, z=z_i[0], zdir = 'z')


"""
#Obtain the axes size (in axpos) in Points
class scatter():
    def __init__(self,x_i,y_i,z_i, col, ax,size=rad_i,**kwargs):
        self.n = len(x_i)
        self.ax = ax
        self.xl = x_i
        self.yl = y_i
        self.zl = z_i
        self.col = col
        self.ax.figure.canvas.draw()
        self.size_data=size
        self.size = size
        self.sc = {}
        self.sc[0] = ax.scatter(x_i,y_i,z_i,s=self.size,**kwargs)
        self.c = 0
        self.resize(event= None)             #initial resize
        self.cid = ax.figure.canvas.mpl_connect('resize_event', self.resize)  #resize on rescale

    def resize(self,event):
        xdiff_a = abs(self.xl[1] - self.xl[0])         # dx in Angstrom
        trans = self.ax.transData.transform
        pixp_0 = trans((self.xl[0],self.yl[0]))
        pixp_1 = trans((self.xl[1],self.yl[1]))
        xdiff_p = abs(pixp_1[0] - pixp_0[0])       #dx in pixels
        ratio = xdiff_p/xdiff_a                      #units = pixels/Angstrom
        rad_s = [0]*N
        for i in range(N):
            rad_s[i] = (ratio * 2*rad_i[i]) ** 2
        if rad_s != self.size:
            self.size = rad_s
            self.sc[self.c].remove()
            self.c += 1
            self.sc[self.c] = ax.scatter(self.xl,self.yl, self.zl, s=self.size, zdir = 'z', color= self.col)
            fig.canvas.draw()
            fig.canvas.flush_events()
            self._redraw_later()

    def _redraw_later(self):
        self.timer = self.ax.figure.canvas.new_timer(interval=10)
        self.timer.single_shot = True
        self.timer.add_callback(lambda : self.ax.figure.canvas.draw_idle())
        self.timer.start()


sc = scatter(x_i,y_i,z_i, col, ax)
"""
plt.show()

f_init.close()
