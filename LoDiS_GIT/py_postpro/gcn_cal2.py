import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as pll
from collections import Counter
from mpl_toolkits.mplot3d import Axes3D

def gcn_atop(atom_tot, cn, nn_list, col):
    cn_bulk = 12
    gcn = [0]*atom_tot
    gcn_col = ['grey']*atom_tot
    for i in range(atom_tot):
        nni = nn_list[i]             #list of NN atoms for i
        if cn[i] < 11:               #surface atoms
            for n in nni:              #discreet approach
                gcn[i] += cn[n]                #GCNi = sum(CNj)/12
            gcn[i] = gcn[i]/cn_bulk
            gcn[i] = float("{0:.1f}".format(gcn[i]))

            #Color coding for 3D bar chart
            if gcn[i] >= 0 and gcn[i]< 2:
                gcn_col[i] = col[0]
            elif gcn[i] >= 2 and gcn[i] <4:
                gcn_col[i] = col[1]
            elif gcn[i] >= 4 and gcn[i] < 6:
                gcn_col[i] = col[2]
            elif gcn[i] >=6 and gcn[i] < 8:
                gcn_col[i] = col[3]
            elif gcn[i]>= 8 and gcn[i] <10:
                gcn_col[i] = col[4]
            elif gcn[i]>= 10:
                gcn_col[i] = col[5]
    gcn_occur = dict(Counter(gcn))
    return [gcn, gcn_col, gcn_occur]

def r_vs_cni_plot(tc_list, snap_no, cn_dict,
                  atom_tot, rd_d, min_ad,
                  bin_max, bin_width, plot_unit):
    plt.figure()
    print()
    if snap_no > 1:                                          #Reduce the number of CN snapshots to 2
        cn_snaps = [min(tc_list), max(tc_list)]
        cn_snap_no = 2
        print('Plotting r vs CNi graphs for %s%s and %s%s' % (cn_snaps[0], plot_unit, cn_snaps[1], plot_unit))
    else:                                                 #If only one snapshot was inputted
        cn_snaps = tc_list
        cn_snap_no = 1
        print('Plotting r vs CNi graph for %s' % cn_snaps[0])
    pr_dict, pc_dict = {}, {}
    max_rads, max_cnis = [], []
    for snap_id in cn_snaps:
        rad_plot, cn_plot = [], []
        for i in range(atom_tot):
            rad_i, cn_i = rd_d[snap_id][i], cn_dict[snap_id][0][i]
            if rad_i >= min_ad and rad_i <= bin_max+0.5*bin_width:
                rad_plot.append(rad_i)
                cn_plot.append(cn_i)
        pr_dict[snap_id] = rad_plot
        pc_dict[snap_id] = cn_plot
        max_cnis.append(max(cn_plot))
        max_rads.append(max(rad_plot))
    cnmax = max(max_cnis)
    rmax = max(max_rads)
    for i in range(cn_snap_no):
        snap_id = cn_snaps[i]
        radial = pr_dict[snap_id]
        cn = pc_dict[snap_id]
        ax = plt.subplot(cn_snap_no, 1,i+1,)
        ax.scatter(radial, cn, marker = 'x')
        plt.annotate(('%.2f%s' %(snap_id, plot_unit)), xy = (0.1, 0.20), xycoords = 'axes fraction', size = 10)
        if cn_snap_no > 1:
            if i+1 < cn_snap_no:
                plt.setp(ax.get_xticklabels(), visible=False)
                plt.setp(ax.get_yticklabels(), visible=False)
            else:
                ax.set_xlabel('Radial distance (Angstrom)')
                ax.set_ylabel('CNi')
        else:
            ax.set_xlabel('Radial distance (Angstrom)')
            ax.set_ylabel('CNi')
        ax.set_xlim(min_ad,rmax)
        ax.set_ylim(0, cnmax)

    plt.tight_layout()
    plt.savefig( 'r_vs_cn.png')
    plt.close()



def gcn_atop_plot(cn_dict, tc_list, snap_no,
            rd_d, plot_unit, atom_tot,
            min_ad, bin_max, bin_width, col):
    print()
    print('Plotting 3D Scatter for r vs CNi vs GCNi')
    if snap_no > 1:                                          #Reduce the number of CN snapshots to 2
        cn_snaps = [min(tc_list), max(tc_list)]
        cn_snap_no = 2
    else:                                                 #If only one snapshot was inputted
        cn_snaps = tc_list
        cn_snap_no = 1

    """
    Obtain new sets of r, CNi and GCNi data for surface atoms (GCNi >0) which fall within
    the predetermined range set by the user
    """
    plotr_dict, plotc_dict, plotg_dict = {}, {}, {}
    plot_col = {}
    min_rads, min_cnis, max_cnis, min_gcnis, max_gcnis = [], [], [], [], []         #lists of min and max r, CI and GCNi
    for snap_id in cn_snaps:
        rad_plot, cn_plot, gcn_plot = [], [], []
        col_plot = []
        for i in range(atom_tot):
            rad_i, gcn_i, cn_i = rd_d[snap_id][i], cn_dict[snap_id][1][i], cn_dict[snap_id][0][i]
            cn_coli = cn_dict[snap_id][2][i]
            if rad_i >= min_ad and rad_i <= bin_max+0.5*bin_width and gcn_i != 0: #Set range and remove gcni = 0
                rad_plot.append(rad_i)
                cn_plot.append(cn_i)
                gcn_plot.append(gcn_i)
                col_plot.append(cn_coli)
        if rad_plot == []:
            print('No surface atoms detected, selected upper bound radial distance is likely too low')
            sys.exit()
        plotr_dict[snap_id] = rad_plot
        plotc_dict[snap_id] = cn_plot
        plotg_dict[snap_id] = gcn_plot
        plot_col[snap_id] = col_plot
        min_rads.append(min(rad_plot))
        min_cnis.append(min(cn_plot))
        max_cnis.append(max(cn_plot))
        min_gcnis.append(min(gcn_plot))
        max_gcnis.append(max(gcn_plot))

    min_r = min(min_rads)                      #lower bound of x-axes where  GCNi > 0
    cnmin, cnmax = min(min_cnis), max(max_cnis)     #the min and max value for CNi
    gcnmin, gcnmax = min(min_gcnis), max(max_gcnis)                     #upper bound of GCNi
    fig = plt.figure()             #dpi = 600
    zplot = [0]*atom_tot
    for i in range(cn_snap_no):
        snap_id = cn_snaps[i]
        cn = plotc_dict[snap_id]
        gcn = plotg_dict[snap_id]
        cn_col = plot_col[snap_id]
        radial_d = plotr_dict[snap_id]

        ax = plt.subplot(cn_snap_no,1,i+1, projection = '3d')
        #ax.bar3d(radial_d, cn, zplot, dx = 0.1, dy = 0.2, dz = gcn, color = cn_col)
        ax.scatter(radial_d, cn, gcn, zdir = 'z', s =20, color = cn_col)
        ax.set_zlim3d([gcnmin, gcnmax])
        ax.set_xlim3d([min_r, bin_max + 0.5*bin_width])
        ax.set_ylim3d([cnmin,cnmax])
        plt.tight_layout()
        plt.annotate(('%.2f%s' %(snap_id, plot_unit)), xy = (0.04, 0.80), xycoords = 'axes fraction') #size = 4
        ax.set_xlabel('Radial distance (A)')   #fontsize =5
        leg_prox = []
        for i in range(len(col)):
            lp = pll.Line2D([0], [0], c=col[i], linestyle = 'none', marker = 'o')
            leg_prox.append(lp)
        leg_div =  [u'GCNi <2', u'2 \u2264 GCNi < 4 ', u'4 \u2264 GCNi < 6',
                   u'6 \u2264 GCNi < 8', u'8 \u2264 GCNi < 10', u'GCNi \u2265 10']
        #set which proxies are visible given the presence of their repective colours in the graph
        leg_prox = [leg_prox[c] for c in range(len(col)) if col[c] in cn_col]
        leg_div = [leg_div[c] for c in range(len(col)) if col[c] in cn_col]
        ax.legend(leg_prox, leg_div, prop = {'size': 8}, numpoints = 1)
        ax.xaxis._axinfo["grid"].update({"linewidth":2, "color" : "black"})
        ax.yaxis._axinfo["grid"].update({"linewidth":2, "color" : "black"})
        ax.zaxis._axinfo["grid"].update({"linewidth":2, "color" : "black"})
        ax.set_ylabel('CNi')       #fontsize = 5
        ax.set_zlabel('GCNi')
        #ax.tick_params(labelsize=4)
        ax.azim, ax.elev = [-43,27]

    #plt.savefig('rd_vs_cn_vs_gcn.png', dpi=600)
    #If fort.56 exits, do not show figures till meV vs I graph is plotted
    if os.path.isfile('fort.56'):
        return
    else:
        print('VACF file *fort.56* not detected within this directory')
        plt.show()


def gcn_genome(cn_dict, tc_list, snap_no, atoms, atom_no, plot_unit):

    fout = open('gcn_genome.txt', 'w+')
    fout.write('%s %s  %s %s \n\n' % (atoms[0], atom_no[0], atoms[1], atom_no[1]))
    if plot_unit == 'K':
        id_txt = 'T (K)'
    else:
        id_txt = 't (ps)'
    for snap_id in tc_list:
        gcn_occur = cn_dict[snap_id][3]
        fout.write('%s = %s \n' % (id_txt, snap_id))

        gcn_keys = list(gcn_occur.keys())
        gcn_keys.sort()
        for k in gcn_keys:
            if k != 0:
                fout.write('%s \t %s \n' % (k, gcn_occur[k]))
        fout.write('\n')
    fout.close()
