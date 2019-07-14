import os
import sys
if sys.version_info[0] < 3:
    import Tkinter as tk
    import tkFont as tkfont
    import tkMessageBox
    import Tkconstants, tkFileDialog
else:
    import tkinter as tk
    import tkinter.font as tkfont
    import tkinter.messagebox as tkMessagebox
    import tkinter.constants as Tkconstants
    import tkinter.filedialog as tkFileDialog
import run_check as rnc

"""
***** Growth Frame *****
Growth parameters including the type of growth to be run.

        self.growth = {'ndepmax': tk.IntVar(),          #number of deposited atoms
                       'lcs': tk.IntVar,                #1=mono growth, 2=mixed-shell, 3=core-shell
                       'at_tipo2': tk.IntVar(),         #number of species 2 deposited
                       'elemd': tk.StringVar(),         #Chemical species of deposited atom (species 2)
                       'prob': tk.DoubleVar(),          #Probability of deposition of species 1
                       'tsorg': tk.DoubleVar(),         #temperature of the source
                       'rad': tk.DoubleVar()}           #radius of source
"""

class Growth_params(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        self.title = tk.Label(self, text = 'Growth parameters')
        self.title.grid(row=0, column=0, columnspan = 3, sticky = tk.W, pady = (0,20))
        self.title.config(font=("Courier", 18))

        lcs_label = tk.Label(self, text = 'Growth type:', padx = 10)
        lcs1 = tk.Checkbutton(self, text = 'Monometallic', variable = controller.growth['lcs'][0],
                                      anchor = tk.W, padx = 10, command = lambda:self.lcs_switch(controller, 0))
        lcs2 = tk.Checkbutton(self, text = 'Mixed-shell', variable = controller.growth['lcs'][1],
                                      anchor = tk.W, padx = 10, command = lambda:self.lcs_switch(controller, 1))
        lcs3 = tk.Checkbutton(self, text = 'Core-shell', variable = controller.growth['lcs'][2],
                                      anchor = tk.W, padx = 10, command = lambda:self.lcs_switch(controller, 2))
        lcs_label.grid(row=1, column=0, sticky=tk.W, pady = (0,10))
        lcs1.grid(row=1, column=1, sticky=tk.W, pady=(0,20))
        lcs2.grid(row=1, column=2, sticky=tk.W, pady=(0,20))
        lcs3.grid(row=1, column=3, sticky=tk.W, pady=(0,20))

        ndep_label = tk.Label(self, text = 'ndepmax', padx = 10)
        ndep_entry = tk.Entry(self, textvariable = controller.growth['ndepmax'])
        ndep_label.grid(row=2, column=0, sticky = tk.W, pady = (0,10))
        ndep_entry.grid(row=2, column=1, sticky = tk.W)

        atipo_label = tk.Label(self, text = 'at_tipo2', padx = 10)
        atipo_entry = tk.Entry(self, textvariable = controller.growth['at_tipo2'])
        atipo_label.grid(row=2, column=2, sticky = tk.W, pady = (0,10))
        atipo_entry.grid(row=2, column=3, sticky = tk.W)

        elemd_label = tk.Label(self, text = 'elemd', padx = 10)
        elemd_entry = tk.Entry(self, textvariable = controller.growth['elemd'])
        elemd_label.grid(row=3, column=0, sticky = tk.W, pady = (0,10))
        elemd_entry.grid(row=3, column=1, sticky = tk.W)

        prob_label = tk.Label(self, text = 'prob', padx = 10)
        prob_entry = tk.Entry(self, textvariable = controller.growth['prob'])
        prob_label.grid(row=3, column=2, sticky = tk.W, pady = (0,10))
        prob_entry.grid(row=3, column=3, sticky = tk.W)

        tsorg_label = tk.Label(self, text = 'tsorg', padx = 10)
        tsorg_entry = tk.Entry(self, textvariable = controller.growth['tsorg'])
        tsorg_label.grid(row=4, column=0, sticky = tk.W, pady = (0,10))
        tsorg_entry.grid(row=4, column=1, sticky = tk.W)

        rad_label = tk.Label(self, text = 'rad', padx = 10)
        rad_entry = tk.Entry(self, textvariable = controller.growth['rad'])
        rad_label.grid(row=4, column=2, sticky = tk.W, pady = (0,20))
        rad_entry.grid(row=4, column=3, sticky = tk.W)

        nextb = tk.Button(self, text = 'Next',padx =30, command = lambda: self.fcheck_proceed(controller))
        nextb.grid(row = 5, column = 4, sticky = tk.E, pady = (30, 0))
        prevb = tk.Button(self, text = 'Back',padx=30, command = lambda: controller.show_frame('Sys_params'))
        prevb.grid(row = 5, column = 3, sticky = tk.E, pady = (30, 0))

    def lcs_switch(self, controller, j):                              #Diables other buttons when one process is selected
        if controller.growth['lcs'][j].get():
            for l in range(3):
                if l != j:
                    controller.growth['lcs'][l].set(False)

    def fcheck_proceed(self, controller):
        gpass = {}
        for g in controller.growth.keys():
            if g != 'lcs':
                try:
                    controller.growth[g].get()                     #checks whether input is correct variable type
                except:
                    tkMessageBox.showinfo('Error', 'Incorrect characters for '+controller.sys_vn[i])
                    return
            else:
                for i in range(3):
                    if True not in [l.get() for l in controller.growth['lcs']]:
                        tkMessageBox.showinfo('Error', 'Must choose a growth process to run')
                        return
        if controller.growth['ndepmax'].get() <= 0:
            tkMessageBox.showinfo('Error', 'The number of deposited atoms cannot be equal to or less than zero')
        elif controller.growth['prob'].get() == 0. and controller.growth['lcs'][1].get() == True:
            tkMessageBox.showinfo('Error', 'Mixed shell growth cannot have a zero *prob*')
        elif controller.growth['at_tipo2'].get() == 0. and controller.growth['lcs'][1].get() == True:
            tkMessageBox.showinfo('Error', 'Mixed shell growth cannot have a zero atoms deposited of species 2')
        elif controller.growth['tsorg'].get() < 0:
            tkMessageBox.showinfo('Error', 'The temperature of the source cannot be zero or less')
        elif controller.growth['rad'].get() < 5.:
            tkMessageBox.showinfo('Error', 'Source radius should be greater than 5 to prevent clipping')
        else:
            if controller.growth['elemd'].get() not in [controller.sys_out['elem1'].get(), controller.sys_out['elem2'].get()]:
                if controller.growth['elemd'].get() not in controller.atom_rad.keys():
                    tkMessageBox.showinfo('Error', 'Unknown element selected for *elemd*')
                else:                                                  # monometallic cluster and system currently set as monometallic, need to change elem2=elemd
                    e_list = controller.atom_rad.keys()
                    e_list.sort()
                    controller.sys_out['elem1'].set(e_list[0])
                    controller.sys_out['elem2'].set(e_list[1])
                    controller.bimon.set(True)

            if controller.ppvar.get():
                controller.show_frame('Postpro_params')
            else:
                rnc.FinalPage.params_record(controller.frames['FinalPage'], controller)               #Generate parameter labels
                controller.show_frame('FinalPage')
