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
***** Post Analysis Frame *****
Parameters necessary to carry out analysis after LoDiS simulations.

        self.postpro = {'T_check': tk.BooleanVar(),
                        'snapshots': tk.StringVar(),
                        'd_min': tk.DoubleVar(),
                        'cutoff': tk.DoubleVar(),
                        'bin_mult': tk.DoubleVar()}
"""
class Postpro_params(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        self.title = tk.Label(self, text = 'Post-processing parameters')
        self.title.grid(row=0, column=0, columnspan = 2, sticky = tk.W, pady = (0,20))
        self.title.config(font=("Courier", 18))

    def create_ppwidgets(self, controller):
        self.tvcheck = tk.Checkbutton(self, text = 'ID = T(K)', variable = controller.postpro['T_check'],anchor = tk.W, padx = 10)
        self.tvcheck.grid(row = 1, column = 0, pady = 10, sticky = tk.W)
        if controller.run_process != 'itMD':
            controller.postpro['T_check'].set(False)
            self.tvcheck.config(state=tk.DISABLED)


        self.snapl = tk.Label(self, text= 'Snapshots t(ps)/T(K)', padx = 10)
        self.snapv = tk.Entry(self,textvariable = controller.postpro['snapshots'])
        self.snapl.grid(row=2, column = 0, sticky= tk.W, pady = (0, 10))
        self.snapv.grid(row = 2,column= 1, columnspan = 2, sticky = tk.W)


        self.mindl = tk.Label(self, text= 'min distance', padx = 10)
        self.mindv = tk.Entry(self,textvariable = controller.postpro['d_min'])
        self.mindl.grid(row=2, column = 2, sticky= tk.W, pady = (0, 10))
        self.mindv.grid(row = 2,column= 3,  sticky = tk.W)

        self.cutml = tk.Label(self, text= 'cut-off distance', padx = 10)
        self.cutmv = tk.Entry(self,textvariable = controller.postpro['cutoff'])
        self.cutml.grid(row=3, column = 0, sticky= tk.W, pady = (0, 10))
        self.cutmv.grid(row = 3,column= 1,  sticky = tk.W)

        self.bwl = tk.Label(self, text= 'bin width multiplier', padx = 10)
        self.bwv = tk.Entry(self,textvariable = controller.postpro['bin_mult'])
        self.bwl.grid(row=3, column = 2, sticky= tk.W, pady = (0, 10))
        self.bwv.grid(row = 3,column= 3,  sticky = tk.W)

        prevb = tk.Button(self, text = 'Back',padx=30, command = lambda: controller.show_frame('Sys_params'))
        prevb.grid(row = 4, column = 3, sticky = tk.E, pady = (30, 0))
        nextb = tk.Button(self, text = 'Next',padx =30, command = lambda: self.f_proceed(controller))
        nextb.grid(row = 4, column = 4, sticky = tk.E, pady = (30, 0))

    def f_proceed(self, controller):
        for i in controller.postpro.keys():
            try:
                controller.postpro[i].get()
            except:
                tkMessageBox.showinfo('Error', '*'+str(i)+'* is unknown')
            if i not in ['snapshots', 'T_check']:
                if controller.postpro[i].get() in [0,0.0, 0.]:
                    tkMessageBox.showinfo('Error', '*'+i + '* cannot be zero')
                    return
        snaps = controller.postpro['snapshots'].get()
        snaps = snaps.replace(',', ' ').split()                                #Make list
        if snaps in [[''], ['()'], ['auto']]:        #Empty or auto list
            controller.snapout = ['auto']
            rnc.FinalPage.params_record(controller.frames['FinalPage'], controller)               #Generate parameter labels
            controller.show_frame('FinalPage')
        else:
            if isinstance(snaps, list):
                try:
                    controller.snapout = [float(s) for s in snaps]                   #Check float values are given
                except:
                    tkMessageBox.showinfo('Error', 'Non-float t/T values')
                    return
                if len(controller.snapout) > 4:
                    tkMessageBox.showinfo('Error', 'Give only 4 snapshots max')
                    return
                else:
                    rnc.FinalPage.params_record(controller.frames['FinalPage'], controller)               #Generate parameter labels
                    controller.show_frame('FinalPage')
            else:
                tkMessageBox.showinfo('Error', 'Unknown value(s)')
                print (snaps)
