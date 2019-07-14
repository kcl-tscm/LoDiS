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
***** Stacking Fault Number Frame *****
Parameters for the SFN collective variable for metadynamics runs.

        self.sfn_on = tk.BooleanVar()
        self.sfn_on.set(False)
        self.sfn = {'nsf_pwr': tk.IntVar(),                      #n
                   'msf_pwr': tk.IntVar(),                       #m
                   'dsf': tk.DoubleVar(),                        #d0 bulk lattice ref
                   'rsf': tk.DoubleVar(),                        # r0 bulk lattice ref
                   'gsfwidth': tk.DoubleVar()}                   # Guassian width

        #SFN            #C2N / d3.4N      <= A subset of SFN
        n = 6           n = 6
        m = 12          m = 12
        d = 1.354       d = 3.4
        r = 0.05        r = 0.1
        g = 25.0        g = 50.0
"""

class SFN_params(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        self.title = tk.Label(self, text = 'SFN parameters')
        self.title.grid(row=0, column=0, columnspan = 2, sticky = tk.W, pady = (0,20))
        self.title.config(font=("Courier", 18))

        c2n_check = tk.Checkbutton(self, text = 'Use C2N values', variable = controller.c2n_on,
                                      anchor = tk.W, padx = 10, command = lambda:self.sfn_switch(controller))
        c2n_check.grid(row =1, column=0, sticky = tk.W, pady = (0,20))

        n_label = tk.Label(self, text = 'n', padx = 10)
        n_entry = tk.Entry(self, textvariable = controller.sfn['nsf_pwr'])
        n_label.grid(row=2, column=0, sticky = tk.W, pady = (0,10))
        n_entry.grid(row=2, column=1, sticky = tk.W)

        m_label = tk.Label(self, text = 'm', padx = 10)
        m_entry = tk.Entry(self, textvariable = controller.sfn['msf_pwr'])
        m_label.grid(row=2, column=2, sticky = tk.W, pady = (0,10))
        m_entry.grid(row=2, column=3, sticky = tk.W)

        d_label = tk.Label(self, text = 'd0', padx = 10)
        d_entry = tk.Entry(self, textvariable = controller.sfn['dsf'])
        d_label.grid(row=3, column=0, sticky = tk.W, pady = (0,10))
        d_entry.grid(row=3, column=1, sticky = tk.W)


        r_label = tk.Label(self, text = 'r0', padx = 10)
        r_entry = tk.Entry(self, textvariable = controller.sfn['rsf'])
        r_label.grid(row=3, column=2, sticky = tk.W, pady = (0,10))
        r_entry.grid(row=3, column=3, sticky = tk.W)

        gw_label = tk.Label(self, text = 'gwidth', padx = 10)
        gw_entry = tk.Entry(self, textvariable = controller.sfn['gsfwidth'])
        gw_label.grid(row=4, column=0, sticky = tk.W, pady = (0,10))
        gw_entry.grid(row=4, column=1, sticky = tk.W)

        nextb = tk.Button(self, text = 'Next',padx =30, command = lambda: self.fcheck_proceed(controller))
        nextb.grid(row = 5, column = 4, sticky = tk.E, pady = (30, 0))
        prevb = tk.Button(self, text = 'Back',padx=30, command = lambda: self.cv_return(controller))
        prevb.grid(row = 5, column = 3, sticky = tk.E, pady = (30, 0))

    def sfn_switch(self,controller):
        if controller.c2n_on.get():
            controller.sfn['dsf'].set(3.4)
            controller.sfn['rsf'].set(0.1)
            controller.sfn['gsfwidth'].set(50.0)
        else:
            controller.sfn['dsf'].set(1.354)
            controller.sfn['rsf'].set(0.05)
            controller.sfn['gsfwidth'].set(25.0)


    def fcheck_proceed(self, controller):
        for m in controller.cn.keys():
            try:
                controller.cn[m].get()
            except:
                tkMessageBox.showinfo('Error', '*'+m+'* contains unknown character')
            if controller.cn[m].get() <= 0 :
                tkMessageBox.showinfo('Error', '*'+m+'* must be larger than zero')
                return


        if controller.metacore['num_cv'] == 2:
            if controller.metacore['collvar_name(2)'].get() != 'SFN':
                next_page = controller.cvpage_link[controller.metacore['collvar_name(2)'].get()]
                controller.show_frame(next_page)
            else:
                if controller.ppvar.get():
                    controller.show_frame('Postpro_params')
                else:
                    rnc.FinalPage.params_record(controller.frames['FinalPage'], controller)               #Generate parameter labels
                    controller.show_frame('FinalPage')

        else:
            if controller.ppvar.get():
                controller.show_frame('Postpro_params')
            else:
                rnc.FinalPage.params_record(controller.frames['FinalPage'], controller)               #Generate parameter labels
                controller.show_frame('FinalPage')

    def cv_return(self, controller):
        if controller.metacore['num_cv'] == 2:
            if controller.metacore['collvar_name(2)'].get() == 'SFN':
                prev_page = controller.cvpage_link[controller.metacore['collvar_name(1)'].get()]
                controller.show_frame(prev_page)
            else:
                controller.show_frame('Metadyn_params')

        else:
            controller.show_frame('Metadyn_params')