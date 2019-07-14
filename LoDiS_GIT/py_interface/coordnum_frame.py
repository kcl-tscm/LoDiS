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
***** Fermi Coordination Number Frame *****
Frame containing the Fermi CN variables for metadynamics and MgO/implicit environment runs.

        self.cn = {'n_pwr': tk.IntVar(),                       #n
                   'm_pwr': tk.IntVar(),                       #m
                   'rzero': tk.DoubleVar(),                    # r0 bulk lattice ref
                   'gwidth': tk.DoubleVar()}                   # Guassian width
"""

class CN_params(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        self.title = tk.Label(self, text = 'Coordination Number parameters')
        self.title.grid(row=0, column=0, columnspan = 3, sticky = tk.W, pady = (0,20))
        self.title.config(font=("Courier", 18))

        n_label = tk.Label(self, text = 'n', padx = 10)
        n_entry = tk.Entry(self, textvariable = controller.cn['n_pwr'])
        n_label.grid(row=1, column=0, sticky = tk.W, pady = (0,10))
        n_entry.grid(row=1, column=1, sticky = tk.W)

        m_label = tk.Label(self, text = 'm', padx = 10)
        m_entry = tk.Entry(self, textvariable = controller.cn['m_pwr'])
        m_label.grid(row=1, column=2, sticky = tk.W, pady = (0,10))
        m_entry.grid(row=1, column=3, sticky = tk.W)

        r_label = tk.Label(self, text = 'r0', padx = 10)
        r_entry = tk.Entry(self, textvariable = controller.cn['rzero'])
        r_label.grid(row=2, column=0, sticky = tk.W, pady = (0,10))
        r_entry.grid(row=2, column=1, sticky = tk.W)

        gw_label = tk.Label(self, text = 'gwidth', padx = 10)
        gw_entry = tk.Entry(self, textvariable = controller.cn['gwidth'])
        gw_label.grid(row=2, column=2, sticky = tk.W, pady = (0,10))
        gw_entry.grid(row=2, column=3, sticky = tk.W)

        nextb = tk.Button(self, text = 'Next',padx =30, command = lambda: self.fcheck_proceed(controller))
        nextb.grid(row = 3, column = 4, sticky = tk.E, pady = (30, 0))
        prevb = tk.Button(self, text = 'Back',padx=30, command = lambda: self.cv_return(controller))
        prevb.grid(row = 3, column = 3, sticky = tk.E, pady = (30, 0))

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
            if controller.metacore['collvar_name(2)'].get() != 'coord_number':
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
            if controller.metacore['collvar_name(2)'].get() == 'coord_number':
                prev_page = controller.cvpage_link[controller.metacore['collvar_name(1)'].get()]
                controller.show_frame(prev_page)
            else:
                controller.show_frame('Metadyn_params')

        else:
            controller.show_frame('Metadyn_params')