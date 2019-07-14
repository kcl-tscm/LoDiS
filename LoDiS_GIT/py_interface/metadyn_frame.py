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

"""
***** Metadynamics Frame *****
Non_CV specific metadynamics parameters.

        self.cvlist = ['coord_number', 'C2N', 'SFN', 'CN_bim', 'd_com']

        self.metacore = {'collvar_wanted': tk.BooleanVar(),         #Write CV values and average forces at each time step for whichever procedure
                         'metaframe': tk.IntVar(),                             #Measuremnt step increments for movie.xyz
                         'gheight': tk.DoubleVar(),                            #Gaussian height (eV)
                         'metaperiod': tk.IntVar(),                            #After how many step a new Gaussian is added
                         'num_cv': tk.IntVar(),                                # Number of CVs used (1 or 2)
                         'collvar_name(1)': tk.StringVar(),                      #'coord_number', 'C2N', 'SFN', 'CN_bim', 'd_com', 'none'
                         'collvar_name(2)': tk.StringVar()}
"""

class Metadyn_params(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        self.title = tk.Label(self, text = 'Metadynmaics core parameters')
        self.title.grid(row=0, column=0, columnspan = 2, sticky = tk.W, pady = (0,20))
        self.title.config(font=("Courier", 18))

        cv_check = tk.Checkbutton(self, text = 'collvar_wanted', variable = controller.metacore['collvar_wanted'],
                                      anchor = tk.W, padx = 10)
        cv_check.grid(row=1, column=0, sticky = tk.W, pady= (0,20))

        cv1_label = tk.Label(self, text = 'collvar_name(1)', padx = 10)
        cv1_label.grid(row=2, column=0, sticky = tk.E, pady = (0,10))
        cv1_option = tk.OptionMenu(self, controller.metacore['collvar_name(1)'], *controller.cvlist)
        cv1_option.grid(row=2, column=1, sticky = tk.W)

        cv2_label = tk.Label(self, text = 'collvar_name(2)', padx = 10)
        cv2_label.grid(row=2, column=2, sticky = tk.E, pady = (0,10))
        cv2_option = tk.OptionMenu(self, controller.metacore['collvar_name(2)'], *(controller.cvlist + ['none']))
        cv2_option.grid(row=2, column=3, sticky = tk.W)

        mf_label = tk.Label(self, text = 'metaframe', padx = 10)
        mf_entry = tk.Entry(self, textvariable = controller.metacore['metaframe'])
        mf_label.grid(row=3, column=0, sticky = tk.E, pady = (0,10))
        mf_entry.grid(row=3, column=1, sticky = tk.W)

        gh_label = tk.Label(self, text = 'gheight', padx = 10)
        gh_entry = tk.Entry(self, textvariable = controller.metacore['gheight'])
        gh_label.grid(row=3, column=2, sticky = tk.E, pady = (0,10))
        gh_entry.grid(row=3, column=3, sticky = tk.W)

        mp_label = tk.Label(self, text = 'metaperiod', padx = 10)
        mp_entry = tk.Entry(self, textvariable = controller.metacore['metaperiod'])
        mp_label.grid(row=4, column=0, sticky = tk.E, pady = (0,10))
        mp_entry.grid(row=4, column=1, sticky = tk.W)


        nextb = tk.Button(self, text = 'Next',padx =30, command = lambda: self.fcheck_proceed(controller))
        nextb.grid(row = 5, column = 4, sticky = tk.E, pady = (30, 0))
        prevb = tk.Button(self, text = 'Back',padx=30, command = lambda: controller.show_frame('Sys_params'))
        prevb.grid(row = 5, column = 3, sticky = tk.E, pady = (30, 0))

    def fcheck_proceed(self, controller):
        for m in controller.metacore.keys():
            if m != 'num_cv':
                try:
                    controller.metacore[m].get()
                except:
                    tkMessageBox.showinfo('Error', '*'+m+'* contains unknown character')
        if controller.metacore['collvar_name(1)'].get() == controller.metacore['collvar_name(2)'].get():
            tkMessageBox.showinfo('Error', 'collvar_name(1) and collvar_name(2) cannot be the same!')
        elif controller.metacore['metaframe'].get() <= 0:
            tkMessageBox.showinfo('Error', 'Measurement steps must be larger than zero')
        elif controller.metacore['gheight'].get() <= 0:
            tkMessageBox.showinfo('Error', 'Guassian height must be larger than zero')
        elif controller.metacore['metaperiod'].get() == 0:
            tkMessageBox.showinfo('Error', 'metaperiod must be greater than zero')
        else:
            if controller.metacore['collvar_name(2)'].get() == 'none':
                controller.metacore['num_cv'] = 1
            else:
                controller.metacore['num_cv'] = 2

            next_page = controller.cvpage_link[ controller.metacore['collvar_name(1)'].get()]
            controller.show_frame(next_page)


