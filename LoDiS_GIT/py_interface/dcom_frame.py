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
***** Squared Distance of COM Frame *****
D**2 COM parameters for metadynamics.

       self.d_gwidth = tk.DoubleVar()                                     # Gaussian width
"""

class Dcom_params(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        self.title = tk.Label(self, text = 'Squared distance of COM parameter')
        self.title.grid(row=0, column=0, columnspan = 2, sticky = tk.W, pady = (0,20))
        self.title.config(font=("Courier", 18))

        gw_label = tk.Label(self, text = 'gwidth', padx = 10)
        gw_entry = tk.Entry(self, textvariable = controller.dcom['d_gwidth'])
        gw_label.grid(row=1, column=0, sticky = tk.E, pady = (0,10))
        gw_entry.grid(row=1, column=1, sticky = tk.W)

        nextb = tk.Button(self, text = 'Next',padx =30, command = lambda: self.fcheck_proceed(controller))
        nextb.grid(row = 2, column = 2, sticky = tk.E, pady = (30, 0))
        prevb = tk.Button(self, text = 'Back',padx=30, command = lambda: self.cv_return(controller))
        prevb.grid(row = 2, column = 1, sticky = tk.E, pady = (30, 0))

    def fcheck_proceed(self, controller):
        try:
            controller.dcom['d_gwidth'].get()
        except:
            tkMessageBox.showinfo('Error', 'gwidth contains unknown character')
        if controller.dcom['d_gwidth'].get() <= 0 :
            tkMessageBox.showinfo('Error', 'gwidth must be larger than zero')
            return

        if controller.metacore['num_cv'] == 2:
            if controller.metacore['collvar_name(2)'].get() != 'd_com':
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
            if controller.metacore['collvar_name(2)'].get() == 'd_com':
                prev_page = controller.cvpage_link[controller.metacore['collvar_name(1)'].get()]
                controller.show_frame(prev_page)
            else:
                controller.show_frame('Metadyn_params')

        else:
            controller.show_frame('Metadyn_params')