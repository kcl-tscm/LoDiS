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
***** Quench Frame *****
Quenching parameters.

        self.quench = {'itremp': tk.IntVar(),              #Step at which quench starts, itremp < npas unless microcanonical
                       'tmin': tk.DoubleVar()}           #Min temperature to reach
"""
class Quench_params(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        self.title = tk.Label(self, text = 'Quenching parameters')
        self.title.grid(row=0, column=0, columnspan = 2, sticky = tk.W, pady = (0,20))
        self.title.config(font=("Courier", 18))


        itremp_label = tk.Label(self, text = 'itremp', padx = 10)
        itremp_entry = tk.Entry(self, textvariable = controller.quench['itremp'])
        itremp_label.grid(row=1, column=0, sticky = tk.W, pady = (0,20))
        itremp_entry.grid(row=1, column=1, sticky = tk.W)

        tmin_label = tk.Label(self, text = 'tmin', padx = 10)
        tmin_entry = tk.Entry(self, textvariable = controller.quench['tmin'])
        tmin_label.grid(row=2, column=0, sticky = tk.W, pady = (0,20))
        tmin_entry.grid(row=2, column=1, sticky = tk.W)

        nextb = tk.Button(self, text = 'Next',padx =30, command = lambda: self.fcheck_proceed(controller))
        nextb.grid(row = 3, column = 2, sticky = tk.E, pady = (30, 0))
        prevb = tk.Button(self, text = 'Back',padx=30, command = lambda: controller.show_frame('Sys_params'))
        prevb.grid(row = 3, column = 1, sticky = tk.E, pady = (30, 0))

    def fcheck_proceed(self, controller):
        try:
            controller.quench['itremp'].get()
            if controller.quench['itremp'].get() > 0:
                try:
                    controller.quench['tmin'].get()
                    if controller.quench['tmin'].get() > 0:
                        if controller.ppvar.get():
                            controller.show_frame('Postpro_params')
                        else:
                            rnc.FinalPage.params_record(controller.frames['FinalPage'], controller)               #Generate parameter labels
                            controller.show_frame('FinalPage')
                    else:
                        tkMessageBox.showinfo('Error', 'Final temperature is equal or less than zero')
                except:
                    tkMessageBox.showinfo('Error', '*tmin* input is unknown')
            else:
                tkMessageBox.showinfo('Error', 'Quench cannot start before process')
        except:
            tkMessageBox.showinfo('Error', '*itremp* input is unknown')