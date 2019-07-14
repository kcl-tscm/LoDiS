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
import read_files as rf
import system_params_gui as spg

"""
***** Process Type Frame *****
A frame that lists a set of Checkboxes corresponding to each process available to LoDiS.
Selecting a process will disable the others boxes, must deselect your choice to choose another.
If NVT is chosen, vel_af checkboxes in System_params frame is activated.
"""
class StartPage(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        self.title = tk.Label(self, text = 'Which process do you wish to run?')
        self.title.grid(row=0, column=0, columnspan = 2, sticky = tk.W, pady = (0,20))
        self.title.config(font=("Courier", 18))
        self.varl = [tk.BooleanVar() for l in range(len(controller.label_list))]
        for l in self.varl:
            l.set(False)
        self.checks = []
        for l in range(len(controller.label_list)):
            self.checks.append(tk.Checkbutton(self, text=controller.label_list[l], variable = self.varl[l], anchor= tk.N,
                                 padx = 30, pady = 10,command = lambda: self.check_on()))
            self.checks[l].grid(row=l+1, column = 0, sticky = tk.W)

        nextb = tk.Button(self, text='Next', padx = 30, command = lambda: self.inputf_proceed(controller))
        nextb.grid(row=len(controller.label_list) +2, column = 1, sticky = tk.E, pady = (20, 0))

    def check_on(self):                              #Diables other buttons when one process is selected
        if True in [l.get() for l in self.varl]:
            for l in range(len(self.varl)):
                self.checks[l].config(state = tk.DISABLED if not self.varl[l].get() else tk.NORMAL)
        else:
            for l in range(len(self.varl)):
                self.checks[l].config(state = tk.NORMAL)

    def inputf_proceed(self, controller):
        varl_get = [l.get() for l in self.varl]
        if True in varl_get:
            pro_idx = varl_get.index(True)
            controller.run_process = controller.label_list[pro_idx]
            if controller.run_process == 'NVT':
                controller.vacf_on = True
                spg.Sys_params.vacf_create(controller.frames['Sys_params'], controller)
            else:
                controller.vacf_on = False
                spg.Sys_params.vacf_create(controller.frames['Sys_params'], controller)

            rf.Input_files.mgo_compatability(controller.frames['Input_files'], controller)     # Lock use of substrate for Growth and Coalescence
            spg.Sys_params.ppro_compatability(controller.frames['Sys_params'], controller)     # Switch off post-analysis calculations for Growth/Coalescence
            spg.Sys_params.growth_envoff(controller.frames['Sys_params'], controller)           # Switch off implicit environment for Growth
            controller.show_frame('Input_files')
        else:
            tkMessageBox.showinfo('Error', 'You have not selected a process to run')
