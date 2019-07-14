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
***** Iterative Temperature MD Frame *****
Parameters for itMD runs.
For melting **tcaloric** must be larger than **tinit**, and vice versa for freezing.

        self.calor = {'deltat':tk.DoubleVar(),            # Temperature step
                      'tcaloric':tk.DoubleVar()}            # Final temperature
"""
class IT_params(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        self.title = tk.Label(self, text = 'Iterative Temperature MD parameters')
        self.title.grid(row=0, column=0, columnspan = 2, sticky = tk.W, pady = (0,20))
        self.title.config(font=("Courier", 18))

        descrip = tk.Label(self, text = 'For temperature increase, set *tcaloric* higher than *tinit*.\n'+
                                        'For temperature decrease, do the opposite')
        descrip.grid(row=1, column=0, columnspan = 2, sticky = tk.W, pady = (0,20))

        del_label = tk.Label(self, text = 'deltat', padx = 10)
        del_entry = tk.Entry(self, textvariable = controller.calor['deltat'])
        del_label.grid(row=2, column=0, sticky = tk.W, pady = (0,20))
        del_entry.grid(row=2, column=1, sticky = tk.W)

        tcal_label = tk.Label(self, text = 'tcaloric', padx = 10)
        tcal_entry = tk.Entry(self, textvariable = controller.calor['tcaloric'])
        tcal_label.grid(row=3, column=0, sticky = tk.W, pady = (0,20))
        tcal_entry.grid(row=3, column=1, sticky = tk.W)

        nextb = tk.Button(self, text = 'Next',padx =30, command = lambda: self.fcheck_proceed(controller))
        nextb.grid(row = 4, column = 2, sticky = tk.E, pady = (30, 0))
        prevb = tk.Button(self, text = 'Back',padx=30, command = lambda: controller.show_frame('Sys_params'))
        prevb.grid(row = 4, column = 1, sticky = tk.E, pady = (30, 0))

    def fcheck_proceed(self, controller):
        for m in controller.calor.keys():
            try:
                controller.calor[m].get()
            except:
                tkMessageBox.showinfo('Error', '*'+m+'* contains unknown character')
                return
            if controller.calor[m].get() <= 0 :
                tkMessageBox.showinfo('Error', '*'+m+'* must be larger than zero')
            else:
                if controller.ppvar.get():
                    controller.show_frame('Postpro_params')
                else:
                    rnc.FinalPage.params_record(controller.frames['FinalPage'], controller)               #Generate parameter labels
                    controller.show_frame('FinalPage')

