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
import coal_visual as cvis

"""
***** Coalescence Frame *****
Frame listing the parameters for the Coalescence process.

        self.coal = {'natom2': tk.IntVar(),
                     'somedist': tk.DoubleVar()}
"""
class Coalescence_params(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        self.title = tk.Label(self, text = 'Coalescence parameters')
        self.title.grid(row=0, column=0, columnspan = 2, sticky = tk.W, pady = (0,20))
        self.title.config(font=("Courier", 18))

        self.poslabel2 = tk.Label(self, text='2nd position file (.xyz format)', pady = 20)
        self.poslabel2.grid(row=1, column= 0, sticky=tk.E)
        self.browsepos2 = tk.Button(self,text= 'Browse', command = lambda: self.pos2_search(controller))
        self.browsepos2.grid(row=1, column=1, sticky = tk.W)


        label_sd = tk.Label(self, text='somedist', padx =10)
        entry_sd = tk.Entry(self, textvariable = controller.coal['somedist'])
        label_sd.grid(row=3, column=0, pady = (0,10))
        entry_sd.grid(row=3, column=1, pady= (0,10), sticky = tk.E)

        self.nextb = tk.Button(self, text = 'Next',padx =30, command = lambda: self.fcheck_proceed(controller))
        self.nextb.grid(row = 5, column = 2, sticky = tk.E, pady = (30, 0))
        self.prevb = tk.Button(self, text = 'Back',padx=30, command = lambda: controller.show_frame('Sys_params'))
        self.prevb.grid(row = 5, column = 1, sticky = tk.E, pady = (30, 0))


    def pos2_search(self, controller):
        controller.infiles['posfile2'] = tkFileDialog.askopenfilename(initialdir = "/home/",title = "Select file",filetypes = (("xyz files","*.xyz"),("all files","*.*")))
        self.pos2l = tk.Label(self, text = controller.infiles['posfile2'])
        self.pos2l.grid(row = 2, column =0, columnspan = 2, sticky = tk.W)
        self.pos2l.config(font=("Courier", 8))
        self.c_elem1, self.c_elem2, natom2, self.err_out = cvis.pos2_visual(controller.infiles['posfile2'], controller.atom_rad, controller.sys_out['elem1'].get(), controller.sys_out['elem2'].get())
        print ('%s, %s'%(natom2, self.err_out))
        controller.coal['natom2'].set(natom2)

    def fcheck_proceed(self, controller):
        try:
            controller.infiles['posfile2']
        except:
            tkMessageBox.showinfo('Error', '2nd cluster position file has not been provided')
            return

        if os.path.isfile(controller.infiles['posfile2']) and controller.infiles['posfile2'].endswith('.xyz'):        #Checks .xyz file was selected
            try:
                controller.coal['natom2'].get()
            except:
                tkMessageBox.showinfo('Error', 'Null or Non-integer value for natom2')
                return

            if controller.coal['natom2'].get() > 0:
                try:
                    controller.coal['somedist'].get()
                except:
                    tkMessageBox.showinfo('Error', 'Separation distance input is unknown')
                    return

                if controller.coal['somedist'].get() > 0:
                    if self.err_out == 0:
                        c_elem = [self.c_elem1, self.c_elem2]
                        f_elem = [controller.sys_out['elem1'].get(), controller.sys_out['elem2'].get()]
                        e_list = controller.atom_rad.keys()
                        e_list.sort()
                        if c_elem != f_elem and controller.bimon.get() == False:          # Cluster 1 = [e1, e1], Cluster 2 = [e2, e2] or [e1, e2]
                            controller.sys_out['elem1'].set(e_list[0])
                            controller.sys_out['elem2'].set(e_list[1])
                            controller.bimon.set(True)

                        if controller.ppvar.get():
                            controller.show_frame('Postpro_params')
                        else:
                            rnc.FinalPage.params_record(controller.frames['FinalPage'], controller)               #Generate parameter labels
                            controller.show_frame('FinalPage')
                    else:
                        tkMessageBox.showinfo('Error', self.err_out)
                else:
                    tkMessageBox.showinfo('Error', 'Separation distance is equal or less than zero')

            else:
                tkMessageBox.showinfo('Error', 'No atoms in the 2nd cluster')
        else:
            tkMessageBox.showinfo('Error', 'File' + controller.infiles['posfile2'] + 'does not exist')
