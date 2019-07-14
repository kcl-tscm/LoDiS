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
import system_params_gui as spg
import cluster_repos_prep as crp
import visualise_cluster as vc


"""
***** Select Input Files Frames *****
Browse input .xyz and .pot files. Option to simulate a supported cluster is present here.

        self.infiles = {'posfile':'N/a', 'potfile':'N/a', 'mgofile':'N/a', 'posfile2':'N/a'}
"""
class Input_files(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        self.title = tk.Label(self, text = 'Provide the nescessary input files for simulation')
        self.title.grid(row=0, column=0, columnspan = 2, sticky = tk.W, pady = (0,20))
        self.title.config(font=("Courier", 18))
        self.poslabel = tk.Label(self, text='Initial positions file (.xyz format)', pady = 20)
        self.poslabel.grid(row=1, column= 0, sticky=tk.W)
        self.browsepos = tk.Button(self,text= 'Browse', command = lambda: self.pos_search(controller))
        self.browsepos.grid(row=1, column=1, sticky = tk.W)

        self.potlabel = tk.Label(self, text='Potential parameters file (.pot format)', pady =20)
        self.potlabel.grid(row=3, column=0,  sticky=tk.W)
        self.browsepot = tk.Button(self,text= 'Browse', command = lambda: self.pot_search(controller))
        self.browsepot.grid(row=3, column=1, sticky = tk.W)


        self.mgocheck = tk.Checkbutton(self, text='Simulate MgO substrate at z=0?', variable= controller.varmgo,
                                  anchor = tk.W, command = lambda: self.mgo_switch(controller))
        self.mgocheck.grid(row=5, column=0, sticky = tk.W, pady = (20,0))
        self.mgolabel = tk.Label(self, text= 'MgO parameters file (.pot format)')
        self.mgolabel.grid(row=6, column =0, sticky = tk.W)
        self.mgobrowse = tk.Button(self,text= 'Browse', command = lambda: self.mgo_search(controller))
        self.mgobrowse.grid(row = 6, column=1, sticky = tk.W)
        self.mgobrowse.config(state = tk.DISABLED)

        self.nextb = tk.Button(self, text = 'Next',padx =30, command = lambda: self.params_proceed(controller))
        self.nextb.grid(row = 8, column = 2, sticky = tk.E, pady = (30, 0))
        self.prevb = tk.Button(self, text = 'Back',padx=30, command = lambda: controller.show_frame('StartPage'))
        self.prevb.grid(row = 8, column = 1, sticky = tk.E, pady = (30, 0))

    def pos_search(self, controller):
        controller.infiles['posfile'] = tkFileDialog.askopenfilename(initialdir = "/home/",title = "Select file",
                                                                     filetypes = (("xyz files","*.xyz"),("all files","*.*")))
        self.pospath = tk.Label(self, text = controller.infiles['posfile'])
        self.pospath.grid(row = 2, column =0, columnspan = 2, sticky = tk.W)
        self.pospath.config(font=("Courier", 8))

    def pot_search(self, controller):
        controller.infiles['potfile'] = tkFileDialog.askopenfilename(initialdir = "/home/",title = "Select file",filetypes = (("pot files","*.pot"),("all files","*.*")))
        self.potpath = tk.Label(self, text = controller.infiles['potfile'])
        self.potpath.grid(row = 4, column =0, columnspan = 2, sticky = tk.W)
        self.potpath.config(font=("Courier", 8))

    def mgo_switch(self, controller):
        if controller.varmgo.get() == True:
            self.mgobrowse.config(state=tk.NORMAL)
        else:
            self.mgobrowse.config(state=tk.DISABLED)

    def mgo_compatability(self, controller):
        if controller.run_process in ['Growth', 'Coalescence']:
            controller.varmgo.set(False)
            self.mgocheck.config(state = tk.DISABLED)
        else:
            self.mgocheck.config(state = tk.NORMAL)

    def mgo_search(self, controller):
        controller.infiles['mgofile'] = tkFileDialog.askopenfilename(initialdir = "/home/",title = "Select file",filetypes = (("pot files","*.pot"),("all files","*.*")))
        self.mgopath = tk.Label(self, text = controller.infiles['mgofile'])
        self.mgopath.grid(row = 7, column =0, columnspan = 2, sticky = tk.W)
        self.mgopath.config(font=("Courier", 8))

    def params_proceed(self, controller):
        try:
            controller.infiles['posfile']
            if os.path.isfile(controller.infiles['posfile']) and controller.infiles['posfile'].endswith('.xyz'):        #Checks .xyz file was selected
                postrue = True
            else:
                postrue = False
        except:
            postrue = False
        try:
            controller.infiles['potfile']
            if os.path.isfile(controller.infiles['potfile']) and controller.infiles['potfile'].endswith('.pot'):
                pottrue = True
            else:
                pottrue = False
        except:
            pottrue = False

        if postrue == True and pottrue == True and controller.varmgo.get() == False:
            f_elem1, f_elem2, numatom, atom_rad, err_out, atom_mass = crp.cluster_repos(controller.infiles['posfile'],'None', controller.infiles['potfile'],
            controller.varmgo.get(), controller.run_process)    #Get elem1, elem2, N
            if err_out == 0:
                controller.sys_var[0].set(f_elem1)            #sets values of elem1, elem2 and natom to values read from .xyz file
                controller.sys_var[1].set(f_elem2)
                controller.sys_var[2].set(numatom)
                controller.atom_rad = atom_rad
                controller.atom_mass = atom_mass
                if f_elem1 == f_elem2:                       # Bimetallic ?
                    controller.bimon.set(False)
                else:
                    controller.bimon.set(True)

                spg.Sys_params.cn_switch(controller.frames['Sys_params'], controller)                 # Turn off cn_cal for bimetallic clusters
                controller.show_frame('Sys_params')            #Run before cluster_visual so the latter appears on top.
                vc.cluster_visual(controller.infiles['posfile'], controller.atom_rad, controller.varmgo.get(), controller.bimon.get())         #show cluster
            else:
                tkMessageBox.showinfo('Error', err_out[1])

        elif postrue == True and pottrue == True and controller.varmgo.get() == True:
            try:
                controller.infiles['mgofile']
                if os.path.isfile(controller.infiles['mgofile']) and controller.infiles['mgofile'].endswith('.pot'):
                    # Viualise and repositioning of cluster should be run here!
                    f_elem1, f_elem2, numatom, atom_rad, err_out, atom_mass = crp.cluster_repos(controller.infiles['posfile'], controller.infiles['mgofile'], controller.infiles['potfile'],
                    controller.varmgo.get(), controller.run_process)    #Repositions cluster, elem1, elem2, N
                    if err_out == 0:
                        controller.sys_var[0].set(f_elem1)            #sets values of elem1, elem2 and natom to values read from .xyz file
                        controller.sys_var[1].set(f_elem2)
                        controller.sys_var[2].set(numatom)
                        controller.atom_rad = atom_rad
                        controller.atom_mass = atom_mass
                        if f_elem1 == f_elem2:
                            controller.bimon.set(False)
                        else:
                            controller.bimon.set(True)

                        if '_repos.xyz' not in controller.infiles['posfile']:
                            controller.infiles['posfile'] = controller.infiles['posfile'].replace('.xyz', '_repos.xyz')           #sets the NEW repos file as the new position file

                        controller.show_frame('Sys_params')            #Run before cluster_visual so the latter appears on top.
                        vc.cluster_visual(controller.infiles['posfile'], controller.atom_rad, controller.varmgo.get(), controller.bimon.get())            #Plot visual of cluster using matplotlib
                    else:
                        tkMessageBox.showinfo('Error', err_out[1])
						
                else:
                    tkMessageBox.showinfo('Error', 'Have not selected a pot file containing the substrate parameters')
            except:
                tkMessageBox.showinfo('Error', 'Have not selected a pot file containing the substrate parameters')
        elif postrue == True and pottrue == False:
            tkMessageBox.showinfo('Error', 'Have not selected a pot file containing the interaction parameters')
        elif postrue == False and pottrue == True:
            tkMessageBox.showinfo('Error', 'Have not selected a xyz file with atom positions')
        else:
            tkMessageBox.showinfo('Error', 'Insufficient files to run simulation')
