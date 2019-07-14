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
import postpro_frame as ppf
import run_check as rnc

"""
***** System Parameters Frame *****
Non-process specific system parameters.

        self.sys_vn = ['elem1', 'elem2', 'natom',
                       'type_potential','tstep', 'npas',
                       'irand','tinit', 'vnu',
                       'npast','scrivo','sticky_atoms',
                       'sticky_k','cn_cal', 'output_xyz', 'impl_env',
                       'pot_a','pot_b', 'eta_a', 'eta_b']
"""

class Sys_params(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        self.title = tk.Label(self, text = 'System parameters')
        self.title.grid(row=0, column=0, columnspan = 2, sticky = tk.W, pady = (0,20))
        self.title.config(font=("Courier", 18))

        v =3                #controller.sys_var index, start from 3 as elem1, elem2 and natom are already known.
        self.k = 1                #line number
        while controller.sys_vn[v] != 'cn_cal':
            if controller.sys_vn[v] == 'type_potential':
                if v%2 == 0:
                    tp_label = tk.Label(self, text = 'type_potential', padx = 10)
                    tp_label.grid(row=self.k, column=0, sticky = tk.E, pady = (0,10))
                    tp_option = tk.OptionMenu(self, controller.sys_var[v], *controller.typot)
                    tp_option.grid(row=self.k, column=1, sticky = tk.W)
                else:
                    tp_label = tk.Label(self, text = 'type_potential', padx = 10)
                    tp_label.grid(row=self.k, column=2, sticky = tk.E, pady = (0,10))
                    tp_option = tk.OptionMenu(self, controller.sys_var[v], *controller.typot)
                    tp_option.grid(row=self.k, column=3, sticky = tk.W)
                    self.k += 1
                v += 1
            else:
                if v%2 == 0:
                    label_v = tk.Label(self, text=controller.sys_vn[v], padx = 10)
                    entry_v = tk.Entry(self,textvariable = controller.sys_var[v])
                    label_v.grid(row=self.k, column = 0, sticky = tk.E, pady = (0, 10))
                    entry_v.grid(row = self.k,sticky = tk.W, column=1)
                else:
                    label_v = tk.Label(self, text=controller.sys_vn[v], padx = 10)
                    entry_v = tk.Entry(self,textvariable = controller.sys_var[v])
                    label_v.grid(row=self.k, column = 2, sticky = tk.E, pady = (0, 10))
                    entry_v.grid(row = self.k,sticky = tk.W, column=3)
                    self.k +=1
                v += 1
        else:                                               #Boolean checks
            self.k +=1
            self.cncheck = tk.Checkbutton(self, text = 'cn_cal', variable = controller.sys_var[v],
                                      anchor = tk.W, padx = 10, command = lambda:self.env_switch(controller))
            self.cncheck.grid(row = self.k, column = 0, pady = 10, sticky = tk.W)
            v +=1

            xyzout = tk.Checkbutton(self, text = 'Output xyz file', variable = controller.sys_var[v],
                                     anchor = tk.W, padx =10)
            xyzout.grid(row = self.k, column = 2, pady = 10, sticky = tk.W)
            v += 1
            self.k += 1
            self.envcheck = tk.Checkbutton(self, text = 'Implicit Environment', variable = controller.sys_var[v],
                                      anchor = tk.W, padx = 10, command = lambda:self.env_switch(controller))
            self.envcheck.grid(row = self.k, column = 0, pady = 10, sticky = tk.W)
            v +=1
            self.ppcheck = tk.Checkbutton(self, text = 'Post-processing', variable = controller.ppvar,anchor = tk.W, padx = 10)
            self.ppcheck.grid(row = self.k, column = 2, pady = 10, sticky = tk.W)
        self.k += 1
        #implicit environment variables, Entry default = DISABLED
        self.impl_vars = []
        for j in range(4):
            if j%2 == 0:
                label_k = tk.Label(self, text=controller.sys_vn[v], padx = 10)
                self.impl_vars.append(tk.Entry(self,textvariable = controller.sys_var[v]))
                label_k.grid(row = self.k, sticky = tk.E, column = 0, pady = (0, 10))
                self.impl_vars[j].grid(row =self.k, column=1, sticky = tk.W)
                self.impl_vars[j].config(state = tk.DISABLED)
            else:
                label_k = tk.Label(self, text=controller.sys_vn[v], padx = 10)
                self.impl_vars.append(tk.Entry(self,textvariable = controller.sys_var[v]))
                label_k.grid(row=self.k , column = 2, sticky = tk.E, pady = (0, 10))
                self.impl_vars[j].grid(row = self.k, column=3, sticky = tk.W)
                self.impl_vars[j].config(state = tk.DISABLED)
                self.k +=1
            v += 1
        self.k += 1
        prevb = tk.Button(self, text = 'Back',padx=30, command = lambda: controller.show_frame('Input_files'))
        prevb.grid(row = self.k+1, column = 3, sticky = tk.E, pady = (30, 0))
        nextb = tk.Button(self, text = 'Next',padx =30, command = lambda: self.process_params(controller))
        nextb.grid(row = self.k+1, column = 4, sticky = tk.E, pady = (30, 0))

    def vacf_create(self, controller):

        vacf_check = tk.Checkbutton(self, text = 'vel_af', variable = controller.vacf_var['vel_af'],
                                    anchor = tk.W, padx = 10)
        vacf_check.grid(row = self.k, column = 0, pady = (0,10), sticky = tk.W)
        vacf_check.config(state= tk.DISABLED)
        if controller.vacf_on:
            vacf_check.config(state = tk.NORMAL)


    def env_switch(self, controller):
        if controller.sys_var[controller.envind].get() == True:
            for i in range(4):
                self.impl_vars[i].config(state=tk.NORMAL)
        else:
            for i in range(4):
                self.impl_vars[i].config(state=tk.DISABLED)

    def growth_envoff(self, controller):
        if controller.run_process == 'Growth':
            controller.sys_var[controller.envind].set(False)
            self.envcheck.config(state = tk.DISABLED)
        else:
            self.envcheck.config(state = tk.NORMAL)

    def cn_switch(self, controller):
        cn_ind = controller.sys_vn.index('cn_cal')
        if controller.bimon.get():
            controller.sys_var[cn_ind].set(False)
            self.cncheck.config(state = tk.DISABLED)
        else:
            self.cncheck.config(state = tk.NORMAL)

    def ppro_compatability(self, controller):
        if controller.run_process in ['Growth', 'Coalescence']:
            controller.ppvar.set(False)
            self.ppcheck.config(state = tk.DISABLED)
        else:
            self.ppcheck.config(state = tk.NORMAL)


    def process_params_choose(self, controller):
        if controller.run_process == controller.label_list[0]:
            if controller.ppvar.get():
                controller.show_frame('Postpro_params')
            else:
                rnc.FinalPage.params_record(controller.frames['FinalPage'], controller)               #Generate parameter labels
                controller.show_frame('FinalPage')
        if controller.run_process == controller.label_list[1]:
            controller.show_frame('IT_params')
        elif controller.run_process == controller.label_list[2]:
            controller.show_frame('Coalescence_params')
        elif controller.run_process == controller.label_list[3]:
            controller.show_frame('Quench_params')
        elif controller.run_process == controller.label_list[4]:
            controller.show_frame('Growth_params')
        elif controller.run_process == controller.label_list[5]:
            controller.show_frame('Metadyn_params')

    def process_params(self, controller):
        for i in range(len(controller.sys_var)):
            try:
                controller.sys_var[i].get()                     #checks whether input is correct variable type
            except:
                tkMessageBox.showinfo('Error', 'Incorrect characters for '+controller.sys_vn[i])
        default_check = ['', ' ', '()', ',']
        if controller.sys_var[0].get() in default_check:
            tkMessageBox.showinfo('Error', 'Element 1 not given')
        elif controller.sys_var[1].get() in default_check:
            tkMessageBox.showinfo('Error', 'Element 2 not given')
        elif controller.sys_var[2].get() == 0:
            tkMessageBox.showinfo('Error', 'Number of atoms not specified')
        elif controller.sys_var[3].get() not in ['rgl', 'lj1', 'gir', 'par']:
            tkMessageBox.showinfo('Error', 'Unknown interaction type chosen')
        elif controller.sys_var[4].get() == 0.0:
            tkMessageBox.showinfo('Error', 'Zero time step given')
        elif controller.sys_var[5].get() == 0:
            tkMessageBox.showinfo('Error', 'Null run selected, please give a simulation period')
        elif controller.sys_var[7].get() == 0.0:
            tkMessageBox.showinfo('Error', 'Initial T(k) must be greater than zero ')
        elif controller.sys_var[8].get() == 0.0 and controller.run_process != 'NVT':
            tkMessageBox.showinfo('Error','Thermostat switched off for non NVE run')
        elif controller.sys_var[9].get() == 0:
            tkMessageBox.showinfo('Error', 'Please provide non zero value for *npast*')
        elif controller.sys_var[10].get() == 0:
            tkMessageBox.showinfo('Error', 'Measurement increments, *scrivo*, not provided')
        elif controller.sys_var[controller.envind].get() == True:
            etaind = controller.sys_vn.index('eta_a')
            if (controller.sys_var[etaind].get() or controller.sys_var[etaind+1].get()) == 0.0:
                tkMessageBox.showinfo('Error', 'Give non zero values for eta a and eta b')
            else:
                controller.sys_out = {controller.sys_vn[s]:controller.sys_var[s] for s in range(len(controller.sys_vn))}
                if controller.ppvar.get():                                                                              #Create widgets for Postpro frame
                    ppf.Postpro_params.create_ppwidgets(controller.frames['Postpro_params'], controller)
                self.process_params_choose(controller)
        else:
            controller.sys_out = {controller.sys_vn[s]:controller.sys_var[s] for s in range(len(controller.sys_vn))}
            if controller.ppvar.get():
                ppf.Postpro_params.create_ppwidgets(controller.frames['Postpro_params'], controller)
            self.process_params_choose(controller)