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
from PIL import Image

"""
***** Input Listings Frame *****
Creates a canvas listing the parameter values selected by the user as a final check.
Note: Booleans are given as either 0 (False) or 1 (True).
"""
class FinalPage(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller

        self.canvas = tk.Canvas(self, background="white", width=850,height=400,scrollregion=(0,0,600,1100))
        self.scrollf = tk.Canvas(self.canvas, background="white")
        self.vsb = tk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(yscrollcommand=self.vsb.set)

        self.vsb.pack(side="right", fill="y")
        self.canvas.pack(side="left", fill="both", expand=True)
        self.canvas.create_window((0,0), window=self.scrollf, anchor="nw",
                                  tags="self.scrollf")

        #self.scrollf.bind("<Configure>", self.onFrameConfigure)


    def params_record(self, controller):

        self.title = tk.Label(self.scrollf, text = 'Simulation Overview', bg = 'white')
        self.title.grid(row=0, column=0, columnspan = 4, sticky = tk.W, pady = (0,20))
        self.title.config(font=("Courier", 18))

        descrip = 'A '+ controller.run_process+' simulation will be carried out in the background with the following inputs:'
        ov_label = tk.Label(self.scrollf, text = descrip , padx = 10, bg='white')
        ov_label.grid(row=1, column=0, columnspan = 3, sticky = tk.W, pady = (0,10))

        #Input files
        inp_label = tk.Label(self.scrollf, text = 'Input files' , padx = 10, bg='white')
        inp_label.grid(row=2, column=0, columnspan = 3, sticky = tk.W, pady = (0,10))
        inp_label.config(font=("Courier", 14))

        pos_label = tk.Label(self.scrollf, text = 'posfile:' , padx = 10, bg='white')
        pos_label.grid(row=3, column=0, sticky = tk.W, pady = (0,10))
        pos_val = tk.Label(self.scrollf, text = controller.infiles['posfile'] , padx = 10, bg='white')
        pos_val.grid(row=3, column=1, columnspan = 3, sticky = tk.W, pady = (0,10))

        pot_label = tk.Label(self.scrollf, text = 'potfile:' , padx = 10, bg='white')
        pot_label.grid(row=4, column=0, sticky = tk.W, pady = (0,10))
        pot_val = tk.Label(self.scrollf, text = controller.infiles['potfile'] , padx = 10, bg='white')
        pot_val.grid(row=4, column=1, columnspan =3, sticky = tk.W, pady = (0,10))

        rn = 5
        if controller.varmgo.get():
            mgo_label = tk.Label(self.scrollf, text = 'mgofile:' , padx = 10, bg='white')
            mgo_label.grid(row=rn, column=0, sticky = tk.W, pady = (0,10))
            mgo_val = tk.Label(self.scrollf, text = controller.infiles['mgofile'] , padx = 10, bg='white')
            mgo_val.grid(row=rn, column=1, columnspan = 3, sticky = tk.W, pady = (0,10))
            rn +=1
        if controller.run_process == 'Coalescence':
            pos2_label = tk.Label(self.scrollf, text = 'posfile2:' , padx = 10, bg='white')
            pos2_label.grid(row=rn, column=0, sticky = tk.W, pady = (0,10))
            pos2_val = tk.Label(self.scrollf, text = controller.infiles['posfile2'] , padx = 10, bg='white')
            pos2_val.grid(row=rn, column=1, columnspan=3, sticky = tk.W, pady = (0,10))
            rn +=1

        #System Parameters
        if controller.bimon.get() in [1, True]:
            sys_label = tk.Label(self.scrollf, text = 'System - Bimetallic' , padx = 10, bg='white')
            sys_label.grid(row=rn, column=0, columnspan=3, sticky = tk.W, pady = (20,10))
            sys_label.config(font=("Courier", 14))
            rn +=1
        else:
            sys_label = tk.Label(self.scrollf, text = 'System - Monometallic' , padx = 10, bg='white')
            sys_label.grid(row=rn, column=0, columnspan=3, sticky = tk.W, pady = (20,10))
            sys_label.config(font=("Courier", 14))
            rn +=1

        el1_label = tk.Label(self.scrollf, text = 'elem1' , padx = 10, bg='white')
        el1_label.grid(row=rn, column=0, sticky = tk.W, pady = (0,10))
        el1_val = tk.Label(self.scrollf, text = controller.sys_out['elem1'].get() , padx = 10, bg='white')
        el1_val.grid(row=rn, column=1, sticky = tk.W, pady = (0,10))

        el2_label = tk.Label(self.scrollf, text = 'elem2' , padx = 10, bg='white')
        el2_label.grid(row=rn, column=2, sticky = tk.W, pady = (0,10))
        el2_val = tk.Label(self.scrollf, text = controller.sys_out['elem2'].get() , padx = 10, bg='white')
        el2_val.grid(row=rn, column=3, sticky = tk.W, pady = (0,10))
        rn +=1

        i = 0
        for key in controller.sys_vn[2:]:
            if key != 'impl_env':                  # Doesn't write value for impl, instead writes env values if impl_env =1
                if controller.sys_out['impl_env'].get() in [0,False] and i == controller.envind +1:
                    break
                else:
                    if i%2 ==0:
                        sysi_label = tk.Label(self.scrollf, text = key , padx = 10, bg='white')
                        sysi_label.grid(row=rn, column=0, sticky = tk.W, pady = (0,10))
                        sysi_val = tk.Label(self.scrollf, text = controller.sys_out[key].get() , padx = 10, bg='white')
                        sysi_val.grid(row=rn, column=1, sticky = tk.W, pady = (0,10))
                    else:
                        sysi_label = tk.Label(self.scrollf, text = key , padx = 10, bg='white')
                        sysi_label.grid(row=rn, column=2, sticky = tk.W, pady = (0,10))
                        sysi_val = tk.Label(self.scrollf, text = controller.sys_out[key].get() , padx = 10, bg='white')
                        sysi_val.grid(row=rn, column=3, sticky = tk.W, pady = (0,10))
                        rn +=1
            i += 1
        rn +=1

        pro_name = tk.Label(self.scrollf, text = controller.run_process , padx = 10, bg='white')
        pro_name.grid(row=rn, column=0, columnspan=3, sticky = tk.W, pady = (20,10))
        pro_name.config(font=("Courier", 14))
        rn +=1
        if controller.run_process not in  ['Metadynamics', 'Growth']:

            if controller.run_process == 'NVT':
                pro_vars = controller.vacf_var
            elif controller.run_process == 'itMD':
                pro_vars = controller.calor
            elif controller.run_process == 'Coalescence':
                pro_vars = controller.coal
            elif controller.run_process == 'Quenching':
                pro_vars = controller.quench
            else:
                print ('Error unknown process detected')
                sys.exit()
            k = 0
            for key,value in pro_vars.iteritems():
                if k%2 == 0:
                    pro_label = tk.Label(self.scrollf, text = key , padx = 10, bg='white')
                    pro_label.grid(row=rn, column=0, sticky = tk.W, pady = (0,10))
                    pro_val = tk.Label(self.scrollf, text = value.get() , padx = 10, bg='white')
                    pro_val.grid(row=rn, column=1, sticky = tk.W, pady = (0,10))
                    k +=1
                else:
                    pro_label = tk.Label(self.scrollf, text = key , padx = 10, bg='white')
                    pro_label.grid(row=rn, column=2, sticky = tk.W, pady = (0,10))
                    pro_val = tk.Label(self.scrollf, text = value.get() , padx = 10, bg='white')
                    pro_val.grid(row=rn, column=3, sticky = tk.W, pady = (0,10))
                    k +=1
                    rn +=1
            rn += 1
        elif controller.run_process == 'Metadynamics':
            k =0
            for key,value in controller.metacore.iteritems():                   #core metadyn values
                if key != 'num_cv':
                    if k%2 == 0:
                        pro_label = tk.Label(self.scrollf, text = key , padx = 10, bg='white')
                        pro_label.grid(row=rn, column=0, sticky = tk.W, pady = (0,10))
                        pro_val = tk.Label(self.scrollf, text = value.get() , padx = 10, bg='white')
                        pro_val.grid(row=rn, column=1, sticky = tk.W, pady = (0,10))
                        k +=1
                    else:
                        pro_label = tk.Label(self.scrollf, text = key , padx = 10, bg='white')
                        pro_label.grid(row=rn, column=2, sticky = tk.W, pady = (0,10))
                        pro_val = tk.Label(self.scrollf, text = value.get() , padx = 10, bg='white')
                        pro_val.grid(row=rn, column=3, sticky = tk.W, pady = (0,10))
                        k += 1
                        rn +=1
                else:
                    if k%2 == 0:
                        pro_label = tk.Label(self.scrollf, text = key , padx = 10, bg='white')
                        pro_label.grid(row=rn, column=0, sticky = tk.W, pady = (0,10))
                        pro_val = tk.Label(self.scrollf, text = value , padx = 10, bg='white')
                        pro_val.grid(row=rn, column=1, sticky = tk.W, pady = (0,10))
                        k +=1
                    else:
                        pro_label = tk.Label(self.scrollf, text = key , padx = 10, bg='white')
                        pro_label.grid(row=rn, column=2, sticky = tk.W, pady = (0,10))
                        pro_val = tk.Label(self.scrollf, text = value, padx = 10, bg='white')
                        pro_val.grid(row=rn, column=3, sticky = tk.W, pady = (0,10))
                        k += 1
                        rn +=1
            rn += 1

            #Get CV1 values
            n1 = controller.metacore['collvar_name(1)'].get()            #CV1 name
            if controller.metacore['collvar_name(1)'].get() == 'SFN' and controller.c2n_on.get() in [1, 1.0, True]:
                cv1_label = tk.Label(self.scrollf, text = 'C2N' , padx = 10, bg='white')
            else:
                cv1_label = tk.Label(self.scrollf, text = n1 , padx = 10, bg='white')
            cv1_label.grid(row=rn, column=0, sticky = tk.W, pady = (20,10))
            cv1_label.config(font=("Courier", 14))
            rn +=1
            k=0
            for key,value in controller.cv_dict[n1].iteritems():
                if k%2 == 0:
                    pro_label = tk.Label(self.scrollf, text = key , padx = 10, bg='white')
                    pro_label.grid(row=rn, column=0, sticky = tk.W, pady = (0,10))
                    pro_val = tk.Label(self.scrollf, text = value.get() , padx = 10, bg='white')
                    pro_val.grid(row=rn, column=1, sticky = tk.W, pady = (0,10))
                    k +=1
                else:
                    pro_label = tk.Label(self.scrollf, text = key , padx = 10, bg='white')
                    pro_label.grid(row=rn, column=2, sticky = tk.W, pady = (0,10))
                    pro_val = tk.Label(self.scrollf, text = value.get() , padx = 10, bg='white')
                    pro_val.grid(row=rn, column=3, sticky = tk.W, pady = (0,10))
                    k += 1
                    rn +=1
            rn += 1
            if controller.metacore['num_cv'] ==2:                      #If CV2 isn't 'none' then:
                n2 = controller.metacore['collvar_name(2)'].get()
                if controller.metacore['collvar_name(2)'].get() == 'SFN' and controller.c2n_on.get() in [1, 1.0, True]:
                    cv2_label = tk.Label(self.scrollf, text = 'C2N' , padx = 10, bg='white')
                else:
                    cv2_label = tk.Label(self.scrollf, text = n2 , padx = 10, bg='white')
                cv2_label.grid(row=rn, column=0, sticky = tk.W, pady = (20,10))
                cv2_label.config(font=("Courier", 14))
                rn +=1
                k=0
                for key,value in controller.cv_dict[n2].iteritems():
                    if k%2 == 0:
                        pro_label = tk.Label(self.scrollf, text = key , padx = 10, bg='white')
                        pro_label.grid(row=rn, column=0, sticky = tk.W, pady = (0,10))
                        pro_val = tk.Label(self.scrollf, text = value.get() , padx = 10, bg='white')
                        pro_val.grid(row=rn, column=1, sticky = tk.W, pady = (0,10))
                        k +=1
                    else:
                        pro_label = tk.Label(self.scrollf, text = key , padx = 10, bg='white')
                        pro_label.grid(row=rn, column=2, sticky = tk.W, pady = (0,10))
                        pro_val = tk.Label(self.scrollf, text = value.get() , padx = 10, bg='white')
                        pro_val.grid(row=rn, column=3, sticky = tk.W, pady = (0,10))
                        k += 1
                        rn +=1
                rn += 1
        elif controller.run_process == 'Growth':
            controller.lcs_ind = [l.get() for l in controller.growth['lcs']].index(True) + 1     #lcs = 1, 2,3
            k =0
            for key,value in controller.growth.iteritems():
                if key != 'lcs':
                    if k%2 == 0:
                        pro_label = tk.Label(self.scrollf, text = key , padx = 10, bg='white')
                        pro_label.grid(row=rn, column=0, sticky = tk.W, pady = (0,10))
                        pro_val = tk.Label(self.scrollf, text = value.get() , padx = 10, bg='white')
                        pro_val.grid(row=rn, column=1, sticky = tk.W, pady = (0,10))
                        k +=1
                    else:
                        pro_label = tk.Label(self.scrollf, text = key , padx = 10, bg='white')
                        pro_label.grid(row=rn, column=2, sticky = tk.W, pady = (0,10))
                        pro_val = tk.Label(self.scrollf, text = value.get() , padx = 10, bg='white')
                        pro_val.grid(row=rn, column=3, sticky = tk.W, pady = (0,10))
                        k += 1
                        rn +=1
                else:
                    if k%2 == 0:
                        pro_label = tk.Label(self.scrollf, text = key , padx = 10, bg='white')
                        pro_label.grid(row=rn, column=0, sticky = tk.W, pady = (0,10))
                        pro_val = tk.Label(self.scrollf, text = controller.lcs_ind , padx = 10, bg='white')
                        pro_val.grid(row=rn, column=1, sticky = tk.W, pady = (0,10))
                        k +=1
                    else:
                        pro_label = tk.Label(self.scrollf, text = key , padx = 10, bg='white')
                        pro_label.grid(row=rn, column=2, sticky = tk.W, pady = (0,10))
                        pro_val = tk.Label(self.scrollf, text = controller.lcs_ind, padx = 10, bg='white')
                        pro_val.grid(row=rn, column=3, sticky = tk.W, pady = (0,10))
                        k += 1
                        rn +=1

            rn += 1
        if controller.ppvar.get():
            pp_label = tk.Label(self.scrollf, text = 'Run post-simulation analysis: True' , padx = 10, bg='white')
            pp_label.grid(row=rn, column=0, columnspan = 3, sticky = tk.W, pady = (0,10))
            rn += 1
        else:
            pp_label = tk.Label(self.scrollf, text = 'Run post-simulation analysis: False' , padx = 10, bg='white')
            pp_label.grid(row=rn, column=0, columnspan = 3, sticky = tk.W, pady = (0,10))
            rn += 1

        prevb = tk.Button(self.scrollf, text = 'Back',padx=30, command = lambda: self.destroy_canvas(controller))
        prevb.grid(row = rn, column = 3, sticky = tk.E, pady = (30, 0))
        nextb = tk.Button(self.scrollf, text = 'Confirm and run',padx =30, command = lambda: self.save_run_lodis(controller))
        nextb.grid(row = rn, column = 4, sticky = tk.E, pady = (30, 0))


    def destroy_canvas(self, controller):
        self.vsb.destroy()
        self.scrollf.destroy()
        self.canvas.destroy()

        #Reconstruct the canvas
        self.canvas = tk.Canvas(self, background="white", width=850,height=400,scrollregion=(0,0,600, 1100))
        self.scrollf = tk.Canvas(self.canvas, background="white")
        self.vsb = tk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(yscrollcommand=self.vsb.set)

        self.vsb.pack(side="right", fill="y")
        self.canvas.pack(side="left", fill="both", expand=True)
        self.canvas.create_window((0,0), window=self.scrollf, anchor="nw",
                                  tags="self.scrollf")
        #Return to the system params page
        controller.show_frame('Sys_params')

    def save_run_lodis(self, controller):
        #Write into a file all the selected values for the input parameters as a record! <== another function

        controller.run_lodisf90()