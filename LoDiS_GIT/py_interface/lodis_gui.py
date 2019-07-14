import os
import sys
import subprocess
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


import numpy as np
import choose_process as cp
import read_files as rf
import system_params_gui as spg
import run_check as rc
import coalesce_frame as coal
import itMD_frame as itmd
import quench_frame as qf
import growth_frame as gf
import metadyn_frame as mf
import coordnum_frame as cnf
import sfn_frame as sf
import cnbim_frame as cb
import dcom_frame as dcf
import postpro_frame as ppf
import gen_input_in as gin
pd = sys.path[0]                                       #sys.path is a list of directories python will search for imported modules, sys.path[0] is the directory lodis_gui is in.
sys.path.insert(0, os.path.join(pd,r'../py_postpro'))
import lodispp as lpp                                    #Ignore warning/error


"""
***** PYTHON INTERFACE for LoDiS Fortran software *****

Creates a tinker GUI to generate the input parameters for LoDiS to run.
After simulation the option to run post-processing analysis code is available.
Related inputs are coded into individual frames for ease of modification.
To add new process/features create new frames.
NOTE: All frames are generated at the beginning, to modify frames (i.e. use variables with dependencies on other
frames) make use of functions (and the destroy() function to clean frames).
NOTE: To make use of a function from another frame/class you must call the function in it's current situation.
NOTE: input variables not created during the runtime of the GUI MUST BE stated in the class LoDiSApp
preferrably in a dictionary format for ease of modification and grouping.
"""

class LoDiSApp(tk.Tk):

    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        self.title_font = tkfont.Font(family='Helvetica', size=25, weight="bold", slant="italic")
        self.title('LoDiS')
        # the container is where we'll stack a bunch of frames
        # on top of each other, then the one we want visible
        # will be raised above the others
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        # System variables
        #self.atom_rad = {'Pt': 1.385, 'Ag':1.44, 'Pd':1.375, 'Au':1.46}
        self.bimon = tk.BooleanVar()
        self.bimon.set(False)
        self.varmgo = tk.BooleanVar()
        self.varmgo.set(False)
        self.infiles = {'posfile':'N/a', 'potfile':'N/a', 'mgofile':'N/a', 'posfile2':'N/a'}
        self.label_list = ['NVT','itMD', 'Coalescence', 'Quenching', 'Growth',
                      'Metadynamics']
        #self.run_process = 'chosen process'   #Generated in choose_process.py
        #list format instead of dictionary for sys_vn and sys_var due to importance of ordering (for loop used)
        #Additional variables should be added before impl_env index in both lists so that one doesn't have
        #to modify the algorithm
        self.typot = ['rgl', 'lj1', 'gir', 'par']
        self.sys_vn = ['elem1', 'elem2', 'natom',
                       'type_potential','tstep', 'npas',
                       'irand','tinit', 'vnu',
                       'npast','scrivo','sticky_atoms',
                       'sticky_k','cn_cal', 'output_xyz', 'impl_env',
                       'pot_a','pot_b', 'eta_a', 'eta_b']

        self.sys_var = [tk.StringVar(), tk.StringVar(), tk.IntVar(),
                        tk.StringVar(), tk.DoubleVar(),tk.IntVar(),
                        tk.IntVar(),tk.DoubleVar(), tk.DoubleVar(),
                        tk.IntVar(),tk.IntVar(), tk.IntVar(),
                        tk.IntVar(),tk.BooleanVar(),tk.BooleanVar(), tk.BooleanVar(),
                        tk.DoubleVar(),tk.DoubleVar(), tk.DoubleVar(), tk.DoubleVar()]
        tpind = self.sys_vn.index('type_potential')
        self.sys_var[tpind].set('rgl')
        tsind = self.sys_vn.index('tstep')
        self.sys_var[tsind].set(5.0)
        vnind = self.sys_vn.index('vnu')
        self.sys_var[vnind].set(5.0)
        self.envind = self.sys_vn.index('impl_env')
        self.sys_var[self.envind].set(False)
        self.sys_var[self.envind -2].set(False)                  #cn_cal
        self.sys_var[self.envind-1].set(True)                     #output_xyz
        # Note self.sys_out = {key:value for key,value in zip(self.sys_vn, self.sys_var)} is created

        #NVT config
        self.vacf_var = {'vel_af':tk.BooleanVar()}
        self.vacf_var['vel_af'].set(True)
        self.vacf_on = False

        #Coalescence variables
        self.coal = {'natom2': tk.IntVar(),
                     'somedist': tk.DoubleVar()}

        #iTMD
        self.calor = {'deltat':tk.DoubleVar(),            # Temperature step
                      'tcaloric':tk.DoubleVar()}            # Final temperature

        #Quenching
        self.quench = {'itremp': tk.IntVar(),              #Step at which quench starts, itremp < npas unless microcanonical
                       'tmin': tk.DoubleVar()}           #Min temperature to reach

        #Growth
        self.growth = {'ndepmax': tk.IntVar(),                            #number of deposited atoms
                       'lcs': [tk.BooleanVar(), tk.BooleanVar(), tk.BooleanVar()],      #[0]=mono growth, [1]=mixed-shell, [2]=core-shell
                       'at_tipo2': tk.IntVar(),                           #number of species 2 deposited
                       'elemd': tk.StringVar(),                           #Chemical species of deposited atom (species 2)
                       'prob': tk.DoubleVar(),                            #Probability of deposition of species 1
                       'tsorg': tk.DoubleVar(),                           #temperature of the source
                       'rad': tk.DoubleVar()}                             #radius of source
        self.growth['lcs'][0].set(False)
        self.growth['lcs'][1].set(False)
        self.growth['lcs'][2].set(False)
        self.lcs_ind = 1
        # self.lcs_ind = 1,2,3    #Generated in run_check.py

        #Metadyn core variables
        self.cvlist = ['coord_number', 'SFN', 'CN_bim', 'd_com']
        self.cvpage_link = {'coord_number': 'CN_params', 'SFN':'SFN_params',
                            'CN_bim':'CNbim_params', 'd_com':'Dcom_params'}

        self.metacore = {'collvar_wanted': tk.BooleanVar(),         #Write CV values and average forces at each time step for whichever procedure
                         'metaframe': tk.IntVar(),                             #Measuremnt step increments for movie.xyz
                         'gheight': tk.DoubleVar(),                            #Gaussian height (eV)
                         'metaperiod': tk.IntVar(),                            #After how many step a new Gaussian is added
                         'num_cv': 1,                                # Number of CVs used (1 or 2)
                         'collvar_name(1)': tk.StringVar(),
                         'collvar_name(2)': tk.StringVar()}
        self.metacore['collvar_wanted'].set(False)
        self.metacore['collvar_name(1)'].set('coord_number')
        self.metacore['collvar_name(2)'].set('none')

        #Coord number
        self.cn = {'n_pwr': tk.IntVar(),                       #n
                   'm_pwr': tk.IntVar(),                       #m
                   'rzero': tk.DoubleVar(),                    # r0 bulk lattice ref
                   'gwidth': tk.DoubleVar()}                   # Guassian width
        self.cn['n_pwr'].set(6)
        self.cn['m_pwr'].set(12)
        self.cn['rzero'].set(0.147)
        self.cn['gwidth'].set(2)

        #SFN
        self.sfn = {'nsf_pwr': tk.IntVar(),                      #n
                   'msf_pwr': tk.IntVar(),                       #m
                   'dsf': tk.DoubleVar(),                        #d0 bulk lattice ref
                   'rsf': tk.DoubleVar(),                        # r0 bulk lattice ref
                   'gsfwidth': tk.DoubleVar()}                   # Guassian width
        self.sfn['nsf_pwr'].set(6)
        self.sfn['msf_pwr'].set(12)
        self.sfn['dsf'].set(1.354)
        self.sfn['rsf'].set(0.05)
        self.sfn['gsfwidth'].set(25.0)

        #C2N / d3.4N      <= A subset of SFN
        self.c2n_on = tk.BooleanVar()
        self.c2n_on.set(False)
        self.sfn_or_c2n = 'SFN'
        #n2n_pwr = 6,        m2n_pwr = 12
        #d2n = 3.4,          r2n = 0.1
        #g2nwidth = 50.0


        #CN bim
        self.cnbim = {'cn_aa': tk.BooleanVar(),                      # species 1 pairs present?
                      'cn_bb': tk.BooleanVar(),                      # species 2 pairs present?
                      'cn_ab': tk.BooleanVar(),                      # mixed pairs present?
                      'cn_n_pwr': tk.IntVar(),                       # n
                      'cn_m_pwr': tk.IntVar(),                       # m
                      'cn_rzero': tk.DoubleVar(),                    # r0 bulk lattice ref
                      'cn_gwidth': tk.DoubleVar()}                   # Guassian width
        self.cnbim['cn_aa'].set(False), self.cnbim['cn_bb'].set(False)
        self.cnbim['cn_ab'].set(True), self.cnbim['cn_n_pwr'].set(6)
        self.cnbim['cn_m_pwr'].set(12), self.cnbim['cn_rzero'].set(0.147)
        self.cnbim['cn_gwidth'].set(2.0)

        #d COM
        self.dcom = {'d_gwidth': tk.DoubleVar()}                                     # Gaussian width
        self.dcom['d_gwidth'].set(2.0)

        self.cv_dict ={'coord_number':self.cn, 'SFN': self.sfn, 'CN_bim':self.cnbim, 'd_com':self.dcom}

        #Post-processing parameters
        self.postpro = {'T_check': tk.BooleanVar(),
                        'snapshots': tk.StringVar(),
                        'd_min': tk.DoubleVar(),
                        'cutoff': tk.DoubleVar(),
                        'bin_mult': tk.DoubleVar()}
        self.postpro['snapshots'].set('auto')
        self.postpro['T_check'].set(False)
        self.postpro['bin_mult'].set(0.04)
        self.ppvar = tk.BooleanVar()                        #Run postprocessing?
        self.ppvar.set(True)
        #self.snapout <== Snapshots inlist form, generated in postpro_frame

        self.frames = {}
        self.geodict = {}
        #width x height
        frame_dict = {cp.StartPage: '500x350', rf.Input_files: '800x350', spg.Sys_params: '800x500',
                      rc.FinalPage: '850x400', coal.Coalescence_params: '500x250', itmd.IT_params: '600x300',
                      qf.Quench_params: '400x250', gf.Growth_params: '700x250', mf.Metadyn_params: '800x250',
                      cnf.CN_params: '700x200', sf.SFN_params: '700x250', cb.CNbim_params: '700x300',
                      dcf.Dcom_params: '600x150', ppf.Postpro_params: '800x230'}
        #Startpage, Input_files, spg.Sys_params, coal.Coalescence_params
        if sys.version_info[0] < 3:
            for F, geo in frame_dict.iteritems():
                page_name = F.__name__
                frame = F(parent=container, controller=self, *args)
                self.frames[page_name] = frame
                self.geodict[page_name] = geo
                # put all of the pages in the same location;
                # the one on the top of the stacking order
                # will be the one that is visible.
                frame.grid(row=0, column=0, sticky="nsew")
        else:
            for F, geo in frame_dict.items():
                page_name = F.__name__
                frame = F(parent=container, controller=self, *args)
                self.frames[page_name] = frame
                self.geodict[page_name] = geo
                # put all of the pages in the same location;
                # the one on the top of the stacking order
                # will be the one that is visible.
                frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame("StartPage")

    def show_frame(self, page_name):
        #Show a frame for the given page name
        frame = self.frames[page_name]
        geo = self.geodict[page_name]
        self.update_idletasks()
        self.geometry(geo)
        frame.tkraise()

    def run_lodisf90(self):
        #print(mid.moldyn.__doc__)


        # Determine cluster composition and convert into arrays
        elem1 = np.array(list(self.sys_out['elem1'].get()), dtype= 'c')
        elem2 = np.array(list(self.sys_out['elem2'].get()), dtype= 'c')
        if self.bimon.get() in [1, True, '1', '1.0']:
            self.ctype = 'bim'
        else:
            self.ctype = 'mon'


        if self.run_process != 'Metadynamics':
            # Metadynamic param are used in other processes, therefore the default values must be called.
            self.metacore['metaframe'].set(1000)
            self.metacore['metaperiod'].set(1000)
            self.metacore['gheight'].set(0.3)

        if self.run_process == 'Coalescence':            #natom = N + N2
            natval = self.sys_out['natom'].get()
            self.sys_out['natom'].set(natval + self.coal['natom2'].get())

        #Changes CV name to C2N if checkbutton is selected
        if self.c2n_on.get() in [1, 1.0, True]:
            self.sfn_or_c2n = 'C2N'
        else:
            self.sfn_or_c2n = 'SFN'
        if self.metacore['collvar_name(1)'].get() == 'SFN':
            self.metacore['collvar_name(1)'].set(self.sfn_or_c2n)
        elif self.metacore['collvar_name(2)'].get() == 'SFN':
            self.metacore['collvar_name(2)'].set(self.sfn_or_c2n)

        #Generate input.in file
        gin.gen_script(self)

        self.destroy()
        self.run_base_lodis()                         # Run the lodis program via subprocess
        print('main_MD> Process completed')

        #	input args = type_process, tstep, snapID, T_check, d_min, cutoff, bin_mult, natom, ndepmax
        if self.run_process not in ['Coalescence', 'Growth']:
            if self.ppvar.get():
                lpp.md_postpro(self.run_process, self.sys_out['tstep'].get(), self.snapout,
                               self.postpro['T_check'].get(), self.postpro['d_min'].get(), self.postpro['cutoff'].get(),
                               self.postpro['bin_mult'].get(), self.sys_out['natom'].get(), self.growth['ndepmax'].get(),
                               self.atom_rad, self.atom_mass)

    def run_base_lodis(self):
        bashcompath = os.path.join(pd, '../base/LODIS_all')
        print ('bashcommand = %s'%bashcompath)
        in_read = open('input.in', 'r')
        subprocess.check_call(bashcompath, stdin=in_read)
        in_read.close()
        return


if __name__ == "__main__":
    app = LoDiSApp()
    app.mainloop()