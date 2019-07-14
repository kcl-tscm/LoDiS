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
Generates input.in file.
"""
def gen_script(controller):
    f = open('input.in', 'w+')

    # &simul namelist
    f.write('&simul\n')
    f.write('  type_process = \'%s\' \n'%controller.run_process)
    f.write('  %s = \'%s\' \n'%('filepos', controller.infiles['posfile']))
    f.write('  %s = \'%s\' \n'%('filepot', controller.infiles['potfile']))
    if controller.varmgo.get() in [1, 1.0, True]:
        f.write('  mgo_substrate = .True. \n')
        f.write('  %s = \'%s\' \n'%('mgo_pot', controller.infiles['mgofile']))
    else:
        f.write('  mgo_substrate = .False. \n')
        f.write('  mgo_pot = \'None\' \n')

    # System params
    if sys.version_info[0] < 3:
        for key in controller.sys_vn:                      #Use list to maintain order
            if key not in ['type_potential', 'elem1', 'elem2', 'natom']:
                if key not in ['cn_cal', 'output_xyz', 'impl_env']:
                    if isinstance(controller.sys_out[key], tk.DoubleVar):
                        if key == 'tstep':
                            f.write('  tstep = %sd-15 \n'%controller.sys_out[key].get())
                        elif key == 'vnu':
                            f.write('  vnu = %sd11 \n'%controller.sys_out[key].get())
                        else:
                            f.write('  %s = %s \n'%(key, str(controller.sys_out[key].get())+'d0'))

                    elif isinstance(controller.sys_out[key], tk.StringVar):              #writes string with qoutation marks
                        f.write('  %s = \'%s\' \n'%(key, controller.sys_out[key].get()))

                    else:
                        f.write('  %s = %s \n'%(key, controller.sys_out[key].get()))

                else:
                    if controller.sys_out[key].get() in [1, 1.0, True]:
                        bvar = '.True.'
                    else:
                        bvar = '.False.'
                    f.write('  %s = %s \n'%(key, bvar))
    else:
        for key in controller.sys_vn:
            if key not in ['type_potential', 'elem1', 'elem2', 'natom']:
                if key not in ['cn_cal', 'output_xyz', 'impl_env']:
                    if isinstance(controller.sys_out[key], tk.DoubleVar):
                        if key == 'tstep':
                            f.write('  tstep = %sd-15 \n'%controller.sys_out[key].get())
                        elif key == 'vnu':
                            f.write('  vnu = %sd11 \n'%controller.sys_out[key].get())
                        else:
                            f.write('  %s = %s \n'%(key, str(controller.sys_out[key].get())+'d0'))

                    elif isinstance(controller.sys_out[key], tk.StringVar):
                        f.write('  %s = \'%s\' \n'%(key, controller.sys_out[key].get()))
                    else:
                        f.write('  %s = %s \n'%(key, controller.sys_out[key].get()))

                else:
                    if controller.sys_out[key].get() in [1, 1.0, True]:
                        bvar = '.True.'
                    else:
                        bvar = '.False.'
                    f.write('  %s = %s \n'%(key, bvar))


    f.write('/ \n\n')
    f.write('&system\n')
    for vn in ['type_potential', 'elem1', 'elem2']:            #strings
        f.write('  %s = \'%s\' \n'%(vn, controller.sys_out[vn].get()))
    f.write('  natom = %s \n'%controller.sys_out['natom'].get())
    f.write('  fattor = 1.d0 \n')
    f.write('  %s = \'%s\' \n'%('sys', controller.ctype))
    f.write('/ \n\n')



    if controller.run_process not in ['Metadynamics', 'Growth']:
        if controller.run_process == 'NVT':
            pro_vars = controller.vacf_var
            f.write('&canon\n')
        elif controller.run_process == 'itMD':
            pro_vars = controller.calor
            f.write('&calor\n')
        elif controller.run_process == 'Coalescence':
            pro_vars = controller.coal
            f.write('&coal\n')
            f.write('  filepos2 = \'%s\' \n'%controller.infiles['posfile2'])
        elif controller.run_process == 'Quenching':
            pro_vars = controller.quench
            f.write('&quench\n')
        else:
            print ('Error unknown process detected')
            sys.exit()
        if sys.version_info[0] < 3:                  #python 2
            for key, value in pro_vars.iteritems():
                if isinstance(value, tk.DoubleVar):
                    f.write('  %s = %s \n'%(key, str(value.get())+'d0'))
                elif isinstance(value, tk.StringVar):
                    f.write('  %s = \'%s\' \n'%(key, value.get()))
                else:
                    f.write('  %s = %s \n'%(key, value.get()))

        else:                                              #python 3
            for key, value in pro_vars.items():
                if isinstance(value, tk.DoubleVar):
                    f.write('  %s = %s \n'%(key, str(value.get())+'d0'))
                elif isinstance(value, tk.StringVar):
                    f.write('  %s = \'%s\' \n'%(key, value.get()))
                else:
                    f.write('  %s = %s \n'%(key, value.get()))
        f.write('/ \n\n')

    elif controller.run_process == 'Growth':
        f.write('&growth\n')
        pro_vars = controller.growth
        lcs_ind = [l.get() for l in controller.growth['lcs']].index(True) + 1     #lcs = 1, 2,3
        if sys.version_info[0] < 3:
            for key, value in pro_vars.iteritems():
                if key != 'lcs':
                    if isinstance(value, tk.DoubleVar):
                        f.write('  %s = %s \n'%(key, str(value.get())+'d0'))
                    elif isinstance(value, tk.StringVar):
                        f.write('  %s = \'%s\' \n'%(key, value.get()))
                    else:
                        f.write('  %s = %s \n'%(key, value.get()))

                else:
                    f.write('  lcs = %s \n'%lcs_ind)

        else:                    #python 3
            for key, value in pro_vars.items():
                if key != 'lcs':
                    if isinstance(value, tk.DoubleVar):
                        f.write('  %s = %s \n'%(key, str(value.get())+'d0'))
                    elif isinstance(value, tk.StringVar):
                        f.write('  %s = \'%s\' \n'%(key, value.get()))
                    else:
                        f.write('  %s = %s \n'%(key, value.get()))
                else:

                     f.write('  lcs = %s \n'%lcs_ind)
        f.write('/ \n\n')


    #Metdynamics
    pro_vars = controller.metacore
    f.write('&metalist\n')
    if sys.version_info[0] < 3:              #python 2
        for key, value in pro_vars.iteritems():
            if key not in ['collvar_wanted', 'num_cv']:
                if isinstance(value, tk.DoubleVar):
                    f.write('  %s = %s \n'%(key, str(value.get())+'d0'))
                elif isinstance(value, tk.StringVar):
                    f.write('  %s = \'%s\' \n'%(key, value.get()))
                else:
                    f.write('  %s = %s \n'%(key, value.get()))

            elif key == 'collvar_wanted':
                if value.get() in [1, 1.0, True]:
                    cw = '.True.'
                else:
                    cw = '.False.'
                f.write('  collvar_wanted = %s \n'%cw)
            else:
                f.write('  num_cv = %s \n'%pro_vars['num_cv'])
        f.write('\n')


        # CV params
        f.write('  !===== %s ===== \n'%(pro_vars['collvar_name(1)'].get()))
        if pro_vars['collvar_name(1)'].get() == 'coord_number':
            cv1_vars = controller.cn
        elif pro_vars['collvar_name(1)'].get() in ['SFN', 'C2N']:
            cv1_vars = controller.sfn
        elif pro_vars['collvar_name(1)'].get() == 'CN_bim':
            cv1_vars = controller.cnbim
        else:
            cv1_vars = controller.dcom

        for key, value in cv1_vars.iteritems():
            if isinstance(value, tk.DoubleVar):
                f.write('  %s = %s \n'%(key, str(value.get())+'d0'))
            else:
                f.write('  %s = %s \n'%(key, value.get()))

        f.write('\n')

        if pro_vars['collvar_name(2)'].get() != 'none':
            f.write('  !===== %s ===== \n'%(pro_vars['collvar_name(2)'].get()))

            if pro_vars['collvar_name(2)'].get() == 'coord_number':
                cv2_vars = controller.cn
            elif pro_vars['collvar_name(2)'].get() in ['SFN', 'C2N']:
                cv2_vars = controller.sfn
            elif pro_vars['collvar_name(2)'].get() == 'CN_bim':
                cv2_vars = controller.cnbim
            else:
                cv2_vars = controller.dcom

            for key, value in cv2_vars.iteritems():
                if isinstance(value, tk.DoubleVar):
                    f.write('  %s = %s \n'%(key, str(value.get())+'d0'))
                else:
                    f.write('  %s = %s \n'%(key, value.get()))


    else:                                     #python 3
        for key, value in pro_vars.items():
            if key not in ['collvar_wanted', 'num_cv']:
                if isinstance(value, tk.DoubleVar):
                    f.write('  %s = %s \n'%(key, str(value.get())+'d0'))
                elif isinstance(value, tk.StringVar):
                    f.write('  %s = \'%s\' \n'%(key, value.get()))
                else:
                    f.write('  %s = %s \n'%(key, value.get()))

            elif key == 'collvar_wanted':
                if value.get() in [1, 1.0, True]:
                    cw = '.True.'
                else:
                    cw = '.False.'
                f.write('  collvar_wanted: %s \n'%cw)
            else:
                f.write('  num_cv = %s \n'%pro_vars['num_cv'])
        f.write('\n')

        # CV params
        v = 0
        f.write('  !===== %s ====='%(pro_vars['collvar_name(1)'].get()))

        if pro_vars['collvar_name(1)'].get() == 'coord_number':
            cv1_vars = controller.cn
        elif pro_vars['collvar_name(1)'].get() in ['SFN', 'C2N']:
            cv1_vars = controller.sfn
        elif pro_vars['collvar_name(1)'].get() == 'CN_bim':
            cv1_vars = controller.cnbim
        else:
            cv1_vars = controller.dcom

        for key, value in cv1_vars.items():
            if isinstance(value, tk.DoubleVar):
                f.write('  %s = %s \n'%(key, str(value.get())+'d0'))
            else:
                f.write('  %s = %s \n'%(key, value.get()))
        f.write('\n')

        if pro_vars['collvar_name(2)'].get() != 'none':
            f.write('  !===== %s ====='%pro_vars['collvar_name(2)'].get())

            if pro_vars['collvar_name(2)'].get() == 'coord_number':
                cv2_vars = controller.cn
            elif pro_vars['collvar_name(2)'].get() in ['SFN', 'C2N']:
                cv2_vars = controller.sfn
            elif pro_vars['collvar_name(2)'].get() == 'CN_bim':
                cv2_vars = controller.cnbim
            else:
                cv2_vars = controller.dcom

            for key, value in cv2_vars.items():
                if isinstance(value, tk.DoubleVar):
                    f.write('  %s = %s \n'%(key, str(value.get())+'d0'))
                else:
                    f.write('  %s = %s \n'%(key, value.get()))
    f.write('/ \n')

    f.close()
    return