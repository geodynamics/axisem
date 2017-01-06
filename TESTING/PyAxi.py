#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  PyAxi.py
#   Purpose:   Python interface for AXISEM
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#-------------------------------------------------------------------

#for debugging: import ipdb; ipdb.set_trace()

#-----------------------------------------------------------------------
#----------------Import required Modules (Python)-----------------------
#-----------------------------------------------------------------------

# Added this line for python 2.5 compatibility
from __future__ import with_statement
import matplotlib as mpl
mpl.rcParams['backend'] = 'agg'
import os
import sys
import time
import subprocess
import warnings
import ConfigParser
import shutil
import glob

print '\n   -----------------------------------------'
bold = "\033[1m"
reset = "\033[0;0m"
print '\t\t' + bold + '   PyAxi' + reset
print '\t' + '   Python interface for Axisem'
#print '\n   Developed by Kasra Hosseini'
#print '   email: hosseini@geophysik.uni-muenchen.de'
print '   -----------------------------------------\n'


global obspy_error
obspy_error = 'N'

try:
    import numpy as np
    import matplotlib.pyplot as plt
    from obspy.core import read, Trace, Stream, UTCDateTime
    from obspy.signal import util
    from obspy.signal.invsim import cosine_taper
except Exception, error:
    obspy_error = 'obspy (http://obspy.org/) should be installed'
    obspy_error += 'for using MSEED, Filter, STF and TEST functionalities.'
    print "\n======================================"
    print 'Error in importing:'
    print error
    print "------------------"
    print "No MSEED, Filter, STF and TEST functionalities (refer to inpython.cfg)"
    print obspy_error
    print "======================================"
    print "\n\n-----------------------------------"
    print 'AXISEM will start running in 10 sec'
    print "-----------------------------------\n"
    time.sleep(10)

########################################################################
############################# Main Program #############################
########################################################################

def PyAxi(**kwargs):

    """
    PyAxi is the function dedicated to the main part of the code.

    To run this code:
    1. change the "inpython.cfg" file based on the inputs that you want
    2. type: python PyAxi.py
    3. if you have run the code once before with the same name, then PyAxi
       will ask whether it should remove the old directory or not!
    """

    # global variables
    global input
    try:
        if sys.argv[1] == '--check':
            import copy
            print 'Check the Basic, Processing and Visualization requirements:\n'
            p = []
            print repr('Compiler').rjust(20).replace("'", '') + "    |    " +   repr('Installed').rjust(11).replace("'", '')
            print '  ' + 42*'-'

            for compilers in ['gfortran', 'ifort', 'mpif90.openmpi', 'mpif90',
                                'gnuplot', 'taup', 'paraview', 'matlab', 'googleearth']:
                if manual_which(compilers) == None: p.append([compilers, 'N'])
                else: p.append([compilers, 'Y'])

            p_print = copy.deepcopy(p)

            if p_print[3][1] == 'N' and p_print[2][1] != 'N': del p_print[3]
            elif p_print[3][1] != 'N' and p_print[2][1] == 'N': del p_print[2]

            if p_print[1][1] == 'N' and p_print[0][1] != 'N': del p_print[1]
            elif p_print[1][1] != 'N' and p_print[0][1] == 'N': del p_print[0]

            for comp in range(0, len(p_print)):
                print repr(p_print[comp][0]).rjust(20).replace("'", '') + "    |    " +  repr(p_print[comp][1]).rjust(7).replace("'", '')

            print '\nSummary:\n'
            if p[0][1] != 'N' or p[1][1] != 'N':
                if p[2][1] != 'N' or p[3][1] != 'N':
                    print 'Basic functionality requirement...CHECK'
                else:
                    print 'Basic functionality requirement...ERROR'
            else: print 'Basic functionality requirement...ERROR'
            if p[4][1] != 'N' and p[5][1] != 'N':
                print 'Processing requirement........CHECK'
            else: print 'Processing requirement........ERROR'
            if p[6][1] != 'N' and p[7][1] != 'N' and p[8][1] != 'N':
                print 'Visualization tools...CHECK'
            else: print 'Visualization tools...ERROR'
            print '\n'
            sys.exit(2)
    except Exception, error:
        print '',
        sys.exit(2)

    # ------------------Read INPUT file (Parameters)--------------------
    read_input_file()

    if input['test'] != 'N':
        subprocess.check_call(['cp',
            os.path.join(input['test_folder'], 'STATIONS'),
            os.path.join(input['axi_address'], 'SOLVER', 'STATIONS')])
        read_input_file()

    # Defining error message
    output = 0
    output_print = 20*'*' + '\n' + str(output) + '\n' + 20*'*'


    if os.path.exists(os.path.join(input['axi_address'], 'SOLVER', input['solver_name'])) == True:
        print '--------------------------------------------------------'
        user_raw_input = raw_input('Directory with the same name as ' + input['solver_name'] +
                ' \n(SOLVER_NAME in ipython.cfg) exists in:\n' +
                os.path.join(input['axi_address'], 'SOLVER', input['solver_name']) + '\n\n' +
                'You could:\n' +
                'S: Stop the program\n' +
                'R: Remove the directory and continue the program\n' +
                'C: Continue the program without removing the directory\n\n').upper()

        print '--------------------------------------------------------'
        if user_raw_input == 'C':
            print "Continue the program!"
            print '--------------------------------------------------------'
        elif user_raw_input == 'R':
            print "Removing the directory and continue the program!"
            print '--------------------------------------------------------'
            subprocess.check_call(['rm', '-rf',
                os.path.join(input['axi_address'], 'SOLVER', input['solver_name'])])
        elif user_raw_input == 'S':
            print 'EXIT AXISEM!'
            sys.exit()

    ##############################################################
    ######################COPY TEMPLATES #########################
    ##############################################################
    stdout_param = None
    os.chdir(input['axi_address'])
    output = subprocess.check_call(['./copytemplates.sh'], stdout = stdout_param)
    #output = subprocess.check_call(['cp', os.path.join('.', 'MESHER', 'Makefile.TEMPLATE'), \
    #                                    os.path.join('.', 'MESHER', 'Makefile')], stdout = stdout_param)
    #if output != 0: print output_print
    #output = subprocess.check_call(['cp', os.path.join('.', 'SOLVER', 'Makefile.TEMPLATE'), \
    #                                    os.path.join('.', 'SOLVER', 'Makefile')], stdout = stdout_param)
    #if output != 0: print output_print

    ##############################################################
    ################CREATE make_axisem.macros ####################
    ##############################################################
    print '\nCreating the make_axisem.macros...',
    os.chdir(input['axi_address'])
    make_axisem_fio = open(os.path.join('./make_axisem.macros'), 'w')
    make_axisem_ls = []
    path_mpirun = os.path.expandvars(input['make_axisem_mpirun'])
    make_axisem_ls.append('MPIRUN = %s \n' % path_mpirun)
    make_axisem_ls.append('USE_NETCDF = %s \n' %(input['make_axisem_use_netcdf']))
    path_netcdf = os.path.expandvars(input['make_axisem_netcdf_path'])
    make_axisem_ls.append('NETCDF_PATH = %s \n' % path_netcdf)
    make_axisem_ls.append('SERIAL = %s \n' %(input['make_axisem_serial']))
    make_axisem_ls.append('INCLUDE_MPI = %s \n' %(input['make_axisem_include_mpi']))
    path_cc = os.path.expandvars(input['make_axisem_CC'])
    make_axisem_ls.append('CC = %s \n' % path_cc)
    path_fc = os.path.expandvars(input['make_axisem_FC'])
    make_axisem_ls.append('FC = %s \n' % path_fc)
    make_axisem_ls.append('FFLAGS = %s \n' %(input['make_axisem_FFLAGS']))
    make_axisem_ls.append('CFLAGS = %s \n' %(input['make_axisem_CFLAGS']))
    make_axisem_ls.append('LDFLAGS = %s \n' %(input['make_axisem_LDFLAGS']))
    for make_axi_item in make_axisem_ls:
        if input['verbose'] == 'Y':
            print(make_axi_item)
        make_axisem_fio.writelines(make_axi_item)
    make_axisem_fio.close()
    print 'DONE'

    ##############################################################
    ############################# MESHER #########################
    ##############################################################
    t1_mesher = time.time()

    if input['mesher'] != 'N':
        print "\n======"
        print "MESHER"
        print "======"
        mesher_dir = os.path.join(input['axi_address'], 'MESHER')
        os.chdir(mesher_dir)
        print os.getcwd()

        # Delete previous mesh_params.h, meshdb.dat*, Diags/*
        if input['verbose'] != 'N':
            print "\n====================================================="
            print 'Removing old mesh_params.h*, meshdb.dat* and ./Diags/*'
            print "=====================================================\n"
            stdout_param = None
        else:
            sys.stdout.write('Remove old files...')
            sys.stdout.flush()
            stdout_param = subprocess.PIPE
        output = subprocess.check_call(['rm', '-rf', './mesh_params.h*'], stdout = stdout_param)
        if output != 0: print output_print
        output = subprocess.check_call(['rm', '-rf', './meshdb.dat*'], stdout = stdout_param)
        if output != 0: print output_print
        output = subprocess.check_call(['rm', '-rf', './Diags/*'], stdout = stdout_param)
        if output != 0: print output_print
        #output = subprocess.check_call(['rm', '-rf', './unrolled_loops.f90'], stdout = stdout_param)
        #if output != 0: print output_print
        print 'DONE'

        # Change the input files + make clean; make
        if input['mesher_make'] != 'N':
            if input['verbose'] != 'N':
                print "==================="
                print "Change inparam_mesh"
                print "==================="
            else:
                sys.stdout.write('Change inparam_mesh...')
                sys.stdout.flush()
            if os.path.isfile('./inparam_mesh'):
                subprocess.check_call(['rm', './inparam_mesh'])
            #subprocess.check_call(['cp', './inparam_mesh.TEMPLATE', './inparam_mesh'])
            inparam_mesh_input = []
            inparam_mesh_input.append('BACKGROUND_MODEL     %s\n' %(input['mesher_bg_model']))
            inparam_mesh_input.append('EXT_MODEL     %s\n' %(input['mesher_ext_model']))
            inparam_mesh_input.append('DOMINANT_PERIOD     %s\n' %(input['mesher_dominant_period']))
            inparam_mesh_input.append('NTHETA_SLICES     %s\n' %(input['mesher_ntheta']))
            inparam_mesh_input.append('NRADIAL_SLICES     %s\n' %(input['mesher_nradial']))
            inparam_mesh_input.append('WRITE_VTK     %s\n' %(input['mesher_write_vtk']))
            inparam_mesh_input.append('COARSENING_LAYERS     %s\n' %(input['mesher_coarsening_layers']))
            inparam_mesh_input.append('IC_SHEAR_WAVE     %s\n' %(input['mesher_ic_shear_wave']))
            inparam_mesh_input.append('NPOL     %s\n' %(input['mesher_npol']))
            inparam_mesh_input.append('EL_PER_LAMBDA     %s\n' %(input['mesher_el_per_lambda']))
            inparam_mesh_input.append('COURANT_NR     %s\n' %(input['mesher_courant_nr']))
            inparam_mesh_input.append('RADIUS     %s\n' %(input['mesher_radius']))
            inparam_mesh_input.append('SAVE_MESH     %s\n' %(input['mesher_save_mesh']))
            inparam_mesh_input.append('VERBOSE     %s\n' %(input['mesher_verbose']))
            inparam_mesh_input.append('NCPU      4')

            inparam_mesh_open = open('./inparam_mesh', 'w')
            for i in range(0, len(inparam_mesh_input)):
                inparam_mesh_open.write(inparam_mesh_input[i])
            inparam_mesh_open.close()

            if input['verbose'] != 'N':
                for inp_item in inparam_mesh_input:
                    print inp_item,
            else:
                print 'DONE'

            if input['verbose'] != 'N':
                print "\n==============="
                print "RUN submit.csh"
                print "================="
                stdout_param = None
            else:
                sys.stdout.write('submit.csh...')
                sys.stdout.flush()
                stdout_param = subprocess.PIPE

            output = subprocess.check_call(os.path.join(mesher_dir,
                                                        'submit.csh'),
                                           stdout=stdout_param)
            if output != 0: print output_print

            # check the xmesh whether it is finished or not!
            process_name= "xmesh"
            check_process="running"
            t_check = 0
            time.sleep(2)
            print '--------------------------'
            print 'submit.csh is running for:'
            while check_process=="running":
                tmp = os.popen("ps -Af").read()
                if process_name in tmp[:]:
                    print str(t_check) + ', ',
                    time.sleep(1)
                    t_check+=2
                else:
                    output_file_open =  open('OUTPUT', 'r')
                    output_file_read = output_file_open.readlines()
                    test = output_file_read[-1].find('....DONE WITH MESHER !')
                    if test==-1:
                        print "\n============================="
                        print '\n       MESHER crashed'
                        print "\n============================="
                        raise RuntimeError('Mesher crashed')
                        return
                    check_process='DONE'
            print 'DONE'
            print '--------------------------'

        # move the mesh to the SOLVER folder
        if input['mesher_move'] != 'N':
            if input['verbose'] != 'N':
                print "==================="
                print "Move mesh files to:"
                print os.path.join('..', 'SOLVER', 'MESHES', input['mesh_name'])
                print "==================="
                stdout_param = None
            else:
                stdout_param = subprocess.PIPE

            if os.path.isdir(os.path.join('..', 'SOLVER', 'MESHES', input['mesh_name'])):
                if input['verbose'] != 'N':
                    print '\n#######################################################################'
                    print "Remove the mesh files from MESHES dir (before running the movemesh.csh)"
                    subprocess.check_call(['rm', '-rf',
                        os.path.join('..', 'SOLVER', 'MESHES', input['mesh_name'])])
                    print os.path.join('..', 'SOLVER', 'MESHES', input['mesh_name'])
                    print "is removed"
                    print '#######################################################################\n'
                else:
                    sys.stdout.write('Remove the old mesh files...')
                    sys.stdout.flush()
                    subprocess.check_call(['rm', '-rf',
                        os.path.join('..', 'SOLVER', 'MESHES', input['mesh_name'])])
                    print 'DONE'

            sys.stdout.write('Move the mesh files...')
            sys.stdout.flush()
            output = subprocess.check_call(['./movemesh.csh', input['mesh_name']], stdout = stdout_param)
            if output != 0: print output_print
            print 'DONE'

    t2_mesher = time.time()

    ##############################################################
    ############################# SOLVER #########################
    ##############################################################
    t1_solver = time.time()

    if input['solver'] != 'N':
        print "\n======"
        print "SOLVER"
        print "======"
        solver_dir = os.path.join(input['axi_address'], 'SOLVER')
        os.chdir(solver_dir)

        # Change the input files + make clean; make
        if input['solver_make'] != 'N':
            if input['verbose'] != 'N':
                print "\n===================="
                print "Change inparam_basic"
                print "===================="
            else:
                sys.stdout.write('Change inparam_basic...')
                sys.stdout.flush()

            if os.path.isfile('inparam_basic'):
                subprocess.check_call(['rm', 'inparam_basic'])

            inparam_basic_input = []
            inparam_basic_input.append('SIMULATION_TYPE     %s\n' %(input['solver_sim_type']))
            inparam_basic_input.append('SEISMOGRAM_LENGTH     %s\n' %(input['solver_seis_length']))
            inparam_basic_input.append('RECFILE_TYPE     %s\n' %(input['solver_recfile_type']))
            inparam_basic_input.append('MESHNAME     %s\n' %(input['mesh_name']))
            inparam_basic_input.append('LAT_HETEROGENEITY     %s\n' %(input['solver_lat_heterogeneity']))
            inparam_basic_input.append('ATTENUATION     %s\n' %(input['solver_attenuation']))
            inparam_basic_input.append('SAVE_SNAPSHOTS     %s\n' %(input['solver_save_snapshots']))
            inparam_basic_input.append('VERBOSITY     %s\n' %(input['solver_verbosity']))

            inparam_solver_open = open('./inparam_basic', 'w')
            for i in range(0, len(inparam_basic_input)):
                inparam_solver_open.write(inparam_basic_input[i])
            inparam_solver_open.close()
            if input['verbose'] != 'N':
                for basic_line in inparam_basic_input: print basic_line,
            else:
                print 'DONE'

            if input['verbose'] != 'N':
                print "\n======================="
                print "Change inparam_advanced"
                print "======================="
            else:
                sys.stdout.write('Change inparam_advanced...')
                sys.stdout.flush()

            if os.path.isfile('inparam_advanced'):
                subprocess.check_call(['rm', 'inparam_advanced'])
            inparam_advanced_input = []
            inparam_advanced_input.append('SAMPLING_PERIOD         %s\n' %(input['solver_sampling_rate']))
            inparam_advanced_input.append('TIME_STEP               %s\n' %(input['solver_time_step']))
            inparam_advanced_input.append('SOURCE_PERIOD           %s\n' %(input['solver_source_period']))
            inparam_advanced_input.append('SOURCE_FUNCTION         %s\n' %(input['source_stf']))
            inparam_advanced_input.append('TIME_SCHEME             %s\n' %(input['solver_time_scheme']))
            inparam_advanced_input.append('DATA_DIR               "%s"\n' %(input['solver_data_dir']))
            inparam_advanced_input.append('INFO_DIR               "%s"\n' %(input['solver_info_dir']))
            inparam_advanced_input.append('DIAGNOSTIC_FILE_OUTPUT "%s"\n' %(input['solver_diag_file_output']))
            inparam_advanced_input.append('MESH_TEST               %s\n' %(input['solver_mesh_test']))
            inparam_advanced_input.append('DEFLATE_LEVEL           %s\n' %(input['solver_deflate_level']))
            inparam_advanced_input.append('SNAPSHOT_DT             %s\n' %(input['solver_snapshot_dt']))
            inparam_advanced_input.append('SNAPSHOTS_FORMAT        %s\n' %(input['solver_snapshots_format']))
            inparam_advanced_input.append('USE_NETCDF              %s\n' %(input['make_axisem_use_netcdf']))
            inparam_advanced_input.append('KERNEL_WAVEFIELDS       %s\n' %(input['solver_kernel_wavefields']))
            inparam_advanced_input.append('KERNEL_SPP              %s\n' %(input['solver_kernel_spp']))
            inparam_advanced_input.append('KERNEL_SOURCE           %s\n' %(input['solver_kernel_source']))
            inparam_advanced_input.append('KERNEL_IBEG             %s\n' %(input['solver_kernel_ibeg']))
            inparam_advanced_input.append('KERNEL_IEND             %s\n' %(input['solver_kernel_iend']))
            inparam_advanced_input.append('NR_LIN_SOLIDS           %s\n' %(input['solver_nr_lin_solids']))
            inparam_advanced_input.append('F_MIN                   %s\n' %(input['solver_fmin']))
            inparam_advanced_input.append('F_MAX                   %s\n' %(input['solver_fmax']))
            inparam_advanced_input.append('F_REFERENCE             %s\n' %(input['solver_fref']))
            inparam_advanced_input.append('SMALL_Q_CORRECTION      %s\n' %(input['solver_small_q_correction']))
            inparam_advanced_input.append('NR_F_SAMPLE             %s\n' %(input['solver_nr_f_sample']))
            inparam_advanced_input.append('MAXINT_SA               %s\n' %(input['solver_maxint_sa']))
            inparam_advanced_input.append('TSTART_SR               %s\n' %(input['solver_tstart_sr']))
            inparam_advanced_input.append('TSTART_AMP              %s\n' %(input['solver_tstart_amp']))
            inparam_advanced_input.append('T_DECAY                 %s\n' %(input['solver_t_decay']))
            inparam_advanced_input.append('FIX_FREQ                %s\n' %(input['solver_fix_freq']))
            inparam_advanced_input.append('DUMP_VTK                %s\n' %(input['solver_dump_vtk']))
            inparam_advanced_input.append('COARSE_GRAINED          %s\n' %(input['solver_coarse_grained']))
            inparam_advanced_input.append('SAVE_ENERGY             %s\n' %(input['solver_save_energy']))
            inparam_advanced_input.append('HOMO_MODEL              %s\n' %(input['solver_homo_model']))
            inparam_advanced_input.append('HOMO_VP                 %s\n' %(input['solver_homo_vp']))
            inparam_advanced_input.append('HOMO_VS                 %s\n' %(input['solver_homo_vs']))
            inparam_advanced_input.append('HOMO_RHO                %s\n' %(input['solver_homo_rho']))
            inparam_advanced_input.append('FORCE_ANISO             %s\n' %(input['solver_force_aniso']))

            inparam_solver_open = open('./inparam_advanced', 'w')
            for i in range(0, len(inparam_advanced_input)):
                inparam_solver_open.write(inparam_advanced_input[i])
            inparam_solver_open.close()
            if input['verbose'] != 'N':
                for advanced_line in inparam_advanced_input: print advanced_line,
            else:
                print 'DONE'

            print 'source_type: ' + input['source_type'];
            if input['sourcefile_type'] == 'sourceparams':
                if input['verbose'] != 'N':
                    print "\n==========================================="
                    print "Change the Source params (inparams_source)"
                    print "==========================================="
                else:
                    sys.stdout.write('Change the Source params (inparams_source)...')
                    sys.stdout.flush()

                source_text = [];
                source_text.append("SOURCE_TYPE  "     + input['source_type'] + "\n")
                source_text.append("SOURCE_DEPTH "     + input['source_depth'] + "\n")
                source_text.append("SOURCE_LAT   "     + input['source_lat'] + "\n")
                source_text.append("SOURCE_LON   "     + input['source_lon'] + "\n")
                source_text.append("SOURCE_AMPLITUDE " + input['source_amp'] + "\n")

                if os.path.isfile('inparam_source'):
                    subprocess.check_call(['rm', 'inparam_source'])
                source_open = open('./inparam_source', 'w')

                for i in range(0, len(source_text)):
                    source_open.write(source_text[i])

                source_open.close()
                if input['verbose'] != 'N':
                    for i in range(0, len(source_text)):
                        print source_text[i];
                else:
                    print 'DONE'

            elif input['sourcefile_type'] == 'cmtsolut':
                if input['verbose'] != 'N':
                    print "\n======================================"
                    print "Change the Source params (CMTSOLUTION)"
                    print "======================================"
                else:
                    sys.stdout.write('Change the Source params (CMTSOLUTION)...')
                    sys.stdout.flush()

                if os.path.isfile('CMTSOLUTION'):
                    subprocess.check_call(['rm', 'CMTSOLUTION'])
                subprocess.check_call(['cp', 'CMTSOLUTION.TEMPLATE', 'CMTSOLUTION'])

                source_open = open('./CMTSOLUTION', 'r')
                source_read = source_open.readlines()

                source_read[0] = "PyAxi generated" + '\n'
                source_read[4] = 'latitude:      ' + input['cmt_lat'] + '\n'
                source_read[5] = 'longitude:     ' + input['cmt_lon'] + '\n'
                source_read[6] = 'depth:         ' + input['cmt_depth'] + '\n'
                source_read[7] = 'Mrr:      ' + input['cmt_Mrr'] + '\n'
                source_read[8] = 'Mtt:       ' + input['cmt_Mtt'] + '\n'
                source_read[9] = 'Mpp:      ' + input['cmt_Mpp'] + '\n'
                source_read[10] = 'Mrt:      ' + input['cmt_Mrt'] + '\n'
                source_read[11] = 'Mrp:       ' + input['cmt_Mrp'] + '\n'
                source_read[12] = 'Mtp:      ' + input['cmt_Mtp'] + '\n'

                source_open.close()
                source_open = open('./CMTSOLUTION', 'w')

                for i in range(0, len(source_read)):
                    source_open.write(source_read[i])

                source_open.close()
                if input['verbose'] != 'N':
                    print source_read[0] + source_read[4] + source_read[5] + \
                                source_read[6] + source_read[7] + \
                                source_read[8] + source_read[9] + \
                                source_read[10] + source_read[11] + \
                                source_read[12]
                else:
                    print 'DONE'

            if input['receiver_type'] == 'colatlon':
                if input['verbose'] != 'N':
                    print "\n=========================================="
                    print "Change the Receiver params (receivers.dat)"
                    print "=========================================="
                else:
                    sys.stdout.write('Change the Receiver params (receivers.dat)...')
                    sys.stdout.flush()

                if not os.path.isfile('receivers.dat'):
                    print "\n'receivers.dat' file does not exist!"
                else:
                    receiver_open = open('./receivers.dat', 'r')
                    receiver_read = receiver_open.readlines()
                    print "Number of receivers: " + receiver_read[0] + \
                            "Number of lines    : " + str(len(receiver_read)-1)
                    if int(receiver_read[0]) != (len(receiver_read) - 1):
                        "Number of stations is not entered correctly!"

            elif input['receiver_type'] == 'stations':
                if input['verbose'] != 'N':
                    print "\n====================================="
                    print "Change the Receiver params (STATIONS)"
                    print "====================================="
                else:
                    sys.stdout.write('Change the Receiver params (STATIONS)...')
                    sys.stdout.flush()

                if input['IO_STATION']:
                    subprocess.check_call(['cp',
                        os.path.join(input['IO_STATION']),
                        os.path.join(input['axi_address'], 'SOLVER', 'STATIONS')])
                if not os.path.isfile('STATIONS'):
                    subprocess.check_call(['cp', 'STATIONS.TEMPLATE', 'STATIONS'])
                receiver_open = open('./STATIONS', 'r')
                receiver_read = receiver_open.readlines()
                print "Number of receivers: " + str(len(receiver_read))

            if input['verbose'] != 'N':
                print "\n=============="
                print "RUN submit.csh"
                print "=============="
                stdout_param = None
            else:
                sys.stdout.write('submit.csh...')
                sys.stdout.flush()
                stdout_param = subprocess.PIPE

            output = subprocess.check_call(['./submit.csh', input['solver_name']], stdout = stdout_param)
            if output != 0: print output_print
            print 'DONE'

            # copy the input file to the folder (for reference!)
            if input['verbose'] != 'N':
                print "=================="
                print "Copy the input in:\n" + os.path.join(input['axi_address'], 'SOLVER', input['solver_name'])
                print "=================="
            else:
                sys.stdout.write('Copy the input...')
                sys.stdout.flush()

            if not os.path.exists(os.path.join(input['axi_address'], 'SOLVER', input['solver_name'])):
                print 'WARNING: $s is not created! Wait for 10 sec and re-check.' \
                        %(os.path.join(input['axi_address'], 'SOLVER', input['solver_name']))
                time.sleep(10)
                if not os.path.exists(os.path.join(input['axi_address'], 'SOLVER', input['solver_name'])):
                    sys.exit('...ERROR...%s is not created.'
                            %(os.path.join(input['axi_address'], 'SOLVER', input['solver_name'])))
            subprocess.check_call(['cp', input['inpython_address'],
                        os.path.join(input['axi_address'], 'SOLVER', input['solver_name'], "inpython.cfg")])
            print 'DONE'

            # check the OUTPUT file and inform the user whenever the program is done!
            print "\n================"
            print "Check the OUTPUT"
            print "================"

            os.chdir(os.path.join(input['axi_address'], 'SOLVER', input['solver_name']))

            test = -1; test_1 = -1; test_2 = -1; test_3 = -1; test_4 = -1
            solverended = 0
            solvercrashed = 0
            same_output = 0
            time.sleep(10)
            #print_output = "Just after 2 seconds!"

            if input['sourcefile_type'] == 'sourceparams':
                print_output = ""
                if not os.path.exists('OUTPUT_' + input['solver_name']):
                    print 'WARNING: %s is not created! wait for 10 seconds and re-check.' \
                                %('OUTPUT_' + input['solver_name'])
                    time.sleep(10)
                    if not os.path.exists('OUTPUT_' + input['solver_name']):
                        sys.exit('...ERROR...%s is not created.' %('OUTPUT_' + input['solver_name']))
                while (test == -1 and solverended == 0):
                    if howmanyofthisprocess('axisem') == 0:
                        solverended = 1

                    output_file_open =  open('OUTPUT_' + input['solver_name'], 'r')
                    output_file_read = output_file_open.readlines()
                    test = output_file_read[-1].find('PROGRAM axisem FINISHED')
                    if output_file_read[-1].find('PROGRAM axisem FINISHED') == -1:
                        test = -1
                        last_output = print_output
                        print_output = output_file_read[-1].split('\n')[0]
                        # Check whether Output is stagnant since more than one minute.
                        if (last_output == print_output) :
                            same_output = same_output + 1
                        else :
                            same_output = 0
                        if (same_output >= 30) :
                            print "\n============================="
                            print '\n       SOLVER froze'
                            print '\n      killing axisem'
                            print '\n        I love it!'
                            print "\n============================="
                            temp = killprocess('axisem')
                            return

                    if input['verbose'] == 'Y':
                        print print_output
                    time.sleep(2)

                if  test==-1:
                    solvercrashed = 1

            elif input['sourcefile_type'] == 'cmtsolut':

                print_output_1 = ""
                print_output_2 = ""
                print_output_3 = ""
                print_output_4 = ""
                if not os.path.exists(os.path.join('MXX_P_MYY', 'OUTPUT_MXX_P_MYY')):
                    print 'WARNING: %s is not created! wait for 10 seconds and re-check.' \
                                %(os.path.join('MXX_P_MYY', 'OUTPUT_MXX_P_MYY'))
                    time.sleep(10)
                    if not os.path.exists(os.path.join('MXX_P_MYY', 'OUTPUT_MXX_P_MYY')):
                        sys.exit('...ERROR...%s is not created.'
                                    %(os.path.join('MXX_P_MYY', 'OUTPUT_MXX_P_MYY')))

                if not os.path.exists(os.path.join('MXY_MXX_M_MYY', 'OUTPUT_MXY_MXX_M_MYY')):
                    print 'WARNING: %s is not created! wait for 10 seconds and re-check.' \
                                %(os.path.join('MXY_MXX_M_MYY', 'OUTPUT_MXY_MXX_M_MYY'))
                    time.sleep(10)
                    if not os.path.exists(os.path.join('MXY_MXX_M_MYY', 'OUTPUT_MXY_MXX_M_MYY')):
                        sys.exit('...ERROR...%s is not created.'
                                    %(os.path.join('MXY_MXX_M_MYY', 'OUTPUT_MXY_MXX_M_MYY')))

                if not os.path.exists(os.path.join('MXZ_MYZ', 'OUTPUT_MXZ_MYZ')):
                    print 'WARNING: %s is not created! wait for 10 seconds and re-check.' \
                                %(os.path.join('MXZ_MYZ', 'OUTPUT_MXZ_MYZ'))
                    time.sleep(10)
                    if not os.path.exists(os.path.join('MXZ_MYZ', 'OUTPUT_MXZ_MYZ')):
                        sys.exit('...ERROR...%s is not created.'
                                    %(os.path.join('MXZ_MYZ', 'OUTPUT_MXZ_MYZ')))

                if not os.path.exists(os.path.join('MZZ', 'OUTPUT_MZZ')):
                    print 'WARNING: %s is not created! wait for 10 seconds and re-check.' \
                                %(os.path.join('MZZ', 'OUTPUT_MZZ'))
                    time.sleep(10)
                    if not os.path.exists(os.path.join('MZZ', 'OUTPUT_MZZ')):
                        sys.exit('...ERROR...%s is not created.'
                                    %(os.path.join('MZZ', 'OUTPUT_MZZ')))

                while ((test_1 == -1 or test_2 == -1 or test_3 == -1 or test_4 == -1) and solverended == 0):
                    output_file =  open(os.path.join('MXX_P_MYY', 'OUTPUT_MXX_P_MYY'), 'r')
                    output_file_read = output_file.readlines()
                    last_output_1 = print_output_1
                    print_output_1 = output_file_read[-1].split('\n')[0]
                    if output_file_read[-1].find('PROGRAM axisem FINISHED') == -1:
                        if (last_output_1 == print_output_1) :
                            same_output_1 = same_output_1 + 1
                        else :
                            same_output_1 = 0
                        if (same_output_1 >= 30) :
                            print "\n============================="
                            print '\n       SOLVER froze'
                            print '\n         MXX_P_MYY   '
                            print '\n      killing axisem'
                            print '\n        I love it!'
                            print "\n============================="
                            temp = killprocess('axisem')
                            return
                        test_1 = -1
                    else:
                        test_1 = 0
                    print 'MXX_P_MYY:     ' + print_output_1
                    output_file.close

                    output_file =  open(os.path.join('MXY_MXX_M_MYY', 'OUTPUT_MXY_MXX_M_MYY'), 'r')
                    output_file_read = output_file.readlines()
                    last_output_2 = print_output_2
                    print_output_2 = output_file_read[-1].split('\n')[0]
                    if output_file_read[-1].find('PROGRAM axisem FINISHED') == -1:
                        # Check whether Output is stagnant since more than one minute.
                        if (last_output_2 == print_output_2) :
                            same_output_2 = same_output_2 + 1
                        else :
                            same_output_2 = 0
                        if (same_output_2 >= 30) :
                            print "\n============================="
                            print '\n       SOLVER froze'
                            print '\n      MXY_MXX_M_MYY'
                            print '\n      killing axisem'
                            print '\n        I love it!'
                            print "\n============================="
                            temp = killprocess('axisem')
                            return
                        test_2 = -1
                    else:
                        test_2 = 0
                    print 'MXY_MXX_M_MYY: ' + print_output_2
                    output_file.close

                    output_file =  open(os.path.join('MXZ_MYZ', 'OUTPUT_MXZ_MYZ'), 'r')
                    output_file_read = output_file.readlines()
                    last_output_3 = print_output_3
                    print_output_3 = output_file_read[-1].split('\n')[0]
                    if output_file_read[-1].find('PROGRAM axisem FINISHED') == -1:
                        test_3 = -1
                        # Check whether Output is stagnant since more than one minute.
                        if (last_output_3 == print_output_3) :
                            same_output_3 = same_output_3 + 1
                        else :
                            same_output_3 = 0
                        if (same_output_3 >= 30) :
                            print "\n============================="
                            print '\n       SOLVER froze'
                            print '\n          MXZ_MYZ   '
                            print '\n      killing axisem'
                            print '\n        I love it!'
                            print "\n============================="
                            temp = killprocess('axisem')
                            return
                    else:
                        test_3 = 0
                    print 'MXZ_MYZ:       ' + print_output_3
                    output_file.close

                    output_file =  open(os.path.join('MZZ', 'OUTPUT_MZZ'), 'r')
                    output_file_read = output_file.readlines()
                    last_output_4 = print_output_4
                    print_output_4 = output_file_read[-1].split('\n')[0]
                    if output_file_read[-1].find('PROGRAM axisem FINISHED') == -1:
                        # Check whether Output is stagnant since more than one minute.
                        if (last_output_4 == print_output_4) :
                            same_output_4 = same_output_4 + 1
                        else :
                            same_output_4 = 0
                        if (same_output_4 >= 30) :
                            print "\n============================="
                            print '\n       SOLVER froze'
                            print '\n            MZZ'
                            print '\n      killing axisem'
                            print '\n        I love it!'
                            print "\n============================="
                            temp = killprocess('axisem')
                            return
                        test_4 = -1
                    else:
                        test_4 = 0
                    print 'MZZ:           ' + print_output_4
                    output_file.close

                    if howmanyofthisprocess('axisem') == 0:
                        solverended = 1
                    time.sleep(2)
                    print '--------------------------------'

                if (test_1 == -1 or test_2 == -1 or test_3 == -1 or test_4 == -1):
                    solvercrashed = 1
    t2_solver = time.time()

    if solvercrashed == 1:
        print "\n============================="
        print '\n       SOLVER crashed'
        print "\n============================="
        raise RuntimeError('Solver crashed')
        return

    ##############################################################
    ######################## Post-Processing #####################
    ##############################################################
    t1_post = time.time()

    if input['post_processing'] != 'N':
        if input['verbose'] != 'N':
            print "\n==============="
            print "Post Processing"
            print "===============\n"
            stdout_param = None
        else:
            sys.stdout.write('\nPost Processing...')
            sys.stdout.flush()
            stdout_param = subprocess.PIPE

        if os.path.isdir(os.path.join(input['axi_address'], 'SOLVER',
                            input['solver_name'], 'Data_Postprocessing')):
            subprocess.check_call(['rm', '-rf', os.path.join(input['axi_address'], 'SOLVER',
                            input['solver_name'], 'Data_Postprocessing')])

        if input['sourcefile_type'] == 'sourceparams':
            os.chdir(os.path.join(input['axi_address'], 'SOLVER', input['solver_name']))
            post_process_open = open('./param_post_processing', 'r')
            post_process_read_fio = post_process_open.readlines()
            post_process_read = []
            for post_r in post_process_read_fio:
                if not post_r.startswith('#'):
                    if not post_r.startswith('\n'):
                        post_process_read.append(post_r)
            # Call edit_param_post_processing to change the post_processing options
            post_process_read = edit_param_post_processing(post_process_read = post_process_read)
            post_process_open.close()
            post_process_open = open('./param_post_processing', 'w')
            for i in range(0, len(post_process_read)):
                post_process_open.write(post_process_read[i])
                if input['verbose'] != 'N': print post_process_read[i].split('\n')[0]
            if input['verbose'] != 'N': print 2*"======================================"
            post_process_open.close()
            output = subprocess.check_call(['./post_processing.csh'], stdout = stdout_param)
            if output != 0: print output_print

        elif input['sourcefile_type'] == 'cmtsolut':
            os.chdir(os.path.join(input['axi_address'], 'SOLVER', input['solver_name']))
            post_process_open = open('./param_post_processing', 'r')
            post_process_read_fio = post_process_open.readlines()
            post_process_read = []
            for post_r in post_process_read_fio:
                if not post_r.startswith('#'):
                    if not post_r.startswith('\n'):
                        post_process_read.append(post_r)
            # Call edit_param_post_processing to change the post_processing options
            post_process_read = edit_param_post_processing(post_process_read = post_process_read)
            post_process_open.close()
            post_process_open = open('./param_post_processing', 'w')
            for i in range(0, len(post_process_read)):
                post_process_open.write(post_process_read[i])
                if input['verbose'] != 'N': print post_process_read[i].split('\n')[0]
            if input['verbose'] != 'N': print 2*"======================================"
            post_process_open.close()
            output = subprocess.check_call(['./post_processing.csh'], stdout = stdout_param)
            if output != 0: print output_print

        print 'DONE'

    t2_post = time.time()

    ##############################################################
    ############################# MSEED ##########################
    ##############################################################
    t1_misc = time.time()

    if input['mseed'] != 'N':
        if input['verbose'] != 'N':
            print "\n============================="
            print "Creating MSEED files"
            print "one MSEED for each seismogram"
            print "=============================\n"
        else:
            sys.stdout.write('Creating MSEED...')
            sys.stdout.flush()
        path = os.path.join(input['axi_address'], 'SOLVER', input['solver_name'],
                                'Data_Postprocessing', 'SEISMOGRAMS')
        axisem2mseed(path = path)
        print 'DONE'

    if input['mseed_all'] != 'N':
        if input['verbose'] != 'N':
            print "\n================================================"
            print "Creating MSEED files"
            print "one MSEED for ALL seismogram (seismograms.mseed)"
            print "================================================\n"
        else:
            sys.stdout.write('Creating MSEED...')
            sys.stdout.flush()
        path = os.path.join(input['axi_address'], 'SOLVER', input['solver_name'],
                                'Data_Postprocessing', 'SEISMOGRAMS')
        axisem2mseed_all(path = path)
        print 'DONE'

    t2_misc = time.time()

    ##############################################################
    ############################# TEST ###########################
    ##############################################################
    t1_test = time.time()

    if input['test'] != 'N' or input['plot_test'] != 'N':

        os.chdir(os.path.join(input['axi_address'], 'TESTING'))
        chans = eval(input['test_chans'])
        fmin = eval(input['test_fmin'])
        fmax = eval(input['test_fmax'])
        halfduration = eval(input['test_half'])
        folder = os.path.join(input['test_folder'], 'ref_data')

        if os.path.isdir(os.path.join(input['test_folder'], 'new_data')):
            shutil.rmtree(os.path.join(input['test_folder'], 'new_data'))
        os.makedirs(os.path.join(input['test_folder'], 'new_data'))
        folder_new = os.path.join(input['test_folder'], 'new_data')

        subprocess.check_call(['cp',
                os.path.join(input['axi_address'], 'SOLVER', input['solver_name'],
                                'Data_Postprocessing', 'SEISMOGRAMS', 'seismograms.mseed'),
                                os.path.join(folder_new, 'axisem.mseed')])

        sgs = []
        st = read(os.path.join(folder, 'axisem.mseed'))
        sgs.append(st)
        if not os.path.isfile(os.path.join(folder, 'yspec.mseed')):
            print "You shoud provide 'yspec.mseed' for testing!"
        st = read(os.path.join(folder, 'yspec.mseed'))
        sgs.append(st)
        st = read(os.path.join(folder_new, 'axisem.mseed'))
        sgs.append(st)

        sigma =  halfduration / np.sqrt(2.) / 3.5

        for st in sgs:
            convSTF(st, sigma=sigma)
        for sg in sgs:
            sg.filter('lowpass', freq=fmax, corners=2)
            sg.filter('highpass', freq=fmin, corners=2)

        nstat = eval(input['test_nstat'])
        stats = np.arange(nstat) + 1

        labels = ['Axisem-ref', 'Yspec', 'Axisem-new']
        colors = ['r', 'gray', 'b']
        linestyles = ['-', '-', '-']

        t = []

        i = 0
        t0 = sgs[i][0].stats.starttime
        dt = sgs[i][0].stats.delta
        npts = sgs[i][0].stats.npts
        t.append(np.linspace(0., dt * (npts -1), npts) + (t0 -  UTCDateTime(0)))

        i = 1
        t0 = sgs[i][0].stats.starttime
        dt = sgs[i][0].stats.delta
        npts = sgs[i][0].stats.npts
        t.append(np.linspace(0., dt * (npts -1), npts) + (t0 -  UTCDateTime(0)))

        i = 2
        t0 = sgs[i][0].stats.starttime
        dt = sgs[i][0].stats.delta
        npts = sgs[i][0].stats.npts
        t.append(np.linspace(0., dt * (npts -1), npts) + (t0 -  UTCDateTime(0)))

        misfitplot = plt.figure()
        for chan in chans:
            recsec = plt.figure()
            ax = recsec.gca()
            maxi = 0.
            # find global maxium
            for stat in stats:
                str_stat = '*%02d' % stat
                tr = sgs[0].select(station=str_stat, channel='*'+chan)[0]
                dat1 = tr.data
                maxi = max(maxi, np.abs(dat1).max())

            l2misfit = []
            for n in stats:
                stat = '*%02d' % (n,)
                for i, sg in enumerate(sgs):
                    dat = sg.select(station=stat, channel='*'+chan)[0].data
                    #if i == 0:
                    #    maxl = dat.max()
                    dat = dat / maxi * 2.
                    #dat = dat / maxl / 2.
                    if n == 1:
                        ax.plot(t[i], dat + n, colors[i], label=labels[i], ls=linestyles[i])
                    else:
                        ax.plot(t[i], dat + n, colors[i], label='_nolegend_', ls=linestyles[i])
                # compute l2 misfits
                dat1 = sgs[0].select(station=stat, channel='*'+chan)[0].data
                dat2 = sgs[2].select(station=stat, channel='*'+chan)[0].data
                l2misfit.append(((dat1 - dat2)**2).sum()**.5 / maxi /
                       sgs[0][0].stats.npts)

            # write l2 misfits to file
            fl2 = open(os.path.join(folder_new, 'l2misfit_%s.dat' % chan), 'w')
            for l2 in l2misfit:
                fl2.write("%4.2e\n" % l2)
            fl2.close()

            ax.set_xlabel('time / seconds')
            ax.set_xlim(0, 1800)
            ax.set_ylim(0, nstat + 1)
            ax.legend()

            recsec.set_size_inches((16,12))
            recsec.suptitle(chan)
            recsec.subplots_adjust(hspace=0.3, wspace=0.2, left=0.08, right=0.95,
                    top=0.93, bottom=0.07)

            if input['save_plots_test'] != 'N':
                recsec.savefig(os.path.join(folder_new, 'record_section_%s.%s' %
                        (chan, input['plot_format_test'])))

            ax = misfitplot.gca()
            ax.semilogy(np.arange(nstat) + 1, np.array(l2misfit) + 1e-12, 'o', label=chan)
            ax.set_xlabel('trace')
            ax.set_ylabel('l2 - misfit to reference data')
            ax.set_ylim(1e-12, 1.)

            ax.axhline(y=1e-5, color='k', ls='--')
            if np.max(l2misfit) > 1e-5:
                fwarn = open(os.path.join(folder_new, 'warning.dat'), 'a')
                fwarn.write("maximum l2 norm misfit larger then 1e-5 in chan %s trace %d\n"
                             % (chan, np.argmax(l2misfit)))
                fwarn.close()

        ax = misfitplot.gca()
        ax.legend()
        if input['save_plots_test'] != 'N':
            misfitplot.savefig(os.path.join(folder_new, 'l2_misfit.' + input['plot_format_test']))
        if input['plot_test'] != 'N':
            plt.show()

    t2_test = time.time()

    print "\n============================================================"
    print "AXISEM statistics:"
    print "Time for MESHER: " + str(t2_mesher - t1_mesher)
    print "Time for SOLVER: " + str(t2_solver - t1_solver)
    print "Time for POST  : " + str(t2_post - t1_post)
    print "Time for MSEED : " + str(t2_misc - t1_misc)
    if input['test'] != 'N':
        print "Time for TEST  : " + str(t2_test - t1_test)
    else:
        print "Time for TEST  : NA"
    print "============================================================"


###################### read_input_file #################################

def read_input_file():
    """
    Read inputs from inpython.cfg file.
    """

    global input, obspy_error
    config = ConfigParser.RawConfigParser()
    input = {}
    try:
        config.read(os.path.join(sys.argv[1]))
        input['inpython_address'] = os.path.join(sys.argv[1])
        print '\n****************************************'
        print 'Read the input file (inpython.cfg) from:'
        print os.path.join(sys.argv[1])
    except Exception, error:
        config.read(os.path.join(os.getcwd(), 'inpython.cfg'))
        input['inpython_address'] = os.path.join(os.getcwd(), 'inpython.cfg')
        print '\n****************************************'
        print 'Read the input file (inpython.cfg) from:'
        print os.path.join(os.getcwd(), 'inpython.cfg')

    try:
        input['IO_STATION'] = os.path.join(sys.argv[2])
        print '\n****************************************'
        print 'Copy the STATIONS file from:'
        print os.path.join(sys.argv[2])
        if not os.path.isabs(input['IO_STATION']):
            input['IO_STATION'] = os.path.join(os.getcwd(),
                                          input['IO_STATION'])
    except Exception, error:
        print 'STATIONS file is not entered as an input!'
        input['IO_STATION'] = False

    if not os.path.isabs(input['inpython_address']):
        input['inpython_address'] = os.path.join(os.getcwd(),
                                      input['inpython_address'])

    input['axi_address'] = config.get('GENERAL', 'AXISEM_DIR')
    if not os.path.isabs(input['axi_address']):
        input['axi_address'] = os.path.join(os.getcwd(), input['axi_address'])
    print '\nWorking Directory:'
    print input['axi_address']
    print '****************************************\n'

    input['solver_name'] = config.get('GENERAL', 'SOLVER_NAME')
    input['verbose'] = config.get('GENERAL', 'VERBOSE')

    input['new_mesh'] = config.get('GENERAL', 'NEW_MESH')
    input['post_processing'] = config.get('GENERAL', 'POST_PROCESSING')

    input['mesher'] = config.get('GENERAL', 'MESHER')
    input['solver'] = config.get('GENERAL', 'SOLVER')

    input['mesher_makefile'] = config.get('GENERAL', 'MESHER_MAKEFILE')
    input['mesher_make'] = config.get('GENERAL', 'MESHER_MAKE')
    input['mesher_move'] = config.get('GENERAL', 'MESHER_MOVE')

    input['solver_cp'] = config.get('GENERAL', 'SOLVER_COPY')
    input['solver_makefile'] = config.get('GENERAL', 'SOLVER_MAKEFILE')
    input['solver_make'] = config.get('GENERAL', 'SOLVER_MAKE')

    if input['mesher_make'] != 'N':
        input['mesher_move'] = 'Y'
        input['solver_cp'] = 'Y'

    input['mesher_bg_model'] = config.get('MESHER_BASIC', 'BACKGROUND_MODEL')
    input['mesher_ext_model'] = config.get('MESHER_BASIC', 'EXT_MODEL')
    input['mesher_dominant_period'] = config.get('MESHER_BASIC', 'DOMINANT_PERIOD')
    input['mesher_ntheta'] = config.get('MESHER_BASIC', 'NTHETA_SLICES')
    input['mesher_nradial'] = config.get('MESHER_BASIC', 'NRADIAL_SLICES')
    input['mesher_write_vtk'] = config.get('MESHER_BASIC', 'WRITE_VTK')
    input['mesher_coarsening_layers'] = config.get('MESHER_BASIC', 'COARSENING_LAYERS')
    input['mesh_name'] = config.get('MESHER_BASIC', 'MESHNAME')
    input['mesher_ic_shear_wave'] = config.get('MESHER_ADVANCED', 'IC_SHEAR_WAVE')
    input['mesher_npol'] = config.get('MESHER_ADVANCED', 'NPOL')
    input['mesher_el_per_lambda'] = config.get('MESHER_ADVANCED', 'EL_PER_LAMBDA')
    input['mesher_courant_nr'] = config.get('MESHER_ADVANCED', 'COURANT_NR')
    input['mesher_radius'] = config.get('MESHER_ADVANCED', 'RADIUS')
    input['mesher_save_mesh'] = config.get('MESHER_ADVANCED', 'SAVE_MESH')
    input['mesher_verbose'] = config.get('MESHER_ADVANCED', 'VERBOSE')

    input['solver_sim_type'] = config.get('SOLVER_BASIC', 'SIMULATION_TYPE')
    if input['solver_sim_type'] == 'single':
        input['sourcefile_type'] = 'sourceparams'
    elif input['solver_sim_type'] == 'moment':
        input['sourcefile_type'] = 'cmtsolut'
    else:
        print 'Check your simulation type, you entered:'
        print input['solver_simu_type']
    input['solver_seis_length'] = config.get('SOLVER_BASIC', 'SEISMOGRAM_LENGTH')
    input['solver_recfile_type'] = config.get('STATION_INFO', 'RECFILE_TYPE')
    input['solver_lat_heterogeneity'] = config.get('SOLVER_BASIC', 'LAT_HETEROGENEITY')
    input['solver_attenuation'] = config.get('SOLVER_BASIC', 'ATTENUATION')
    input['solver_save_snapshots'] = config.get('SOLVER_BASIC', 'SAVE_SNAPSHOTS')
    input['solver_verbosity'] = config.get('SOLVER_BASIC', 'VERBOSITY')


    input['solver_sampling_rate'] = config.get('SOLVER_ADVANCED', 'SAMPLING_PERIOD')
    input['solver_time_step'] = config.get('SOLVER_ADVANCED', 'TIME_STEP')
    input['solver_source_period'] = config.get('SOLVER_ADVANCED', 'SOURCE_PERIOD')
    input['solver_time_scheme'] = config.get('SOLVER_ADVANCED', 'TIME_SCHEME')
    input['solver_data_dir'] = config.get('SOLVER_ADVANCED', 'DATA_DIR')
    input['solver_info_dir'] = config.get('SOLVER_ADVANCED', 'INFO_DIR')
    input['solver_diag_file_output'] = config.get('SOLVER_ADVANCED', 'DIAGNOSTIC_FILE_OUTPUT')
    input['solver_mesh_test'] = config.get('SOLVER_ADVANCED', 'MESH_TEST')
    input['solver_deflate_level'] = config.get('SOLVER_ADVANCED', 'DEFLATE_LEVEL')
    input['solver_snapshot_dt'] = config.get('SOLVER_ADVANCED', 'SNAPSHOT_DT')
    input['solver_snapshots_format'] = config.get('SOLVER_ADVANCED', 'SNAPSHOTS_FORMAT')
    input['solver_kernel_wavefields'] = config.get('SOLVER_ADVANCED', 'KERNEL_WAVEFIELDS')
    input['solver_kernel_spp'] = config.get('SOLVER_ADVANCED', 'KERNEL_SPP')
    input['solver_kernel_source'] = config.get('SOLVER_ADVANCED', 'KERNEL_SOURCE')
    input['solver_kernel_ibeg'] = config.get('SOLVER_ADVANCED', 'KERNEL_IBEG')
    input['solver_kernel_iend'] = config.get('SOLVER_ADVANCED', 'KERNEL_IEND')
    input['solver_nr_lin_solids'] = config.get('SOLVER_ADVANCED', 'NR_LIN_SOLIDS')
    input['solver_fmin'] = config.get('SOLVER_ADVANCED', 'F_MIN')
    input['solver_fmax'] = config.get('SOLVER_ADVANCED', 'F_MAX')
    input['solver_fref'] = config.get('SOLVER_ADVANCED', 'F_REFERENCE')
    input['solver_small_q_correction'] = config.get('SOLVER_ADVANCED', 'SMALL_Q_CORRECTION')
    input['solver_nr_f_sample'] = config.get('SOLVER_ADVANCED', 'NR_F_SAMPLE')
    input['solver_maxint_sa'] = config.get('SOLVER_ADVANCED', 'MAXINT_SA')
    input['solver_tstart_sr'] = config.get('SOLVER_ADVANCED', 'TSTART_SR')
    input['solver_tstart_amp'] = config.get('SOLVER_ADVANCED', 'TSTART_AMP')
    input['solver_t_decay'] = config.get('SOLVER_ADVANCED', 'T_DECAY')
    input['solver_fix_freq'] = config.get('SOLVER_ADVANCED', 'FIX_FREQ')
    input['solver_dump_vtk'] = config.get('SOLVER_ADVANCED', 'DUMP_VTK')
    input['solver_coarse_grained'] = config.get('SOLVER_ADVANCED', 'COARSE_GRAINED')
    input['solver_save_energy'] = config.get('SOLVER_ADVANCED', 'SAVE_ENERGY')
    input['solver_homo_model'] = config.get('SOLVER_ADVANCED', 'HOMO_MODEL')
    input['solver_homo_vp'] = config.get('SOLVER_ADVANCED', 'HOMO_VP')
    input['solver_homo_vs'] = config.get('SOLVER_ADVANCED', 'HOMO_VS')
    input['solver_homo_rho'] = config.get('SOLVER_ADVANCED', 'HOMO_RHO')
    input['solver_force_aniso'] = config.get('SOLVER_ADVANCED', 'FORCE_ANISO')

    input['receiver_type'] = input['solver_recfile_type']
    input['source_type'] = config.get('SOURCE_INFO', 'SOURCE_TYPE')

    input['source_depth'] = config.get('SOURCE_INFO', 'SOURCE_DEPTH')
    input['source_lat'] = config.get('SOURCE_INFO', 'SOURCE_LATITUDE')
    input['source_lon'] = config.get('SOURCE_INFO', 'SOURCE_LONGITUDE')
    input['source_stf'] = config.get('SOURCE_INFO', 'SOURCE_STF')
    input['source_amp'] = config.get('SOURCE_INFO', 'SOURCE_AMPLITUDE')

    input['cmt_lat']   = config.get('SOURCE_INFO', 'CMT_LAT')
    input['cmt_lon']   = config.get('SOURCE_INFO', 'CMT_LON')
    input['cmt_depth'] = config.get('SOURCE_INFO', 'CMT_DEPTH')
    input['cmt_Mrr']   = config.get('SOURCE_INFO', 'CMT_MRR')
    input['cmt_Mtt']   = config.get('SOURCE_INFO', 'CMT_MTT')
    input['cmt_Mpp']   = config.get('SOURCE_INFO', 'CMT_MPP')
    input['cmt_Mrt']   = config.get('SOURCE_INFO', 'CMT_MRT')
    input['cmt_Mrp']   = config.get('SOURCE_INFO', 'CMT_MRP')
    input['cmt_Mtp']   = config.get('SOURCE_INFO', 'CMT_MTP')

    input['post_rec_comp_sys'] = config.get('POST_PROCESSING', 'REC_COMP_SYS')
    input['post_conv_period'] = config.get('POST_PROCESSING', 'CONV_PERIOD')
    input['post_conv_stf'] = config.get('POST_PROCESSING', 'CONV_STF')
    input['post_seistype'] = config.get('POST_PROCESSING', 'SEISTYPE')
    input['post_load_snaps'] = config.get('POST_PROCESSING', 'LOAD_SNAPS')
    input['post_data_dir'] = config.get('POST_PROCESSING', 'DATA_DIR')
    input['post_negative_time'] = config.get('POST_PROCESSING', 'NEGATIVE_TIME')
    input['post_3D_phi_start'] = config.get('POST_PROCESSING', '3D_PHI_START')
    input['post_3D_phi_end'] = config.get('POST_PROCESSING', '3D_PHI_END')
    input['post_3D_rtop'] = config.get('POST_PROCESSING', '3D_RTOP')
    input['post_3D_rbot'] = config.get('POST_PROCESSING', '3D_RBOT')
    input['post_3D_plot_top'] = config.get('POST_PROCESSING', '3D_PLOT_TOP')
    input['post_3D_plot_bot'] = config.get('POST_PROCESSING', '3D_PLOT_BOT')
    input['post_3D_snap_beg'] = config.get('POST_PROCESSING', '3D_SNAP_BEG')
    input['post_3D_snap_end'] = config.get('POST_PROCESSING', '3D_SNAP_END')
    input['post_3D_snap_stride'] = config.get('POST_PROCESSING', '3D_SNAP_STRIDE')

    input['mseed'] = config.get('MISC', 'MSEED')
    input['mseed_all'] = config.get('MISC', 'MSEED_ALL')
    input['convSTF'] = config.get('MISC', 'CONV_STF')
    input['halfduration'] = eval(config.get('MISC', 'HALF_DURATION'))
    input['filter'] = config.get('MISC', 'FILTER')
    input['fmin'] = eval(config.get('MISC', 'FMIN'))
    input['fmax'] = eval(config.get('MISC', 'FMAX'))

    input['make_axisem_mpirun'] = config.get('MAKE_AXISEM', 'MPIRUN')
    input['make_axisem_use_netcdf'] = config.get('MAKE_AXISEM', 'USE_NETCDF')
    input['make_axisem_netcdf_path'] = config.get('MAKE_AXISEM', 'NETCDF_PATH')
    input['make_axisem_serial'] = config.get('MAKE_AXISEM', 'SERIAL')
    input['make_axisem_include_mpi'] = config.get('MAKE_AXISEM', 'INCLUDE_MPI')
    input['make_axisem_CC'] = config.get('MAKE_AXISEM', 'CC')
    input['make_axisem_FC'] = config.get('MAKE_AXISEM', 'FC')
    input['make_axisem_FFLAGS'] = config.get('MAKE_AXISEM', 'FFLAGS')
    input['make_axisem_CFLAGS'] = config.get('MAKE_AXISEM', 'CFLAGS')
    input['make_axisem_LDFLAGS'] = config.get('MAKE_AXISEM', 'LDFLAGS')

    input['test'] = config.get('TEST', 'TEST')
    input['test_folder'] = config.get('TEST', 'TEST_FOLDER')
    input['plot_test'] = config.get('TEST', 'PLOT')
    input['save_plots_test'] = config.get('TEST', 'SAVE_PLOTS')
    input['plot_format_test'] = config.get('TEST', 'PLOT_FORMAT')
    input['test_chans'] = config.get('TEST', 'CHANS')
    input['test_fmin'] = config.get('TEST', 'FMIN')
    input['test_fmax'] = config.get('TEST', 'FMAX')
    input['test_half'] = config.get('TEST', 'HALF_DURATION')
    input['test_nstat'] = config.get('TEST', 'NSTAT')

    if input['new_mesh'] == 'Y':
        input['mesher'] = 'Y'
        input['solver'] = 'Y'
        input['mesher_makefile'] = 'Y'
        input['mesher_make'] = 'Y'
        input['mesher_move'] = 'Y'
        input['solver_cp'] = 'Y'
        input['solver_makefile'] = 'Y'
        input['solver_make'] = 'Y'

    if input['new_mesh'] == 'N':
        input['mesher'] = 'N'
        input['solver'] = 'Y'
        input['solver_cp'] = 'N'
        input['solver_makefile'] = 'N'
        input['solver_make'] = 'Y'

    if input['test'] == 'N':
        input['plot_test'] = 'N'
        input['save_plots_test'] = 'N'

    if input['plot_test'] != 'N':
        input['mesher'] = 'N'
        input['solver'] = 'N'
        input['solver_cp'] = 'N'
        input['solver_makefile'] = 'N'
        input['solver_make'] = 'N'
        input['post_processing'] = 'N'
        input['test'] = 'N'
        input['mseed'] = 'N'

        print '##################################'
        print "PyAxi tries to copy the data from:"
        print "(solver_name flag in inpython.cfg)"
        print input['solver_name']
        print '##################################'

    if input['receiver_type'] == 'database':
        input['post_processing'] = 'N'

    if obspy_error != 'N':
        input['test'] = 'N'
        input['mseed'] = 'N'

###################### edit_param_post_processing ######################
def edit_param_post_processing(post_process_read):
    """
    edit param_post_processing file
    """
    global input

    if input['post_rec_comp_sys'] != 'DNC':
        post_process_read[0] = \
            'REC_COMP_SYS    %s\n' %(input['post_rec_comp_sys'])
    if input['post_conv_period'] != 'DNC':
        post_process_read[1] = \
            'CONV_PERIOD     %s\n' %(input['post_conv_period'])
    if input['post_conv_stf'] != 'DNC':
        post_process_read[2] = \
            'CONV_STF        %s\n' %(input['post_conv_stf'])
    if input['post_seistype'] != 'DNC':
        post_process_read[3] = \
            'SEISTYPE        %s\n' %(input['post_seistype'])
    if input['post_load_snaps'] != 'DNC':
        post_process_read[4] = \
            'LOAD_SNAPS      %s\n' %(input['post_load_snaps'])
    if input['post_data_dir'] != 'DNC':
        post_process_read[5] = \
            'DATA_DIR        %s\n' %(input['post_data_dir'])
    if input['post_negative_time'] != 'DNC':
        post_process_read[6] = \
            'NEGATIVE_TIME   %s\n' %(input['post_negative_time'])
    if input['post_3D_phi_start'] != 'DNC':
        post_process_read[7] = \
            '3D_PHI_START     %s\n' %(input['post_3D_phi_start'])
    if input['post_3D_phi_end'] != 'DNC':
        post_process_read[8] = \
            '3D_PHI_END      %s\n' %(input['post_3D_phi_end'])
    if input['post_3D_rtop'] != 'DNC':
        post_process_read[9] = \
            '3D_RTOP         %s\n' %(input['post_3D_rtop'])
    if input['post_3D_rbot'] != 'DNC':
        post_process_read[10] = \
            '3D_RBOT         %s\n' %(input['post_3D_rbot'])
    if input['post_3D_plot_top'] != 'DNC':
        post_process_read[11] = \
            '3D_PLOT_TOP     %s\n' %(input['post_3D_plot_top'])
    if input['post_3D_plot_bot'] != 'DNC':
        post_process_read[12] = \
            '3D_PLOT_BOT     %s\n' %(input['post_3D_plot_bot'])
    if input['post_3D_snap_beg'] != 'DNC':
        post_process_read[13] = \
            '3D_SNAP_BEG      %s\n' %(input['post_3D_snap_beg'])
    if input['post_3D_snap_end'] != 'DNC':
        post_process_read[14] = \
            '3D_SNAP_END      %s\n' %(input['post_3D_snap_end'])
    if input['post_3D_snap_stride'] != 'DNC':
        post_process_read[15] = \
            '3D_SNAP_STRIDE   %s\n' %(input['post_3D_snap_stride'])

    return post_process_read

########################## axisem2mseed ################################
def axisem2mseed(path):
    """
    change .dat files into MSEED format
    """
    global input

    if not os.path.isdir(os.path.join(path, 'MSEED')):
        os.mkdir(os.path.join(path, 'MSEED'))
    else:
        print 'Following directory already exists:'
        print os.path.join(path, 'MSEED')
        sys.exit()

    t = UTCDateTime(0)
    traces = []

    for file in glob.iglob(os.path.join(path, '*.dat')):
        stationID = file.split('/')[-1].split('_')[0]
        networkID = file.split('/')[-1].split('_')[1]
        chan = file.split('/')[-1].split('_')[-1].split('.')[0]
        if chan == 'E': chan = 'BHE'
        elif chan == 'N': chan = 'BHN'
        elif chan == 'Z': chan = 'BHZ'
        try:
            dat = np.loadtxt(file)
            npts = len(dat[:,0])
            stats = {'network': networkID,
                     'station': stationID,
                     'location': '',
                     'channel': chan,
                     'npts': npts,
                     'sampling_rate': (npts - 1.)/(dat[-1,0] - dat[0,0]),
                     'starttime': t + dat[0,0],
                     'mseed' : {'dataquality': 'D'}}
            st = Stream(Trace(data=dat[:,1], header=stats))
            if input['convSTF'] == 'Y':
                sigma =  input['halfduration'] / np.sqrt(2.) / 3.5
                convSTF(st, sigma=sigma)
            if input['filter'] == 'Y':
                st.filter('lowpass', freq=input['fmax'], corners=2)
                st.filter('lowpass', freq=input['fmax'], corners=2)
                st.filter('lowpass', freq=input['fmax'], corners=2)
                st.filter('lowpass', freq=input['fmax'], corners=2)
                st.filter('highpass', freq=input['fmin'], corners=2)
                st.filter('highpass', freq=input['fmin'], corners=2)
            fname =  os.path.join(path, 'MSEED', 'dis.' + stationID + '..' + chan)
            st.write(fname, format='MSEED')
        except Exception, e:
            print e
            print networkID + '.' + stationID + '.' + chan + '.mseed'
            print '-------------------------------------------------'

########################## axisem2mseed_all ###########################
def axisem2mseed_all(path):
    """
    change .dat files into MSEED format
    """
    global input

    t = UTCDateTime(0)
    traces = []

    for file in glob.iglob(os.path.join(path, '*.dat')):
        stationID = file.split('/')[-1].split('_')[0]
        chan = file.split('/')[-1].split('_')[-1].split('.')[0]
        dat = np.loadtxt(file)
        npts = len(dat[:,0])
        stats = {'network': 'SG',
                 'station': stationID,
                 'location': '',
                 'channel': chan,
                 'npts': npts,
                 'sampling_rate': (npts - 1.)/(dat[-1,0] - dat[0,0]),
                 'starttime': t + dat[0,0],
                 'mseed' : {'dataquality': 'D'}}
        traces.append(Trace(data=dat[:,1], header=stats))

    st = Stream(traces)
    st.sort()
    fname =  os.path.join(path, 'seismograms.mseed')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        st.write(fname, format='MSEED')

########################## convSTF ################################
def convSTF(st, sigma=30.):

    gauss = lambda (t, s): 1. / (2. * np.pi * s**2.)**.5 \
                           * np.exp(-1*(t**2)/(2*(s**2)))

    df = st[0].stats.sampling_rate
    dt = 1./df

    t = np.linspace(0., sigma * 20., sigma * 20 * df + 1)
    stf = gauss((t-sigma*10, sigma))
    nstf = len(stf)

    for tr in st:
        tr.data *= cosine_taper(len(tr.data), p=0.05)
        nfft = util.next_pow_2(max(nstf, tr.stats.npts)) * 2
        stff = np.fft.rfft(stf, n=nfft) * dt
        trf = np.fft.rfft(tr, n=nfft) * dt
        tr.data = np.fft.irfft(stff * trf)[sigma*10*df:sigma*10*df+len(tr.data)] * df

    return 1

########################## isProcessRunning #############################

def howmanyofthisprocess( process_name ):
    tmp = os.popen("ps ").read()
    proccount = tmp.count(process_name)
    return proccount

########################## killprocess #################################

def killprocess( process_name ):
    os.system("killall "+process_name);
    time.sleep(10)
    os.system("killall -9 "+process_name);

########################## manual_which ################################
def manual_which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

########################################################################
########################################################################
########################################################################

if __name__ == "__main__":

    t1_pro = time.time()
    status = PyAxi()
    t_pro = time.time() - t1_pro
    print "\n==========================="
    print "total time: " + str(t_pro)
    print "===========================\n"
    #sys.exit(status)

