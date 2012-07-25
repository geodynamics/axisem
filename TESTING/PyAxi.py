#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  PyAxi.py
#   Purpose:   Interface for AXISEM in Python
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#-------------------------------------------------------------------

#for debugging: import ipdb; ipdb.set_trace()

#-----------------------------------------------------------------------
#----------------Import required Modules (Python)-----------------------
#-----------------------------------------------------------------------

# Added this line for python 2.5 compatibility
from __future__ import with_statement
import os
import sys
import time
import subprocess
import ConfigParser
import shutil
import glob

global obspy_error
obspy_error = 'N'

try:
    import numpy as np
    import matplotlib.pyplot as plt
    from obspy.core import read, Trace, Stream, UTCDateTime
    from obspy.signal import util, cosTaper, detrend
except Exception, error:
    obspy_error = 'Please install obspy (http://obspy.org/)'
    print "\n*************************************"
    print error
    print "*************************************"
    print "No MISC and TEST functionalities"
    print obspy_error
    print "*************************************\n"

########################################################################
############################# Main Program #############################
########################################################################

def PyAxi(**kwargs):
    
    """
    PyAxi is the function dedicated to the main part of the code.
    
    To run this code:
    1. change the "inpython.cfg" file based on what you want
    2. type: python PyAxi.py
    """
    
    # global variables
    global input
    
    # ------------------Read INPUT file (Parameters)--------------------
    read_input_file()
    
    if input['test'] != 'N':
        subprocess.check_call(['cp', \
            os.path.join(input['test_folder'], 'inpython.cfg'), \
            input['inpython_address']])
        subprocess.check_call(['cp', \
            os.path.join(input['test_folder'], 'STATIONS'), \
            os.path.join(input['axi_address'], 'SOLVER', 'STATIONS')])
        read_input_file()
    
    # Defining error message
    output = 0
    output_print = 20*'*' + '\n' + str(output) + '\n' + 20*'*'
    
    
    if os.path.exists(os.path.join(input['axi_address'], 'SOLVER', input['solver_name'])) == True:
        print '--------------------------------------------------------'
        user_raw_input = raw_input('Folder with the same name as ' + input['solver_name'] + \
                ' \n(solver_name in ipython.cfg) exists in:\n' + \
                os.path.join(input['axi_address'], 'SOLVER', input['solver_name']) + \
                '\n\n' + \
                'You could:' + '\n' + \
                'S: Stop the program' + '\n' + \
                'R: Remove the folder and continue the program' + '\n' + \
                'C: continue the program without removing the folder' + '\n\n').upper()
        
        print '--------------------------------------------------------'
        if user_raw_input == 'C':
            print "Continue the program!"
            print '--------------------------------------------------------'
        elif user_raw_input == 'R':
            print "Removing the directory and continue the program!"
            print '--------------------------------------------------------'
            #shutil.rmtree(os.path.join(input['axi_address'], 'SOLVER', input['solver_name']))
            subprocess.check_call(['rm', '-rf', \
                os.path.join(input['axi_address'], 'SOLVER', input['solver_name'])])
        elif user_raw_input == 'S':
            print '------------'
            print 'EXIT AXISEM!'
            print '------------'
            sys.exit()
        
    ##############################################################
    ############################# MESHER #########################
    ##############################################################
    t1_mesher = time.time()
    
    if input['mesher'] != 'N':
        print "\n======"
        print "MESHER"
        print "======\n"
        os.chdir(os.path.join(input['axi_address'], 'MESHER'))
        
        # Create Mesher Makefile (required just once!)
        if input['mesher_makefile'] != 'N':
            print "========================="
            print "Create Mesher Makefile"
            print "=========================\n"
            output = subprocess.check_call('./makemake.pl')
            if output != 0: print output_print
            
        # Change the input files + make clean; make
        if input['mesher_make'] != 'N':
            
            # Change the mpi_compiler based on the local machine
            if input['mpi_compiler'] != 'N':
                print "==============="
                print "Change F90 flag"
                print "==============="

                makefile_open = open('./Makefile', 'r')
                makefile_read = makefile_open.readlines()
                search = 'F90 = '
                for i in range(0, len(makefile_read)):
                    if makefile_read[i].find(search) != -1:
                        makefile_read[i] = 'F90 = ' + input['mpi_compiler'] + '\n'
                        num = i
                makefile_open.close()
                makefile_open = open('./Makefile', 'w')

                for i in range(0, len(makefile_read)):
                    makefile_open.write(makefile_read[i])

                makefile_open.close()
                print makefile_read[num]
            
            print "==================="
            print "Change inparam_mesh"
            print "==================="

            inparam_mesh_open = open('./inparam_mesh', 'r')
            inparam_mesh_read = inparam_mesh_open.readlines()

            inparam_mesh_read[1] = input['model'] + \
                "  \t    Background model: 'prem','prem_solid' etc (SEE BELOW)\n"
            inparam_mesh_read[2] = input['period'] + \
                "        \t    DOMINANT period [s]\n"
            inparam_mesh_read[3] = input['no_proc'] + \
                "              \t    Number of processors to be used\n"

            inparam_mesh_open.close()
            inparam_mesh_open = open('./inparam_mesh', 'w')

            for i in range(0, len(inparam_mesh_read)):
                inparam_mesh_open.write(inparam_mesh_read[i])

            inparam_mesh_open.close()
            print inparam_mesh_read[1] + inparam_mesh_read[2] + inparam_mesh_read[3]

            print "================"
            print "make clean; make"
            print "================"

            output = subprocess.check_call(['make', 'clean'])
            if output != 0: print output_print
            
            output = subprocess.check_call(['make'])
            if output != 0: print output_print
                
                
            if os.path.isfile('xmesh'):
                print "\n================="
                print "xmesh --- Created"
                print "================="
            
            print "RUN submit.csh"
            print "================="
            output = subprocess.check_call(['./submit.csh'])
            if output != 0: print output_print
        
        # move the mesh to the SOLVER folder
        if input['mesher_move'] != 'N':
            print "==================="
            print "Move mesh files to:"
            print os.path.join('..', 'SOLVER', 'MESHES', input['mesh_name'])
            print "==================="
            if os.path.isfile(os.path.join('..', 'SOLVER', 'MESHES', input['mesh_name'])):
                subprocess.check_call(['rm', '-rf', \
                    os.path.join('..', 'SOLVER', 'MESHES', input['mesh_name'])])
            
            output = subprocess.check_call(['./movemesh', input['mesh_name']])
            if output != 0: print output_print
            
    t2_mesher = time.time()
    ##############################################################
    ############################# SOLVER #########################
    ##############################################################
    t1_solver = time.time()
    
    if input['solver'] != 'N':
        print "\n======"
        print "SOLVER"
        print "======\n"
        os.chdir(os.path.join(input['axi_address'], 'SOLVER'))
        
        # Copy the mesh_params.h in the main directory
        print "============================"
        print "Defines the mesh as a header"
        print "============================\n"
        if input['solver_cp'] != 'N':
            output = subprocess.check_call(['cp', '-f', os.path.join('MESHES', \
                        input['mesh_name'], 'mesh_params.h'), '.'])
            if output != 0: print output_print
        
        # Create SOLVER Makefile
        if input['solver_makefile'] != 'N':
            print "========================="
            print "Create Solver Makefile"
            print "=========================\n"
            if input['netCDF'] == 'N':
                output = subprocess.check_call('./makemake.csh')
                if output != 0: print output_print
            elif input['netCDF'] != 'N':
                output = subprocess.check_call(['./makemake.csh', '-netcdf'])
                if output != 0: print output_print
        # Change the input files + make clean; make
        if input['solver_make'] != 'N':
            
            if input['mpi_compiler'] != 'N':
                print "\n==============="
                print "Change F90 flag"
                print "==============="

                makefile_open = open('./Makefile', 'r')
                makefile_read = makefile_open.readlines()
                search = 'F90 = '
                for i in range(0, len(makefile_read)):
                    if makefile_read[i].find(search) != -1:
                        makefile_read[i] = 'F90 = ' + input['mpi_compiler'] + '\n'
                        num = i
                makefile_open.close()
                makefile_open = open('./Makefile', 'w')

                for i in range(0, len(makefile_read)):
                    makefile_open.write(makefile_read[i])

                makefile_open.close()
                print makefile_read[num]
            
            if input['netCDF'] != 'N':
                print "\n==================="
                print "Change netCDF flags"
                print "==================="
                
                makefile_open = open('./Makefile', 'r')
                makefile_read = makefile_open.readlines()
                search = 'LIBS = '
                for i in range(0, len(makefile_read)):
                    if makefile_read[i].find(search) != -1:
                        makefile_read[i] = 'LIBS = -lm -lnetcdf -lnetcdff' + '\n'
                        num_1 = i
                search = 'INCLUDE ='
                for i in range(0, len(makefile_read)):
                    if makefile_read[i].find(search) != -1:
                        makefile_read[i] = 'INCLUDE = ' + input['netCDF_include'] + '\n'
                        num_2 = i
                makefile_open.close()
                makefile_open = open('./Makefile', 'w')

                for i in range(0, len(makefile_read)):
                    makefile_open.write(makefile_read[i])

                makefile_open.close()
                print makefile_read[num_1]
                print makefile_read[num_2]
            
            print "\n=============="
            print "Change inparam"
            print "=============="
            
            if not os.path.isfile('inparam'):
                subprocess.check_call(['cp', 'inparam.TEMPLATE', 'inparam'])
            
            inparam_solver_open = open('./inparam', 'r')
            inparam_solver_read = inparam_solver_open.readlines()
            
            inparam_solver_read[1] = input['no_simu'] + \
                '               number of simulations. 1: single Mij/f_i; 2: forces; 4: moment tensor \n'
            inparam_solver_read[4] = input['seis_length'] + \
                '           seismogram length [s]\n'
            inparam_solver_read[5] = input['time_step'] + \
                "            time step [s]. Put to 0.0 to use mesher's suggestion (mesh_params.h)\n"
            inparam_solver_read[10] = input['source_type'] + \
                "    source file type: 'sourceparams','cmtsolut'\n"
            inparam_solver_read[11] = input['receiver_type'] + \
                "        receiver file type: 'colatlon','stations'\n"
            
            if input['netCDF'] == 'N':
                inparam_solver_read[34] = 'binary' + \
                    '          Output format for seismograms and wavefields: binary, netcdf\n'
            elif input['netCDF'] != 'N':
                inparam_solver_read[34] = 'netcdf' + \
                    '          Output format for seismograms and wavefields: binary, netcdf\n'
            
            inparam_solver_open.close()
            inparam_solver_open = open('./inparam', 'w')

            for i in range(0, len(inparam_solver_read)):
                inparam_solver_open.write(inparam_solver_read[i])

            inparam_solver_open.close()
            print inparam_solver_read[1] + inparam_solver_read[4] + \
                        inparam_solver_read[10] + inparam_solver_read[11]
            
            if input['source_type'] == 'sourceparams':
                print "\n========================"
                print "Change the Source params"
                print "========================"
                
                source_open = open('./sourceparams.dat', 'r')
                source_read = source_open.readlines()
                source_read[0] = input['src_Mzz'] + ' ' + input['src_Mxx'] + ' ' + \
                                input['src_Myy'] + ' ' + input['src_Mxz'] + ' ' + \
                                input['src_Myz'] + ' ' + input['src_Mxy'] + ' ' + \
                "  moment tensor (Mzz Mxx Myy Mxz Myz Mxy) [Nm]\n"
                source_read[1] = input['sourceparams_type'] + \
                    "       excitation type: 'monopole', 'dipole', 'quadpole'\n"
                source_read[2] = input['sourceparams_MDQ'] + \
                    "            'explosion','mxx_p_myy','mzz','vertforce' (MONOPOLE) \n"
                source_read[5] = input['source_dp'] + \
                    '              source depth [km]\n'
                source_read[6] = input['source_colat'] + \
                    '              source colatitude [degrees]  \n'
                source_read[7] = input['source_lon'] + \
                    '              source longitude [ldegrees]\n'
                source_read[8] = input['source_stf'] + \
                    "        source time function: \n"

                source_open.close()
                source_open = open('./sourceparams.dat', 'w')

                for i in range(0, len(source_read)):
                    source_open.write(source_read[i])

                source_open.close()
                print source_read[1] + source_read[2] + \
                            source_read[5] + source_read[6] + \
                            source_read[7] + source_read[8]
                
            elif input['source_type'] == 'cmtsolut':
                print "\n========================"
                print "Change the Source params"
                print "========================"
                
                if not os.path.isfile('CMTSOLUTION'):
                    subprocess.check_call(['cp', 'CMTSOLUTION.TEMPLATE', 'CMTSOLUTION'])
                
                source_open = open('./CMTSOLUTION', 'r')
                source_read = source_open.readlines()
                
                source_read[4] = 'latitude:      ' + input['cmt_lat'] + '\n'
                source_read[5] = 'longitude:     ' + input['cmt_lon'] + '\n'
                source_read[6] = 'depth:         ' + input['cmt_dp'] + '\n'
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
                print source_read[4] + source_read[5] + \
                            source_read[6] + source_read[7] + \
                            source_read[8] + source_read[9] + \
                            source_read[10] + source_read[11] + \
                            source_read[12]
                
            if input['receiver_type'] == 'colatlon':
                print "\n=========================="
                print "Change the Receiver params"
                print "=========================="
                
                if not os.path.isfile('receivers.dat'):
                    print "\n'receivers.dat' file does not exist!"
                else:
                    receiver_open = open('./receivers.dat', 'r')
                    receiver_read = receiver_open.readlines()
                    print "Number of receivers: " + receiver_read[0] + \
                            "Number of lines    : " + str(len(receiver_read)-1)
                    if int(receiver_read[0]) != (len(receiver_read) - 1):
                        "\nNumber of stations is not entered correctly!"
                        
            elif input['receiver_type'] == 'stations':
                print "\n=========================="
                print "Change the Receiver params"
                print "=========================="
                
                if not os.path.isfile('STATIONS'):
                    subprocess.check_call(['cp', 'STATIONS.TEMPLATE', 'STATIONS'])
                receiver_open = open('./STATIONS', 'r')
                receiver_read = receiver_open.readlines()
                print "\nNumber of receivers: " + str(len(receiver_read))
                
            print "\n================"
            print "make clean; make"
            print "================"

            output = subprocess.check_call(['make', 'clean'])
            if output != 0: print output_print
                
            output = subprocess.check_call(['make'])
            if output != 0: print output_print
                
                
            if os.path.isfile('xsem'):
                print "\n================="
                print "xsem --- Created"
                print "================="
            
            print "RUN submit.csh"
            print "================="
            
            output = subprocess.check_call(['./submit.csh', input['solver_name']])
            if output != 0: print output_print
            
            # copy the input file to the folder (for reference!)
            print "=================="
            print "Copy the input in:\n" + os.path.join(input['axi_address'], 'SOLVER', input['solver_name'])
            print "=================="
            
            subprocess.check_call(['cp', \
                        input['inpython_address'], \
                        os.path.join(input['axi_address'], 'SOLVER', input['solver_name'], "inpython.cfg")])
            
            # check the OUTPUT file and inform the user whenever the program is done!
            print "================"
            print "Check the OUTPUT"
            print "================"
            
            os.chdir(os.path.join(input['axi_address'], 'SOLVER', input['solver_name']))
            
            test = -1
            test_1 = -1
            test_2 = -1
            test_3 = -1
            test_4 = -1
            
            time.sleep(20)
            print_output = "Just after 20 seconds!"
            
            if input['source_type'] == 'sourceparams':
                while (test == -1):
                    time.sleep(2)
                    output_file_open =  open('OUTPUT_' + input['solver_name'], 'r')
                    output_file_read = output_file_open.readlines()
                    test = output_file_read[-1].find('PROGRAM axisem FINISHED')
                    for k in range(0, int(input['no_proc'])):
                        if output_file_read[-1-k].find('PROGRAM axisem FINISHED') == -1:
                            test = -1
                            print_output = output_file_read[-1-k].split('\n')[0]
                    print print_output
                    
                    print '--------------------------------'
                    
                    
            elif input['source_type'] == 'cmtsolut':
                while (test_1 == -1 or test_2 == -1 or test_3 == -1 or test_4 == -1):
                        time.sleep(2)
                        output_file_open =  open(os.path.join('MXX_P_MYY', 'OUTPUT_MXX_P_MYY'), 'r')
                        output_file_read = output_file_open.readlines()
                        test_1 = output_file_read[-1].find('PROGRAM axisem FINISHED')
                        for k in range(0, int(input['no_proc'])):
                            if output_file_read[-1-k].find('PROGRAM axisem FINISHED') == -1:
                                test_1 = -1
                                print_output = output_file_read[-1-k].split('\n')[0]
                        print 'MXX_P_MYY:     ' + print_output
                        
                        output_file_open =  open(os.path.join('MXY_MXX_M_MYY', 'OUTPUT_MXY_MXX_M_MYY'), 'r')
                        output_file_read = output_file_open.readlines()
                        test_2 = output_file_read[-1].find('PROGRAM axisem FINISHED')
                        for k in range(0, int(input['no_proc'])):
                            if output_file_read[-1-k].find('PROGRAM axisem FINISHED') == -1:
                                test_2 = -1
                                print_output = output_file_read[-1-k].split('\n')[0]
                        print 'MXY_MXX_M_MYY: ' + print_output
                        
                        output_file_open =  open(os.path.join('MXZ_MYZ', 'OUTPUT_MXZ_MYZ'), 'r')
                        output_file_read = output_file_open.readlines()
                        test_3 = output_file_read[-1].find('PROGRAM axisem FINISHED')
                        for k in range(0, int(input['no_proc'])):
                            if output_file_read[-1-k].find('PROGRAM axisem FINISHED') == -1:
                                test_3 = -1
                                print_output = output_file_read[-1-k].split('\n')[0]
                        print 'MXZ_MYZ:       ' + print_output
                        
                        output_file_open =  open(os.path.join('MZZ', 'OUTPUT_MZZ'), 'r')
                        output_file_read = output_file_open.readlines()
                        test_4 = output_file_read[-1].find('PROGRAM axisem FINISHED')
                        for k in range(0, int(input['no_proc'])):
                            if output_file_read[-1-k].find('PROGRAM axisem FINISHED') == -1:
                                test_4 = -1
                                print_output = output_file_read[-1-k].split('\n')[0]
                        print 'MZZ:           ' + print_output
                        
                        print '--------------------------------'
                        
    t2_solver = time.time()
    ##############################################################
    ######################## Post-Processing #####################
    ##############################################################
    t1_post = time.time()
    
    if input['post_processing'] != 'N':
        print "\n==============="
        print "Post Processing"
        print "===============\n"
        
        if input['source_type'] == 'sourceparams':
            os.chdir(os.path.join(input['axi_address'], 'SOLVER', input['solver_name']))
            post_process_open = open('./param_post_processing', 'r')
            post_process_read = post_process_open.readlines()
            
            post_process_read = edit_param_post_processing(post_process_read = post_process_read)
            
            post_process_open.close()
            post_process_open = open('./param_post_processing', 'w')
            
            for i in range(0, len(post_process_read)):
                post_process_open.write(post_process_read[i])
                print post_process_read[i].split('\n')[0]
            print 2*"======================================"
            post_process_open.close()
            
            output = subprocess.check_call(['./post_processing.csh'])
            if output != 0: print output_print
            

        elif input['source_type'] == 'cmtsolut':
            for add in ['MXX_P_MYY', 'MXY_MXX_M_MYY', 'MXZ_MYZ', 'MZZ']:
                os.chdir(os.path.join(input['axi_address'], 'SOLVER', input['solver_name'], add))
                post_process_open = open('./param_post_processing', 'r')
                post_process_read = post_process_open.readlines()
                
                post_process_read = edit_param_post_processing(post_process_read = post_process_read)
                
                post_process_open.close()
                post_process_open = open('./param_post_processing', 'w')
                
                print '\n' + add + ':'
                for i in range(0, len(post_process_read)):
                    post_process_open.write(post_process_read[i])
                    print post_process_read[i].split('\n')[0]
                print 2*"======================================"
                post_process_open.close()
            
            os.chdir(os.path.join(input['axi_address'], 'SOLVER', input['solver_name']))
            output = subprocess.check_call(['./post_processing.csh'])
            if output != 0: print output_print
        
    t2_post = time.time()
    ##############################################################
    ############################# MISC ###########################
    ##############################################################
    t1_misc = time.time()
    
    if input['mseed'] != 'N':
        print "\n===================="
        print "Creating MSEED files"
        print "====================\n"
        
        path = os.path.join(input['axi_address'], 'SOLVER', input['solver_name'], \
                                'Data_Postprocessing', 'SEISMOGRAMS')
                
        axisem2mseed(path = path)
        
    t2_misc = time.time()
    ##############################################################
    ############################# TEST ###########################
    ##############################################################
    t1_test = time.time()
    
    if input['test'] != 'N' or input['plot'] != 'N':        
        
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
        
        subprocess.check_call(['cp', \
                os.path.join(input['axi_address'], 'SOLVER', input['solver_name'], \
                                'Data_Postprocessing', 'SEISMOGRAMS', 'seismograms.mseed'), \
                os.path.join(folder_new, 'axisem.mseed')])
        
        """
        if not os.path.isfile(os.path.join(folder, 'axisem.mseed')):
            subprocess.check_call(['cp', \
                os.path.join(input['axi_address'], 'SOLVER', input['solver_name'], \
                                'Data_Postprocessing', 'SEISMOGRAMS', 'seismograms.mseed'), \
                os.path.join(folder, 'axisem.mseed')])
        """
        
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
            sg.filter('lowpass', freq=fmax, corners=2)
            sg.filter('lowpass', freq=fmax, corners=2)
            sg.filter('lowpass', freq=fmax, corners=2)
            sg.filter('highpass', freq=fmin, corners=2)
            sg.filter('highpass', freq=fmin, corners=2)

        nstat = eval(input['test_nstat'])
        stats = np.arange(nstat) + 1

        labels = ['Axisem-ref', 'Yspec', 'Axisem-new']
        colors = ['k', 'r', 'b']
        linestyles = ['--', '-', '-']

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
        
        for chan in chans:
            plt.figure()
            maxi = 0.

            # find global maxium
            for stat in stats:
                str_stat = '*%02d' % stat
                tr = sgs[0].select(station=str_stat, channel='*'+chan)[0]
                dat1 = tr.data
                maxi = max(maxi, np.abs(dat1).max())

            print maxi


            for n in stats:
                stat = '*%02d' % (n,)
                print stat
                for i, sg in enumerate(sgs):
                    dat = sg.select(station=stat, channel='*'+chan)[0].data
                    if i == 0:
                        maxl = dat.max()
                    dat = dat / maxi
                    #dat = dat / maxl / 2.
                    
                    if n == 1:
                        plt.plot(t[i], dat + n, colors[i], label=labels[i], ls=linestyles[i])
                    else:
                        plt.plot(t[i], dat + n, colors[i], label='_nolegend_', ls=linestyles[i])


            plt.xlabel('time / seconds')
            plt.xlim(0, 1800)
            plt.ylim(0, nstat + 1)

            plt.legend()

            fig = plt.gcf()
            fig.set_size_inches((16,12))

            plt.suptitle(chan)
            plt.subplots_adjust(hspace=0.3, wspace=0.2, left=0.08, right=0.95, top=0.93, bottom=0.07)

            #plt.savefig('record_section_%s.pdf' % (chan))
        plt.show()
    
    t2_test = time.time()
    
    print "============================================================"
    print "AXISEM statistics:"
    print "Time for MESHER: " + str(t2_mesher - t1_mesher)
    print "Time for SOLVER: " + str(t2_solver - t1_solver)
    print "Time for POST  : " + str(t2_post - t1_post)
    print "Time for MISC  : " + str(t2_misc - t1_misc)
    print "Time for TEST  : " + str(t2_test - t1_test)
    print "============================================================"
       
    
###################### read_input_file #################################

def read_input_file():  
    
    """
    Read inputs from inpython.cfg file.
    """
    
    global input, obspy_error
    
    config = ConfigParser.RawConfigParser()
    config.read(os.path.join(os.getcwd(), 'inpython.cfg'))
    
    input = {}
    input['inpython_address'] = os.path.join(os.getcwd(), 'inpython.cfg')
    
    input['axi_address'] = config.get('general', 'axi_address')
    if not os.path.isabs(input['axi_address']):
        input['axi_address'] = os.path.join(os.getcwd(), input['axi_address'])
    input['mesh_name'] = config.get('general', 'mesh_name')
    input['solver_name'] = config.get('general', 'solver_name')
    
    #input['all_steps'] = config.get('general', 'all_steps')
    input['new_mesh'] = config.get('general', 'new_mesh')
    input['post_processing'] = config.get('general', 'post_processing')
    input['test'] = config.get('general', 'test')
    
    input['test_folder'] = config.get('general', 'test_folder')
    
    input['plot'] = config.get('general', 'plot')
    
    input['mesher'] = config.get('general', 'mesher')
    input['solver'] = config.get('general', 'solver')
    
    input['mesher_makefile'] = config.get('general', 'mesher_makefile')
    input['mesher_make'] = config.get('general', 'mesher_make')
    input['mesher_move'] = config.get('general', 'mesher_move')
    
    input['solver_cp'] = config.get('general', 'solver_cp')
    input['solver_makefile'] = config.get('general', 'solver_makefile')
    input['solver_make'] = config.get('general', 'solver_make')
    
    if input['mesher_make'] != 'N':
        input['mesher_move'] = 'Y'
        input['solver_cp'] = 'Y'
    
    
    input['mpi_compiler'] = config.get('mpi_netCDF', 'mpi_compiler')
    input['netCDF'] = config.get('mpi_netCDF', 'netCDF')
    input['netCDF_include'] = config.get('mpi_netCDF', 'netCDF_include')

    
    input['model'] = config.get('mesher', 'model')
    input['period'] = config.get('mesher', 'period')
    input['no_proc'] = config.get('mesher', 'no_proc')
    
    
    input['no_simu'] = config.get('solver', 'no_simu')
    input['seis_length'] = config.get('solver', 'seis_length')
    input['time_step'] = config.get('solver', 'time_step')
    input['source_type'] = config.get('solver', 'source_type')
    input['receiver_type'] = config.get('solver', 'receiver_type')
    
    input['sourceparams_type'] = config.get('solver', 'sourceparams_type')
    input['sourceparams_MDQ'] = config.get('solver', 'sourceparams_MDQ')
    
    input['src_Mzz'] = config.get('solver', 'src_Mzz')
    input['src_Mxx'] = config.get('solver', 'src_Mxx')
    input['src_Myy'] = config.get('solver', 'src_Myy')
    input['src_Mxz'] = config.get('solver', 'src_Mxz')
    input['src_Myz'] = config.get('solver', 'src_Myz')
    input['src_Mxy'] = config.get('solver', 'src_Mxy')
    
    input['source_dp'] = config.get('solver', 'source_dp')
    input['source_colat'] = config.get('solver', 'source_colat')
    input['source_lon'] = config.get('solver', 'source_lon')
    input['source_stf'] = config.get('solver', 'source_stf')
    
    input['cmt_lat'] = config.get('solver', 'cmt_lat')
    input['cmt_lon'] = config.get('solver', 'cmt_lon')
    input['cmt_dp'] = config.get('solver', 'cmt_dp')
    input['cmt_Mrr'] = config.get('solver', 'cmt_Mrr')
    input['cmt_Mtt'] = config.get('solver', 'cmt_Mtt')
    input['cmt_Mpp'] = config.get('solver', 'cmt_Mpp')
    input['cmt_Mrt'] = config.get('solver', 'cmt_Mrt')
    input['cmt_Mrp'] = config.get('solver', 'cmt_Mrp')
    input['cmt_Mtp'] = config.get('solver', 'cmt_Mtp')
    
    input['post_components'] = config.get('post_processing', 'post_components')
    input['post_conv_period'] = config.get('post_processing', 'post_conv_period')
    input['post_rotate'] = config.get('post_processing', 'post_rotate')
    input['post_full_Mij'] = config.get('post_processing', 'post_full_Mij')
    input['post_Mrr'] = config.get('post_processing', 'post_Mrr')
    input['post_Mtt'] = config.get('post_processing', 'post_Mtt')
    input['post_Mpp'] = config.get('post_processing', 'post_Mpp')
    input['post_Mrt'] = config.get('post_processing', 'post_Mrt')
    input['post_Mrp'] = config.get('post_processing', 'post_Mrp')
    input['post_Mtp'] = config.get('post_processing', 'post_Mtp')
    input['post_STF'] = config.get('post_processing', 'post_STF')
    input['post_Scolat'] = config.get('post_processing', 'post_Scolat')
    input['post_Slon'] = config.get('post_processing', 'post_Slon')
    input['post_snap'] = config.get('post_processing', 'post_snap')
    input['post_dv'] = config.get('post_processing', 'post_dv')
    input['post_path'] = config.get('post_processing', 'post_path')
    input['post_negative'] = config.get('post_processing', 'post_negative')
    
    input['mseed'] = config.get('MISC', 'mseed')
    
    input['test_chans'] = config.get('test', 'chans')
    input['test_fmin'] = config.get('test', 'fmin')
    input['test_fmax'] = config.get('test', 'fmax')
    input['test_half'] = config.get('test', 'halfduration')
    input['test_nstat'] = config.get('test', 'nstat')
    
    """
    if input['all_steps'] != 'N':
        input['mesher'] = 'Y'
        input['solver'] = 'Y'
        input['mesher_makefile'] = 'Y'
        input['mesher_make'] = 'Y'
        input['mesher_move'] = 'Y'
        input['solver_cp'] = 'Y'
        input['solver_makefile'] = 'Y'
        input['solver_make'] = 'Y'
    """
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
    
    if input['plot'] != 'N':
        input['mesher'] = 'N'
        input['solver'] = 'N'
        input['solver_cp'] = 'N'
        input['solver_makefile'] = 'N'
        input['solver_make'] = 'N'
        input['post_processing'] = 'N'
        input['test'] = 'N'
        input['mseed'] = 'N'
    
    if obspy_error != 'N':
        input['test'] = 'N'
        input['mseed'] = 'N'

###################### edit_param_post_processing ######################

def edit_param_post_processing(post_process_read):
    
    """
    edit param_post_processing file
    """
    
    global input

    if input['post_rotate'] != 'N': \
        post_process_read[0] = '                   ' + \
        input['post_rotate'] + \
        '                                 rotate receivers?\n'
    if input['post_components'] != 'N': \
        post_process_read[1] = "                    " + \
        input['post_components'] + \
        "          receiver components: enz,sph,cyl,xyz,src\n"
    if input['post_full_Mij'] != 'N': \
        post_process_read[2] = '                   ' + \
        input['post_full_Mij'] + \
        '                                   sum to full Mij\n'
    if input['post_Mrr'] != 'N': \
        post_process_read[3] = '             ' + \
        input['post_Mrr'] + \
        '                                              Mrr \n'
    if input['post_Mtt'] != 'N': \
        post_process_read[4] = '             ' + \
        input['post_Mtt'] + \
        '                                              Mtt \n'
    if input['post_Mpp'] != 'N': \
        post_process_read[5] = '             ' + \
        input['post_Mpp'] + \
        '                                              Mpp \n'
    if input['post_Mrt'] != 'N': \
        post_process_read[6] = '             ' + \
        input['post_Mrt'] + \
        '                                              Mrt \n'
    if input['post_Mrp'] != 'N': \
        post_process_read[7] = '             ' + \
        input['post_Mrp'] + \
        '                                              Mrp \n'
    if input['post_Mtp'] != 'N': \
        post_process_read[8] = '             ' + \
        input['post_Mtp'] + \
        '                                              Mtp \n'
    if input['post_conv_period'] != 'N': \
        post_process_read[9] = '             ' + \
        input['post_conv_period'] + \
        '             convolve period (0. if not convolved)\n'
    if input['post_STF'] != 'N': \
        post_process_read[10] = "                " + \
        input['post_STF'] + \
        "         source time function type for convolution\n"
    if input['post_Scolat'] != 'N': \
        post_process_read[11] = '             ' + \
        input['post_Scolat'] + \
        '                                 Source colatitude\n'
    if input['post_Slon'] != 'N': \
        post_process_read[12] = '             ' + \
        input['post_Slon'] + \
        '                                  Source longitude\n'
    if input['post_snap'] != 'N': \
        post_process_read[13] = '                        ' + \
        input['post_snap'] + \
        '                                plot global snaps?\n'
    if input['post_dv'] != 'N': \
        post_process_read[14] = '                     ' + \
        input['post_dv'] + \
        '                          disp or velo seismograms\n'
    if input['post_path'] != 'N': \
        post_process_read[15] = "    " + \
        input['post_path'] + \
        "                 Directory for post processed data\n"
    if input['post_negative'] != 'N': \
        post_process_read[16] = '                        ' + \
        input['post_negative'] + \
        '   seismograms at negative time (0 at max. of stf)\n'
    
    return post_process_read

########################## axisem2mseed ################################

def axisem2mseed(path):
    
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

    print fname
    st.write(fname, format='MSEED')

########################## axisem2mseed ################################

def convSTF(st, sigma=30.):

    gauss = lambda (t, s): 1. / (2. * np.pi * s**2.)**.5 \
                           * np.exp(-1*(t**2)/(2*(s**2)))

    df = st[0].stats.sampling_rate
    dt = 1./df

    t = np.linspace(0., sigma * 20., sigma * 20 * df + 1)
    stf = gauss((t-sigma*10, sigma))
    nstf = len(stf)

    for tr in st:
        tr.data *= cosTaper(len(tr.data), p=0.05)
        nfft = util.nextpow2(max(nstf, tr.stats.npts)) * 2
        stff = np.fft.rfft(stf, n=nfft) * dt
        trf = np.fft.rfft(tr, n=nfft) * dt
        tr.data = np.fft.irfft(stff * trf)[sigma*10*df:sigma*10*df+len(tr.data)] * df

    return 1
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
