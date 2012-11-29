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
import os
import sys
import time
import subprocess
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
    from obspy.signal import util, cosTaper, detrend
except Exception, error:
    obspy_error = 'Please install obspy (http://obspy.org/)'
    print "\n======================================"
    print 'Error in importing:'
    print error
    print "------------------"
    print "No MISC and TEST functionalities"
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
       will ask whether it should remove the folder or not!
    """
    
    # global variables
    global input
    try:
        if sys.argv[1] == '--check':
            print 'Check the Basic, Processing and Visualization requirements:\n'
            p = []
            
            print repr('Compiler').rjust(20).replace("'", '') + "    |    " +   repr('Installed').rjust(11).replace("'", '')
            print '  ' + 42*'-'
            
            for compilers in ['gfortran', 'ifort', 'mpif90', 'mpif90.openmpi', \
                                'gnuplot', 'taup', 'paraview', 'matlab', 'googleearth']:
                p.append([compilers, subprocess.Popen(['which', compilers], \
                                    stdout = subprocess.PIPE).communicate()])
            for comp in range(0, len(p)):
                if p[comp][1][0] != '':
                    print repr(p[comp][0]).rjust(20).replace("'", '') + "    |    " +  repr('Y').rjust(7).replace("'", '')
                else:
                    print repr(p[comp][0]).rjust(20).replace("'", '') + "    |    " +  repr('N').rjust(7).replace("'", '')
            
            print '\nSummary:\n'
            if p[0][1][0] != '' or p[1][1][0] != '':
                if p[2][1][0] != '' or p[3][1][0] != '':
                    print 'Basic Functionality...CHECK'
                else:
                    print 'Basic Functionality...ERROR'
            else:
                print 'Basic Functionality...ERROR'
            if p[4][1][0] != '' and p[5][1][0] != '': 
                print 'Processing Req........CHECK'
            else:
                print 'Processing Req........ERROR'
            if p[6][1][0] != '' and p[7][1][0] != '' and p[8][1][0] != '': 
                print 'Visualization Tools...CHECK'
            else:
                print 'Visualization Tools...ERROR'
            print '\n'
            sys.exit(2) 
    except Exception, error:
        print '',
    # ------------------Read INPUT file (Parameters)--------------------
    read_input_file()
    
    if input['test'] != 'N':
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
            subprocess.check_call(['rm', '-rf', \
                os.path.join(input['axi_address'], 'SOLVER', input['solver_name'])])
        elif user_raw_input == 'S':
            print 'EXIT AXISEM!'
            sys.exit()
        
    ##############################################################
    ############################# MESHER #########################
    ##############################################################
    t1_mesher = time.time()
    
    if input['mesher'] != 'N':
        print "\n======"
        print "MESHER"
        print "======"
        os.chdir(os.path.join(input['axi_address'], 'MESHER'))
        
        # Delete previous mesh_params.h, meshdb.dat*, Diags/*
        if input['verbose'] != 'N': 
            print "\n=========================================================================="
            print 'Removing old mesh_params.h*, meshdb.dat*, ./Diags/* and unrolled_loops.f90'
            print "==========================================================================\n"
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
        output = subprocess.check_call(['rm', '-rf', './unrolled_loops.f90'], stdout = stdout_param)
        if output != 0: print output_print
        print 'DONE'

        # Create Mesher Makefile (required just once!)
        if input['mesher_makefile'] != 'N':
            if input['verbose'] != 'N':
                print "========================="
                print "Create Mesher Makefile"
                print "=========================\n"
                stdout_param = None
            else:
                sys.stdout.write('Create Mesher Makefile...')
                sys.stdout.flush()
                stdout_param = subprocess.PIPE

            if input['make_flag'] != 'N':
                output = subprocess.check_call(['./makemake.pl', input['make_flag']], stdout = stdout_param)
                if output != 0: print output_print
            else:
                output = subprocess.check_call(['./makemake.pl'], stdout = stdout_param)
                if output != 0: print output_print
            print 'DONE'
            
        # Change the input files + make clean; make
        if input['mesher_make'] != 'N':
            # Change the mpi_compiler based on the local machine
            if input['mpi_compiler'] != 'N':
                if input['verbose'] != 'N':
                    print "==============="
                    print "Change F90 flag"
                    print "==============="
                else:
                    sys.stdout.write('Change F90 flag...')
                    sys.stdout.flush()
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
                if input['verbose'] != 'N':
                   print makefile_read[num]
                else:
                    print 'DONE'
            
            if input['verbose'] != 'N':
                print "==================="
                print "Change inparam_mesh"
                print "==================="
            else:
                sys.stdout.write('Change inparam_mesh...')
                sys.stdout.flush()
            if os.path.isfile('./inparam_mesh'):
                subprocess.check_call(['rm', './inparam_mesh'])
            subprocess.check_call(['cp', './inparam_mesh.TEMPLATE', './inparam_mesh'])
            
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

            if input['verbose'] != 'N':
                print inparam_mesh_read[1] + inparam_mesh_read[2] + inparam_mesh_read[3]
            else:
                print 'DONE'
            
            if input['verbose'] != 'N':
                print "================"
                print "make clean; make"
                print "================"
                stdout_param = None
            else:
                sys.stdout.write('make clean; make...')
                sys.stdout.flush()
                stdout_param = subprocess.PIPE
            
            output = subprocess.check_call(['make', 'clean'], stdout = stdout_param)
            if output != 0: print output_print
            
            output = subprocess.check_call(['make'], stdout = stdout_param)
            if output != 0: print output_print
            print 'DONE'

            if os.path.isfile('xmesh'):
                if input['verbose'] != 'N':
                    print "\n==============="
                    print "xmesh --- CHECK"
                    print "==============="
                    print "RUN submit.csh"
                    print "================="
                    stdout_param = None
                else:
                    print 'xmesh...CHECK'
                    print 'submit.csh...'
                    stdout_param = subprocess.PIPE
            output = subprocess.check_call(['./submit.csh'], stdout = stdout_param)
            if output != 0: print output_print
            
            # check the xmesh whether it is finished or not!
            process_name= "xmesh"
            check_process="running"
            t_check = 0
            
            print '--------------------------'
            print 'submit.csh is running for:'
            while check_process=="running":
                tmp = os.popen("ps -Af").read()
                if process_name in tmp[:]:
                    print str(t_check) + ', ',
                    time.sleep(1)
                    t_check+=2
                    
                else:
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
                    
                    subprocess.check_call(['rm', '-rf', \
                        os.path.join('..', 'SOLVER', 'MESHES', input['mesh_name'])])
                    
                    print os.path.join('..', 'SOLVER', 'MESHES', input['mesh_name'])
                    print "is removed"
                    print '#######################################################################\n'
                else:
                    sys.stdout.write('Remove the old mesh files...')
                    sys.stdout.flush()
                    subprocess.check_call(['rm', '-rf', \
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
        os.chdir(os.path.join(input['axi_address'], 'SOLVER'))
        
        if input['solver_cp'] != 'N':
            # Copy the mesh_params.h in the main directory
            if input['verbose'] != 'N':
                print "\n======================"
                print "Copy the mesh_params.h"
                print "======================\n"
            else:
                sys.stdout.write('copy the mesh_params.h...')
                sys.stdout.flush()
            output = subprocess.check_call(['cp', '-f', os.path.join('MESHES', \
                        input['mesh_name'], 'mesh_params.h'), '.'])
            if output != 0: print output_print
            print 'DONE'
        
        # Create SOLVER Makefile
        if input['solver_makefile'] != 'N':
            if input['verbose'] != 'N':
                print "======================"
                print "Create Solver Makefile"
                print "======================\n"
                stdout_param = None
            else:
                stdout_param = subprocess.PIPE

            if input['netCDF'] == 'N' and input['make_flag'] == 'N':
                output = subprocess.check_call(['./makemake.pl'], stdout = stdout_param)
                if output != 0: print output_print
            elif input['netCDF'] == 'N' and input['make_flag'] != 'N':
                output = subprocess.check_call(['./makemake.pl', input['make_flag']], stdout = stdout_param)
                if output != 0: print output_print
            elif input['netCDF'] != 'N' and input['make_flag'] == 'N':
                output = subprocess.check_call(['./makemake.pl', '-netcdf'], stdout = stdout_param)
                if output != 0: print output_print
            elif input['netCDF'] != 'N' and input['make_flag'] != 'N':
                output = subprocess.check_call(['./makemake.pl', input['make_flag'], '-netcdf'], stdout = stdout_param)
                if output != 0: print output_print
            print 'Create Solver Makefile...DONE'
       
       # Change the input files + make clean; make
        if input['solver_make'] != 'N':
            
            if input['mpi_compiler'] != 'N':
                if input['verbose'] != 'N':
                    print "\n==============="
                    print "Change F90 flag"
                    print "==============="
                else:
                    sys.stdout.write('Change F90 flag...')
                    sys.stdout.flush()

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
                if input['verbose'] != 'N':
                    print makefile_read[num]
                else:
                    print 'DONE'
            
            if input['netCDF'] != 'N':
                if input['verbose'] != 'N':
                    print "\n==================="
                    print "Change netCDF flags"
                    print "==================="
                else:
                    sys.stdout.write('Change netCDF flags...')
                    sys.stdout.flush()
                
                makefile_open = open('./Makefile', 'r')
                makefile_read = makefile_open.readlines()
                search = '#F90FLAGS = -Dunc'
                for i in range(0, len(makefile_read)):
                    if makefile_read[i].find(search) != -1:
                        makefile_read[i] = 'F90FLAGS = -Dunc' + '\n'
                        num_0 = i
                search = 'LIBS = '
                for i in range(0, len(makefile_read)):
                    if makefile_read[i].find(search) != -1:
                        makefile_read[i] = 'LIBS = ' + input['netCDF_LIBS'] + '\n'
                        num_1 = i
                search = 'INCLUDE ='
                for i in range(0, len(makefile_read)):
                    if makefile_read[i].find(search) != -1:
                        makefile_read[i] = 'INCLUDE = ' + input['netCDF_INCLUDE'] + '\n'
                        num_2 = i
                makefile_open.close()
                makefile_open = open('./Makefile', 'w')

                for i in range(0, len(makefile_read)):
                    makefile_open.write(makefile_read[i])

                makefile_open.close()
                if input['verbose'] != 'N':
                    print makefile_read[num_0]
                    print makefile_read[num_1]
                    print makefile_read[num_2]
                else:
                    print 'DONE'

            if input['verbose'] != 'N':
                print "\n=============="
                print "Change inparam"
                print "=============="
            else:
                sys.stdout.write('Change inparam...')
                sys.stdout.flush()
            
            if os.path.isfile('inparam'):
                subprocess.check_call(['rm', 'inparam'])
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
                "        source file type: 'sourceparams','cmtsolut'\n"
            inparam_solver_read[11] = input['receiver_type'] + \
                "        receiver file type: 'colatlon','stations','database'\n"
            inparam_solver_read[18] = input['save_XDMF'] + \
                "         save XDMF files (high resolution 2D wavefields), more options in inparam_xdmf\n"
            if input['netCDF'] == 'N':
                inparam_solver_read[35] = 'binary' + \
                    '          Output format for seismograms and wavefields: binary, netcdf\n'
            elif input['netCDF'] != 'N':
                inparam_solver_read[35] = 'netcdf' + \
                    '          Output format for seismograms and wavefields: binary, netcdf\n'
            inparam_solver_read[36] = input['force_aniso'] + \
                "         force anisotropic model handling"
            inparam_solver_open.close()
            inparam_solver_open = open('./inparam', 'w')

            for i in range(0, len(inparam_solver_read)):
                inparam_solver_open.write(inparam_solver_read[i])

            inparam_solver_open.close()
            if input['verbose'] != 'N':
                print inparam_solver_read[1] + inparam_solver_read[4] + \
                            inparam_solver_read[10] + inparam_solver_read[11]
            else:
                print 'DONE'

            if input['source_type'] == 'sourceparams':
                if input['verbose'] != 'N':
                    print "\n========================"
                    print "Change the Source params"
                    print "========================"
                else:
                    sys.stdout.write('Change the Source params...')
                    sys.stdout.flush()

                if os.path.isfile('sourceparams.dat'):
                    subprocess.check_call(['rm', 'sourceparams.dat'])
                subprocess.check_call(['cp', 'sourceparams.dat.TEMPLATE', 'sourceparams.dat'])
                
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
                if input['verbose'] != 'N':
                    print source_read[1] + source_read[2] + \
                                source_read[5] + source_read[6] + \
                                source_read[7] + source_read[8]
                else:
                    print 'DONE'

            elif input['source_type'] == 'cmtsolut':
                if input['verbose'] != 'N':
                    print "\n========================"
                    print "Change the Source params"
                    print "========================"
                else:
                    sys.stdout.write('Change the Source params...')
                    sys.stdout.flush()

                if os.path.isfile('CMTSOLUTION'):
                    subprocess.check_call(['rm', 'CMTSOLUTION'])
                subprocess.check_call(['cp', 'CMTSOLUTION.TEMPLATE', 'CMTSOLUTION'])
                
                source_open = open('./CMTSOLUTION', 'r')
                source_read = source_open.readlines()
                
                source_read[0] = input['cmt_STF'] + \
                                " PDE 1994  6  9  0 33 16.40 -13.8300  " + \
                                "-67.5600 637.0 6.9 6.8 NORTHERNBOLIVIA" + '\n'
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
                    print "\n=========================="
                    print "Change the Receiver params"
                    print "=========================="
                else:
                    sys.stdout.write('Change the Receiver params...')
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
                    print "\n=========================="
                    print "Change the Receiver params"
                    print "=========================="
                else:
                    sys.stdout.write('Change the Receiver params...')
                    sys.stdout.flush()
                
                if not os.path.isfile('STATIONS'):
                    subprocess.check_call(['cp', 'STATIONS.TEMPLATE', 'STATIONS'])
                receiver_open = open('./STATIONS', 'r')
                receiver_read = receiver_open.readlines()
                print "Number of receivers: " + str(len(receiver_read))
            
            if input['verbose'] != 'N':
                print "\n================"
                print "make clean; make"
                print "================"
                stdout_param = None
            else:
                sys.stdout.write('make clean; make...')
                sys.stdout.flush()
                stdout_param = subprocess.PIPE

            output = subprocess.check_call(['make', 'clean'], stdout = stdout_param)
            if output != 0: print output_print
                
            output = subprocess.check_call(['make'], stdout = stdout_param)
            if output != 0: print output_print
            
            print 'DONE'
                
            if os.path.isfile('xsem'):
                if input['verbose'] != 'N':
                    print "\n=============="
                    print "xsem --- CHECK"
                    print "==============" 
                    print "RUN submit.csh"
                    print "================="
                    stdout_param = None
                else:
                    print 'xsem...CHECK'
                    sys.stdout.write('submit.csh...')
                    sys.stdout.flush()
                    stdout_param = subprocess.PIPE
            else:
                print 'xsem does not exist!'
                sys.exit(2)
            
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

            subprocess.check_call(['cp', \
                        input['inpython_address'], \
                        os.path.join(input['axi_address'], 'SOLVER', input['solver_name'], "inpython.cfg")])
            print 'DONE'

            # check the OUTPUT file and inform the user whenever the program is done!
            print "\n================"
            print "Check the OUTPUT"
            print "================"
            
            os.chdir(os.path.join(input['axi_address'], 'SOLVER', input['solver_name']))
            
            test = -1
            test_1 = -1
            test_2 = -1
            test_3 = -1
            test_4 = -1
            
            time.sleep(10)
            print_output = "Just after 10 seconds!"
            
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
        if input['verbose'] != 'N':
            print "\n==============="
            print "Post Processing"
            print "===============\n"
            stdout_param = None
        else:
            sys.stdout.write('\nPost Processing...')
            sys.stdout.flush()
            stdout_param = subprocess.PIPE
        
        if os.path.isdir(os.path.join(input['axi_address'], 'SOLVER', \
                            input['solver_name'], 'Data_Postprocessing')):
            subprocess.check_call(['rm', '-rf', os.path.join(input['axi_address'], 'SOLVER', \
                            input['solver_name'], 'Data_Postprocessing')])
        
        if input['source_type'] == 'sourceparams':
            os.chdir(os.path.join(input['axi_address'], 'SOLVER', input['solver_name']))
            post_process_open = open('./param_post_processing', 'r')
            post_process_read = post_process_open.readlines()
            
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
                    if input['verbose'] != 'N': print post_process_read[i].split('\n')[0]
                if input['verbose'] != 'N': print 2*"======================================"
                post_process_open.close()
            
            os.chdir(os.path.join(input['axi_address'], 'SOLVER', input['solver_name']))
            output = subprocess.check_call(['./post_processing.csh'], stdout = stdout_param)
            if output != 0: print output_print
        
        print 'DONE'         
    
    t2_post = time.time()
    ##############################################################
    ############################# MISC ###########################
    ##############################################################
    t1_misc = time.time()
    
    if input['mseed'] != 'N':
        if input['verbose'] != 'N':
            print "\n===================="
            print "Creating MSEED files"
            print "====================\n"
        else:
            sys.stdout.write('Creating MSEED...')
            sys.stdout.flush()
        
        path = os.path.join(input['axi_address'], 'SOLVER', input['solver_name'], \
                                'Data_Postprocessing', 'SEISMOGRAMS')
                
        axisem2mseed(path = path)
        print 'DONE'

    if input['mseed_all'] != 'N':
        if input['verbose'] != 'N':
            print "\n===================="
            print "Creating MSEED files"
            print "====================\n"
        else:
            sys.stdout.write('Creating MSEED...')
            sys.stdout.flush()

        path = os.path.join(input['axi_address'], 'SOLVER', input['solver_name'], \
                                'Data_Postprocessing', 'SEISMOGRAMS')
                
        axisem2mseed_all(path = path)
        print 'DONE'

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

        print t0, dt, npts
        
        i = 1
        t0 = sgs[i][0].stats.starttime
        dt = sgs[i][0].stats.delta
        npts = sgs[i][0].stats.npts
        t.append(np.linspace(0., dt * (npts -1), npts) + (t0 -  UTCDateTime(0)))
        
        print t0, dt, npts
        
        i = 2
        t0 = sgs[i][0].stats.starttime
        dt = sgs[i][0].stats.delta
        npts = sgs[i][0].stats.npts
        t.append(np.linspace(0., dt * (npts -1), npts) + (t0 -  UTCDateTime(0)))
        
        print t0, dt, npts
        
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

            print maxi

            l2misfit = []

            for n in stats:
                stat = '*%02d' % (n,)
                print stat
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

            if input['save_plots'] != 'N':
                recsec.savefig(os.path.join(folder_new, 'record_section_%s.%s' %
                        (chan, input['plot_format'])))

            ax = misfitplot.gca()
            
            ax.semilogy(np.arange(nstat) + 1, np.array(l2misfit) + 1e-12, 'o', label=chan)
            ax.set_xlabel('trace')
            ax.set_ylabel('l2 - misfit to reference data')
            ax.set_ylim(1e-12, 1.)

            ax.axhline(y=1e-8, color='k', ls='--')
            
            if np.max(l2misfit) > 1e-8:
                fwarn = open(os.path.join(folder_new, 'warning.dat'), 'a')
                fwarn.write("maximum l2 norm misfit larger then 1e-8 in chan %s trace %d\n" 
                             % (chan, np.argmax(l2misfit)))
                fwarn.close()

        ax = misfitplot.gca()
        ax.legend()
        if input['save_plots'] != 'N':
            misfitplot.savefig(os.path.join(folder_new, 'l2_misfit.' + input['plot_format']))
        if input['plot'] != 'N':
            plt.show()
    
    t2_test = time.time()
    
    print "\n============================================================"
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
    input = {}
    
    try:
        config.read(os.path.join(sys.argv[1], 'inpython.cfg'))
        input['inpython_address'] = os.path.join(sys.argv[1], 'inpython.cfg')
        print '\n****************************************'
        print 'Read the input file (inpython.cfg) from:'
        print os.path.join(sys.argv[1], 'inpython.cfg')
    except Exception, error:
        config.read(os.path.join(os.getcwd(), 'inpython.cfg'))
        input['inpython_address'] = os.path.join(os.getcwd(), 'inpython.cfg')
        print '\n****************************************'
        print 'Read the input file (inpython.cfg) from:'
        print os.path.join(os.getcwd(), 'inpython.cfg')
    
    if not os.path.isabs(input['inpython_address']):
        input['inpython_address'] = os.path.join(os.getcwd(), \
                                            input['inpython_address'])
    
    input['axi_address'] = config.get('general', 'address')
    if not os.path.isabs(input['axi_address']):
        input['axi_address'] = os.path.join(os.getcwd(), input['axi_address'])
    print '\nWorking Directory:'
    print input['axi_address'] 
    print '****************************************\n'
    
    input['mesh_name'] = config.get('general', 'mesh_name')
    input['solver_name'] = config.get('general', 'solver_name')
    input['verbose'] = config.get('general', 'verbose')
    
    #input['all_steps'] = config.get('general', 'all_steps')
    input['new_mesh'] = config.get('general', 'new_mesh')
    input['post_processing'] = config.get('general', 'post_processing')
    
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
    

    input['make_flag'] = config.get('mpi_netCDF', 'make_flag')
    input['mpi_compiler'] = config.get('mpi_netCDF', 'mpi_compiler')
    input['netCDF'] = config.get('mpi_netCDF', 'netCDF')
    input['netCDF_LIBS'] = config.get('mpi_netCDF', 'netCDF_LIBS')
    input['netCDF_INCLUDE'] = config.get('mpi_netCDF', 'netCDF_INCLUDE')

    
    input['model'] = config.get('mesher', 'model')
    input['period'] = config.get('mesher', 'period')
    input['no_proc'] = config.get('mesher', 'no_proc')
    
    
    input['no_simu'] = config.get('solver', 'no_simu')
    input['seis_length'] = config.get('solver', 'seis_length')
    input['time_step'] = config.get('solver', 'time_step')
    input['source_type'] = config.get('solver', 'source_type')
    input['receiver_type'] = config.get('solver', 'receiver_type')
    input['save_XDMF'] = config.get('solver', 'save_XDMF')
    input['force_aniso'] = config.get('solver', 'force_aniso')

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
    
    input['cmt_STF'] = config.get('solver', 'cmt_STF')
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
    input['mseed_all'] = config.get('MISC', 'mseed_all')
    input['convSTF'] = config.get('MISC', 'convSTF')
    input['halfduration'] = eval(config.get('MISC', 'halfduration'))
    input['filter'] = config.get('MISC', 'filter')
    input['fmin'] = eval(config.get('MISC', 'fmin'))
    input['fmax'] = eval(config.get('MISC', 'fmax'))
    
    input['test'] = config.get('testing', 'test')
    input['test_folder'] = config.get('testing', 'test_folder')
    input['plot'] = config.get('testing', 'plot')
    input['save_plots'] = config.get('testing', 'save_plots')
    input['plot_format'] = config.get('testing', 'plot_format')
    input['test_chans'] = config.get('testing', 'chans')
    input['test_fmin'] = config.get('testing', 'fmin')
    input['test_fmax'] = config.get('testing', 'fmax')
    input['test_half'] = config.get('testing', 'halfduration')
    input['test_nstat'] = config.get('testing', 'nstat')
    
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
        input['plot'] = 'N'  
        input['save_plots'] = 'N'
        
    if input['plot'] != 'N':
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
        input['netCDF'] = 'Y'
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

    if input['post_rotate'] != 'N':
        post_process_read[0] = '                   ' + \
        input['post_rotate'] + \
        '                                 rotate receivers?\n'
    if input['post_components'] != 'N':
        post_process_read[1] = "                    " + \
        input['post_components'] + \
        "          receiver components: enz,sph,cyl,xyz,src\n"
    if input['post_full_Mij'] != 'N':
        post_process_read[2] = '                   ' + \
        input['post_full_Mij'] + \
        '                                   sum to full Mij\n'
    if input['post_Mrr'] != 'N':
        post_process_read[3] = '             ' + \
        input['post_Mrr'] + \
        '                                              Mrr \n'
    if input['post_Mtt'] != 'N':
        post_process_read[4] = '             ' + \
        input['post_Mtt'] + \
        '                                              Mtt \n'
    if input['post_Mpp'] != 'N':
        post_process_read[5] = '             ' + \
        input['post_Mpp'] + \
        '                                              Mpp \n'
    if input['post_Mrt'] != 'N':
        post_process_read[6] = '             ' + \
        input['post_Mrt'] + \
        '                                              Mrt \n'
    if input['post_Mrp'] != 'N':
        post_process_read[7] = '             ' + \
        input['post_Mrp'] + \
        '                                              Mrp \n'
    if input['post_Mtp'] != 'N':
        post_process_read[8] = '             ' + \
        input['post_Mtp'] + \
        '                                              Mtp \n'
    if input['post_conv_period'] != 'N':
        if input['post_conv_period'] == '0.':
            post_process_read[9] = '                       ' + \
            '0.             convolve period (0. if not convolved)\n'
        else:
            post_process_read[9] = '             ' + \
            input['post_conv_period'] + \
            '             convolve period (0. if not convolved)\n'
    if input['post_STF'] != 'N':
        post_process_read[10] = "                " + \
        input['post_STF'] + \
        "         source time function type for convolution\n"
    if input['post_Scolat'] != 'N':
        post_process_read[11] = '             ' + \
        input['post_Scolat'] + \
        '                                 Source colatitude\n'
    if input['post_Slon'] != 'N':
        post_process_read[12] = '             ' + \
        input['post_Slon'] + \
        '                                  Source longitude\n'
    if input['post_snap'] != 'N':
        post_process_read[13] = '                        ' + \
        input['post_snap'] + \
        '                                plot global snaps?\n'
    if input['post_dv'] != 'N':
        post_process_read[14] = '                     ' + \
        input['post_dv'] + \
        '                          disp or velo seismograms\n'
    if input['post_path'] != 'N':
        post_process_read[15] = "    " + \
        input['post_path'] + \
        "                 Directory for post processed data\n"
    if input['post_negative'] != 'N':
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
        if chan == 'E':
            chan = 'BHE'
        elif chan == 'N':
            chan = 'BHN'
        elif chan == 'Z':
            chan = 'BHZ'
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

    print fname
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
