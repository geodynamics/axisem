#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Added this line for python 2.5 compatibility
from __future__ import with_statement
import os
import sys
import time
import subprocess
import ConfigParser

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
    global test_param
    # ------------------Read INPUT file (Parameters)--------------------
    read_input_file()

    ##############################################################
    ############################# MESHER #########################
    ##############################################################

    if test_param['mesher'] != 'N':
        # Change the input files + make clean; make
        if test_param['mesher_make'] != 'N':
            inparam_mesh_input = []
            inparam_mesh_input.append('BACKGROUND_MODEL     %s\n' %(test_param['mesher_bg_model']))
            inparam_mesh_input.append('EXT_MODEL     %s\n' %(test_param['mesher_ext_model']))
            inparam_mesh_input.append('DOMINANT_PERIOD     %s\n' %(test_param['mesher_dominant_period']))
            inparam_mesh_input.append('NTHETA_SLICES     %s\n' %(test_param['mesher_ntheta']))
            inparam_mesh_input.append('NRADIAL_SLICES     %s\n' %(test_param['mesher_nradial']))
            inparam_mesh_input.append('WRITE_VTK     %s\n' %(test_param['mesher_write_vtk']))
            inparam_mesh_input.append('COARSENING_LAYERS     %s\n' %(test_param['mesher_coarsening_layers']))
            inparam_mesh_input.append('IC_SHEAR_WAVE     %s\n' %(test_param['mesher_ic_shear_wave']))
            inparam_mesh_input.append('NPOL     %s\n' %(test_param['mesher_npol']))
            inparam_mesh_input.append('EL_PER_LAMBDA     %s\n' %(test_param['mesher_el_per_lambda']))
            inparam_mesh_input.append('COURANT_NR     %s\n' %(test_param['mesher_courant_nr']))
            inparam_mesh_input.append('RADIUS     %s\n' %(test_param['mesher_radius']))
            inparam_mesh_input.append('SAVE_MESH     %s\n' %(test_param['mesher_save_mesh']))
            inparam_mesh_input.append('VERBOSE     %s\n' %(test_param['mesher_verbose']))
            inparam_mesh_input.append('NCPU      4')

            inparam_mesh_open = open('inparam_mesh', 'w')
            for i in range(0, len(inparam_mesh_input)):
                inparam_mesh_open.write(inparam_mesh_input[i])
            inparam_mesh_open.close()

    ##############################################################
    ############################# SOLVER #########################
    ##############################################################

    if test_param['solver'] != 'N':
        # Change the input files + make clean; make
        if test_param['solver_make'] != 'N':
            inparam_basic_input = []
            inparam_basic_input.append('SIMULATION_TYPE     %s\n' %(test_param['solver_sim_type']))
            inparam_basic_input.append('SEISMOGRAM_LENGTH     %s\n' %(test_param['solver_seis_length']))
            inparam_basic_input.append('RECFILE_TYPE     %s\n' %(test_param['solver_recfile_type']))
            inparam_basic_input.append('MESHNAME     %s\n' %(test_param['mesh_name']))
            inparam_basic_input.append('LAT_HETEROGENEITY     %s\n' %(test_param['solver_lat_heterogeneity']))
            inparam_basic_input.append('ATTENUATION     %s\n' %(test_param['solver_attenuation']))
            inparam_basic_input.append('SAVE_SNAPSHOTS     %s\n' %(test_param['solver_save_snapshots']))
            inparam_basic_input.append('VERBOSITY     %s\n' %(test_param['solver_verbosity']))

            inparam_solver_open = open('inparam_basic', 'w')
            inparam_basic_contents = "".join(inparam_basic_input)
            inparam_solver_open.write(inparam_basic_contents)
            inparam_solver_open.close()

            inparam_advanced_input = []
            inparam_advanced_input.append('SAMPLING_PERIOD         %s\n' %(test_param['solver_sampling_rate']))
            inparam_advanced_input.append('TIME_STEP               %s\n' %(test_param['solver_time_step']))
            inparam_advanced_input.append('SOURCE_PERIOD           %s\n' %(test_param['solver_source_period']))
            inparam_advanced_input.append('SOURCE_FUNCTION         %s\n' %(test_param['source_stf']))
            inparam_advanced_input.append('TIME_SCHEME             %s\n' %(test_param['solver_time_scheme']))
            inparam_advanced_input.append('DATA_DIR               "%s"\n' %(test_param['solver_data_dir']))
            inparam_advanced_input.append('INFO_DIR               "%s"\n' %(test_param['solver_info_dir']))
            inparam_advanced_input.append('DIAGNOSTIC_FILE_OUTPUT "%s"\n' %(test_param['solver_diag_file_output']))
            inparam_advanced_input.append('MESH_TEST               %s\n' %(test_param['solver_mesh_test']))
            inparam_advanced_input.append('DEFLATE_LEVEL           %s\n' %(test_param['solver_deflate_level']))
            inparam_advanced_input.append('SNAPSHOT_DT             %s\n' %(test_param['solver_snapshot_dt']))
            inparam_advanced_input.append('SNAPSHOTS_FORMAT        %s\n' %(test_param['solver_snapshots_format']))
            inparam_advanced_input.append('USE_NETCDF              %s\n' %(test_param['make_axisem_use_netcdf']))
            inparam_advanced_input.append('KERNEL_WAVEFIELDS       %s\n' %(test_param['solver_kernel_wavefields']))
            inparam_advanced_input.append('KERNEL_SPP              %s\n' %(test_param['solver_kernel_spp']))
            inparam_advanced_input.append('KERNEL_SOURCE           %s\n' %(test_param['solver_kernel_source']))
            inparam_advanced_input.append('KERNEL_IBEG             %s\n' %(test_param['solver_kernel_ibeg']))
            inparam_advanced_input.append('KERNEL_IEND             %s\n' %(test_param['solver_kernel_iend']))
            inparam_advanced_input.append('NR_LIN_SOLIDS           %s\n' %(test_param['solver_nr_lin_solids']))
            inparam_advanced_input.append('F_MIN                   %s\n' %(test_param['solver_fmin']))
            inparam_advanced_input.append('F_MAX                   %s\n' %(test_param['solver_fmax']))
            inparam_advanced_input.append('F_REFERENCE             %s\n' %(test_param['solver_fref']))
            inparam_advanced_input.append('SMALL_Q_CORRECTION      %s\n' %(test_param['solver_small_q_correction']))
            inparam_advanced_input.append('NR_F_SAMPLE             %s\n' %(test_param['solver_nr_f_sample']))
            inparam_advanced_input.append('MAXINT_SA               %s\n' %(test_param['solver_maxint_sa']))
            inparam_advanced_input.append('TSTART_SR               %s\n' %(test_param['solver_tstart_sr']))
            inparam_advanced_input.append('TSTART_AMP              %s\n' %(test_param['solver_tstart_amp']))
            inparam_advanced_input.append('T_DECAY                 %s\n' %(test_param['solver_t_decay']))
            inparam_advanced_input.append('FIX_FREQ                %s\n' %(test_param['solver_fix_freq']))
            inparam_advanced_input.append('DUMP_VTK                %s\n' %(test_param['solver_dump_vtk']))
            inparam_advanced_input.append('COARSE_GRAINED          %s\n' %(test_param['solver_coarse_grained']))
            inparam_advanced_input.append('SAVE_ENERGY             %s\n' %(test_param['solver_save_energy']))
            inparam_advanced_input.append('HOMO_MODEL              %s\n' %(test_param['solver_homo_model']))
            inparam_advanced_input.append('HOMO_VP                 %s\n' %(test_param['solver_homo_vp']))
            inparam_advanced_input.append('HOMO_VS                 %s\n' %(test_param['solver_homo_vs']))
            inparam_advanced_input.append('HOMO_RHO                %s\n' %(test_param['solver_homo_rho']))
            inparam_advanced_input.append('FORCE_ANISO             %s\n' %(test_param['solver_force_aniso']))

            inparam_solver_open = open('inparam_advanced', 'w')
            for i in range(0, len(inparam_advanced_input)):
                inparam_solver_open.write(inparam_advanced_input[i])
            inparam_solver_open.close()

            if test_param['sourcefile_type'] == 'sourceparams':
                source_text = []
                source_text.append("SOURCE_TYPE  "     + test_param['source_type'] + "\n")
                source_text.append("SOURCE_DEPTH "     + test_param['source_depth'] + "\n")
                source_text.append("SOURCE_LAT   "     + test_param['source_lat'] + "\n")
                source_text.append("SOURCE_LON   "     + test_param['source_lon'] + "\n")
                source_text.append("SOURCE_AMPLITUDE " + test_param['source_amp'] + "\n")

                source_open = open('inparam_source', 'w')

                for i in range(0, len(source_text)):
                    source_open.write(source_text[i])

                source_open.close()
                if test_param['verbose'] != 'N':
                    for i in range(0, len(source_text)):
                        print source_text[i];
                else:
                    print 'DONE'

            elif test_param['sourcefile_type'] == 'cmtsolut':
                if test_param['verbose'] != 'N':
                    print "\n======================================"
                    print "Change the Source params (CMTSOLUTION)"
                    print "======================================"
                else:
                    sys.stdout.write('Change the Source params (CMTSOLUTION)...')
                    sys.stdout.flush()

                subprocess.check_call(['cp', 'CMTSOLUTION.TEMPLATE', 'CMTSOLUTION'])

                source_open = open('./CMTSOLUTION', 'r')
                source_read = source_open.readlines()

                source_read[0] = "PyAxi generated" + '\n'
                source_read[4] = 'latitude:      ' + test_param['cmt_lat'] + '\n'
                source_read[5] = 'longitude:     ' + test_param['cmt_lon'] + '\n'
                source_read[6] = 'depth:         ' + test_param['cmt_depth'] + '\n'
                source_read[7] = 'Mrr:      ' + test_param['cmt_Mrr'] + '\n'
                source_read[8] = 'Mtt:       ' + test_param['cmt_Mtt'] + '\n'
                source_read[9] = 'Mpp:      ' + test_param['cmt_Mpp'] + '\n'
                source_read[10] = 'Mrt:      ' + test_param['cmt_Mrt'] + '\n'
                source_read[11] = 'Mrp:       ' + test_param['cmt_Mrp'] + '\n'
                source_read[12] = 'Mtp:      ' + test_param['cmt_Mtp'] + '\n'

                source_open.close()
                source_open = open('./CMTSOLUTION', 'w')

                for i in range(0, len(source_read)):
                    source_open.write(source_read[i])

                source_open.close()
                if test_param['verbose'] != 'N':
                    print source_read[0] + source_read[4] + source_read[5] + \
                                source_read[6] + source_read[7] + \
                                source_read[8] + source_read[9] + \
                                source_read[10] + source_read[11] + \
                                source_read[12]
                else:
                    print 'DONE'

            if test_param['receiver_type'] == 'colatlon':
                if test_param['verbose'] != 'N':
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

            elif test_param['receiver_type'] == 'stations':
                if test_param['verbose'] != 'N':
                    print "\n====================================="
                    print "Change the Receiver params (STATIONS)"
                    print "====================================="
                else:
                    sys.stdout.write('Change the Receiver params (STATIONS)...')
                    sys.stdout.flush()

                if test_param['IO_STATION']:
                    subprocess.check_call(['cp', test_param['IO_STATION'], '.'])


###################### read_input_file #################################

def read_input_file():
    """
    Read inputs from inpython.cfg file.
    """

    global test_param
    config = ConfigParser.RawConfigParser()
    test_param = {}
    config.read(os.path.join(sys.argv[1]))
    test_param['inpython_address'] = os.path.join(sys.argv[1])
    print '\n****************************************'
    print 'Read the input file (inpython.cfg) from:'
    print os.path.join(sys.argv[1])

    test_param['IO_STATION'] = os.path.join(sys.argv[2])
    print '\n****************************************'
    print 'Copy the STATIONS file from:'
    print os.path.join(sys.argv[2])
    if not os.path.isabs(test_param['IO_STATION']):
        test_param['IO_STATION'] = os.path.join(os.getcwd(), test_param['IO_STATION'])

    test_param['axi_address'] = config.get('GENERAL', 'AXISEM_DIR')
    if not os.path.isabs(test_param['axi_address']):
        test_param['axi_address'] = os.path.join(os.getcwd(), test_param['axi_address'])
    print '\nWorking Directory:'
    print test_param['axi_address']
    print '****************************************\n'

    test_param['solver_name'] = config.get('GENERAL', 'SOLVER_NAME')
    test_param['verbose'] = config.get('GENERAL', 'VERBOSE')

    test_param['new_mesh'] = config.get('GENERAL', 'NEW_MESH')
    test_param['post_processing'] = config.get('GENERAL', 'POST_PROCESSING')

    test_param['mesher'] = config.get('GENERAL', 'MESHER')
    test_param['solver'] = config.get('GENERAL', 'SOLVER')

    test_param['mesher_makefile'] = config.get('GENERAL', 'MESHER_MAKEFILE')
    test_param['mesher_make'] = config.get('GENERAL', 'MESHER_MAKE')
    test_param['mesher_move'] = config.get('GENERAL', 'MESHER_MOVE')

    test_param['solver_cp'] = config.get('GENERAL', 'SOLVER_COPY')
    test_param['solver_makefile'] = config.get('GENERAL', 'SOLVER_MAKEFILE')
    test_param['solver_make'] = config.get('GENERAL', 'SOLVER_MAKE')

    if test_param['mesher_make'] != 'N':
        test_param['mesher_move'] = 'Y'
        test_param['solver_cp'] = 'Y'

    test_param['mesher_bg_model'] = config.get('MESHER_BASIC', 'BACKGROUND_MODEL')
    test_param['mesher_ext_model'] = config.get('MESHER_BASIC', 'EXT_MODEL')
    test_param['mesher_dominant_period'] = config.get('MESHER_BASIC', 'DOMINANT_PERIOD')
    test_param['mesher_ntheta'] = config.get('MESHER_BASIC', 'NTHETA_SLICES')
    test_param['mesher_nradial'] = config.get('MESHER_BASIC', 'NRADIAL_SLICES')
    test_param['mesher_write_vtk'] = config.get('MESHER_BASIC', 'WRITE_VTK')
    test_param['mesher_coarsening_layers'] = config.get('MESHER_BASIC', 'COARSENING_LAYERS')
    test_param['mesh_name'] = config.get('MESHER_BASIC', 'MESHNAME')
    test_param['mesher_ic_shear_wave'] = config.get('MESHER_ADVANCED', 'IC_SHEAR_WAVE')
    test_param['mesher_npol'] = config.get('MESHER_ADVANCED', 'NPOL')
    test_param['mesher_el_per_lambda'] = config.get('MESHER_ADVANCED', 'EL_PER_LAMBDA')
    test_param['mesher_courant_nr'] = config.get('MESHER_ADVANCED', 'COURANT_NR')
    test_param['mesher_radius'] = config.get('MESHER_ADVANCED', 'RADIUS')
    test_param['mesher_save_mesh'] = config.get('MESHER_ADVANCED', 'SAVE_MESH')
    test_param['mesher_verbose'] = config.get('MESHER_ADVANCED', 'VERBOSE')

    test_param['solver_sim_type'] = config.get('SOLVER_BASIC', 'SIMULATION_TYPE')
    if test_param['solver_sim_type'] == 'single':
        test_param['sourcefile_type'] = 'sourceparams'
    elif test_param['solver_sim_type'] == 'moment':
        test_param['sourcefile_type'] = 'cmtsolut'
    else:
        print 'Check your simulation type, you entered:'
        print test_param['solver_simu_type']
    test_param['solver_seis_length'] = config.get('SOLVER_BASIC', 'SEISMOGRAM_LENGTH')
    test_param['solver_recfile_type'] = config.get('STATION_INFO', 'RECFILE_TYPE')
    test_param['solver_lat_heterogeneity'] = config.get('SOLVER_BASIC', 'LAT_HETEROGENEITY')
    test_param['solver_attenuation'] = config.get('SOLVER_BASIC', 'ATTENUATION')
    test_param['solver_save_snapshots'] = config.get('SOLVER_BASIC', 'SAVE_SNAPSHOTS')
    test_param['solver_verbosity'] = config.get('SOLVER_BASIC', 'VERBOSITY')


    test_param['solver_sampling_rate'] = config.get('SOLVER_ADVANCED', 'SAMPLING_PERIOD')
    test_param['solver_time_step'] = config.get('SOLVER_ADVANCED', 'TIME_STEP')
    test_param['solver_source_period'] = config.get('SOLVER_ADVANCED', 'SOURCE_PERIOD')
    test_param['solver_time_scheme'] = config.get('SOLVER_ADVANCED', 'TIME_SCHEME')
    test_param['solver_data_dir'] = config.get('SOLVER_ADVANCED', 'DATA_DIR')
    test_param['solver_info_dir'] = config.get('SOLVER_ADVANCED', 'INFO_DIR')
    test_param['solver_diag_file_output'] = config.get('SOLVER_ADVANCED', 'DIAGNOSTIC_FILE_OUTPUT')
    test_param['solver_mesh_test'] = config.get('SOLVER_ADVANCED', 'MESH_TEST')
    test_param['solver_deflate_level'] = config.get('SOLVER_ADVANCED', 'DEFLATE_LEVEL')
    test_param['solver_snapshot_dt'] = config.get('SOLVER_ADVANCED', 'SNAPSHOT_DT')
    test_param['solver_snapshots_format'] = config.get('SOLVER_ADVANCED', 'SNAPSHOTS_FORMAT')
    test_param['solver_kernel_wavefields'] = config.get('SOLVER_ADVANCED', 'KERNEL_WAVEFIELDS')
    test_param['solver_kernel_spp'] = config.get('SOLVER_ADVANCED', 'KERNEL_SPP')
    test_param['solver_kernel_source'] = config.get('SOLVER_ADVANCED', 'KERNEL_SOURCE')
    test_param['solver_kernel_ibeg'] = config.get('SOLVER_ADVANCED', 'KERNEL_IBEG')
    test_param['solver_kernel_iend'] = config.get('SOLVER_ADVANCED', 'KERNEL_IEND')
    test_param['solver_nr_lin_solids'] = config.get('SOLVER_ADVANCED', 'NR_LIN_SOLIDS')
    test_param['solver_fmin'] = config.get('SOLVER_ADVANCED', 'F_MIN')
    test_param['solver_fmax'] = config.get('SOLVER_ADVANCED', 'F_MAX')
    test_param['solver_fref'] = config.get('SOLVER_ADVANCED', 'F_REFERENCE')
    test_param['solver_small_q_correction'] = config.get('SOLVER_ADVANCED', 'SMALL_Q_CORRECTION')
    test_param['solver_nr_f_sample'] = config.get('SOLVER_ADVANCED', 'NR_F_SAMPLE')
    test_param['solver_maxint_sa'] = config.get('SOLVER_ADVANCED', 'MAXINT_SA')
    test_param['solver_tstart_sr'] = config.get('SOLVER_ADVANCED', 'TSTART_SR')
    test_param['solver_tstart_amp'] = config.get('SOLVER_ADVANCED', 'TSTART_AMP')
    test_param['solver_t_decay'] = config.get('SOLVER_ADVANCED', 'T_DECAY')
    test_param['solver_fix_freq'] = config.get('SOLVER_ADVANCED', 'FIX_FREQ')
    test_param['solver_dump_vtk'] = config.get('SOLVER_ADVANCED', 'DUMP_VTK')
    test_param['solver_coarse_grained'] = config.get('SOLVER_ADVANCED', 'COARSE_GRAINED')
    test_param['solver_save_energy'] = config.get('SOLVER_ADVANCED', 'SAVE_ENERGY')
    test_param['solver_homo_model'] = config.get('SOLVER_ADVANCED', 'HOMO_MODEL')
    test_param['solver_homo_vp'] = config.get('SOLVER_ADVANCED', 'HOMO_VP')
    test_param['solver_homo_vs'] = config.get('SOLVER_ADVANCED', 'HOMO_VS')
    test_param['solver_homo_rho'] = config.get('SOLVER_ADVANCED', 'HOMO_RHO')
    test_param['solver_force_aniso'] = config.get('SOLVER_ADVANCED', 'FORCE_ANISO')

    test_param['receiver_type'] = test_param['solver_recfile_type']
    test_param['source_type'] = config.get('SOURCE_INFO', 'SOURCE_TYPE')

    test_param['source_depth'] = config.get('SOURCE_INFO', 'SOURCE_DEPTH')
    test_param['source_lat'] = config.get('SOURCE_INFO', 'SOURCE_LATITUDE')
    test_param['source_lon'] = config.get('SOURCE_INFO', 'SOURCE_LONGITUDE')
    test_param['source_stf'] = config.get('SOURCE_INFO', 'SOURCE_STF')
    test_param['source_amp'] = config.get('SOURCE_INFO', 'SOURCE_AMPLITUDE')

    test_param['cmt_lat']   = config.get('SOURCE_INFO', 'CMT_LAT')
    test_param['cmt_lon']   = config.get('SOURCE_INFO', 'CMT_LON')
    test_param['cmt_depth'] = config.get('SOURCE_INFO', 'CMT_DEPTH')
    test_param['cmt_Mrr']   = config.get('SOURCE_INFO', 'CMT_MRR')
    test_param['cmt_Mtt']   = config.get('SOURCE_INFO', 'CMT_MTT')
    test_param['cmt_Mpp']   = config.get('SOURCE_INFO', 'CMT_MPP')
    test_param['cmt_Mrt']   = config.get('SOURCE_INFO', 'CMT_MRT')
    test_param['cmt_Mrp']   = config.get('SOURCE_INFO', 'CMT_MRP')
    test_param['cmt_Mtp']   = config.get('SOURCE_INFO', 'CMT_MTP')

    test_param['post_rec_comp_sys'] = config.get('POST_PROCESSING', 'REC_COMP_SYS')
    test_param['post_conv_period'] = config.get('POST_PROCESSING', 'CONV_PERIOD')
    test_param['post_conv_stf'] = config.get('POST_PROCESSING', 'CONV_STF')
    test_param['post_seistype'] = config.get('POST_PROCESSING', 'SEISTYPE')
    test_param['post_load_snaps'] = config.get('POST_PROCESSING', 'LOAD_SNAPS')
    test_param['post_data_dir'] = config.get('POST_PROCESSING', 'DATA_DIR')
    test_param['post_negative_time'] = config.get('POST_PROCESSING', 'NEGATIVE_TIME')
    test_param['post_3D_phi_start'] = config.get('POST_PROCESSING', '3D_PHI_START')
    test_param['post_3D_phi_end'] = config.get('POST_PROCESSING', '3D_PHI_END')
    test_param['post_3D_rtop'] = config.get('POST_PROCESSING', '3D_RTOP')
    test_param['post_3D_rbot'] = config.get('POST_PROCESSING', '3D_RBOT')
    test_param['post_3D_plot_top'] = config.get('POST_PROCESSING', '3D_PLOT_TOP')
    test_param['post_3D_plot_bot'] = config.get('POST_PROCESSING', '3D_PLOT_BOT')
    test_param['post_3D_snap_beg'] = config.get('POST_PROCESSING', '3D_SNAP_BEG')
    test_param['post_3D_snap_end'] = config.get('POST_PROCESSING', '3D_SNAP_END')
    test_param['post_3D_snap_stride'] = config.get('POST_PROCESSING', '3D_SNAP_STRIDE')

    test_param['mseed'] = config.get('MISC', 'MSEED')
    test_param['mseed_all'] = config.get('MISC', 'MSEED_ALL')
    test_param['convSTF'] = config.get('MISC', 'CONV_STF')
    test_param['halfduration'] = eval(config.get('MISC', 'HALF_DURATION'))
    test_param['filter'] = config.get('MISC', 'FILTER')
    test_param['fmin'] = eval(config.get('MISC', 'FMIN'))
    test_param['fmax'] = eval(config.get('MISC', 'FMAX'))

    test_param['make_axisem_use_netcdf'] = config.get('MAKE_AXISEM', 'USE_NETCDF')
    test_param['make_axisem_netcdf_path'] = config.get('MAKE_AXISEM', 'NETCDF_PATH')

    test_param['test'] = config.get('TEST', 'TEST')
    test_param['test_folder'] = config.get('TEST', 'TEST_FOLDER')
    test_param['plot_test'] = config.get('TEST', 'PLOT')
    test_param['save_plots_test'] = config.get('TEST', 'SAVE_PLOTS')
    test_param['plot_format_test'] = config.get('TEST', 'PLOT_FORMAT')
    test_param['test_chans'] = config.get('TEST', 'CHANS')
    test_param['test_fmin'] = config.get('TEST', 'FMIN')
    test_param['test_fmax'] = config.get('TEST', 'FMAX')
    test_param['test_half'] = config.get('TEST', 'HALF_DURATION')
    test_param['test_nstat'] = config.get('TEST', 'NSTAT')

    if test_param['new_mesh'] == 'Y':
        test_param['mesher'] = 'Y'
        test_param['solver'] = 'Y'
        test_param['mesher_makefile'] = 'Y'
        test_param['mesher_make'] = 'Y'
        test_param['mesher_move'] = 'Y'
        test_param['solver_cp'] = 'Y'
        test_param['solver_makefile'] = 'Y'
        test_param['solver_make'] = 'Y'

    if test_param['new_mesh'] == 'N':
        test_param['mesher'] = 'N'
        test_param['solver'] = 'Y'
        test_param['solver_cp'] = 'N'
        test_param['solver_makefile'] = 'N'
        test_param['solver_make'] = 'Y'

    if test_param['test'] == 'N':
        test_param['plot_test'] = 'N'
        test_param['save_plots_test'] = 'N'

    if test_param['plot_test'] != 'N':
        test_param['mesher'] = 'N'
        test_param['solver'] = 'N'
        test_param['solver_cp'] = 'N'
        test_param['solver_makefile'] = 'N'
        test_param['solver_make'] = 'N'
        test_param['post_processing'] = 'N'
        test_param['test'] = 'N'
        test_param['mseed'] = 'N'

        print '##################################'
        print "PyAxi tries to copy the data from:"
        print "(solver_name flag in inpython.cfg)"
        print test_param['solver_name']
        print '##################################'

    if test_param['receiver_type'] == 'database':
        test_param['post_processing'] = 'N'

########################################################################
########################################################################
########################################################################

if __name__ == "__main__":
    status = PyAxi()

