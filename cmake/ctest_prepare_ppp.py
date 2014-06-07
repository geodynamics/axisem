#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Added this line for python 2.5 compatibility
from __future__ import with_statement
import os
import sys
import ConfigParser

########################################################################
############################# Main Program #############################
########################################################################

def PyAxi(**kwargs):

    # global variables
    global test_param
    # ------------------Read INPUT file (Parameters)--------------------
    read_input_file()

    if test_param['post_processing'] != 'N':
        if test_param['sourcefile_type'] == 'sourceparams':
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
                if test_param['verbose'] != 'N': print post_process_read[i].split('\n')[0]
            if test_param['verbose'] != 'N': print 2*"======================================"
            post_process_open.close()

        elif test_param['sourcefile_type'] == 'cmtsolut':
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
                if test_param['verbose'] != 'N': print post_process_read[i].split('\n')[0]
            if test_param['verbose'] != 'N': print 2*"======================================"
            post_process_open.close()


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

    print '****************************************\n'

    test_param['solver_name'] = config.get('GENERAL', 'SOLVER_NAME')
    test_param['verbose'] = config.get('GENERAL', 'VERBOSE')
    test_param['post_processing'] = config.get('GENERAL', 'POST_PROCESSING')

    test_param['solver_sim_type'] = config.get('SOLVER_BASIC', 'SIMULATION_TYPE')
    if test_param['solver_sim_type'] == 'single':
        test_param['sourcefile_type'] = 'sourceparams'
    elif test_param['solver_sim_type'] == 'moment':
        test_param['sourcefile_type'] = 'cmtsolut'
    else:
        print 'Check your simulation type, you entered:'
        print test_param['solver_simu_type']

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

###################### edit_param_post_processing ######################
def edit_param_post_processing(post_process_read):
    """
    edit param_post_processing file
    """
    global test_param

    if test_param['post_rec_comp_sys'] != 'DNC':
        post_process_read[0] = \
            'REC_COMP_SYS    %s\n' %(test_param['post_rec_comp_sys'])
    if test_param['post_conv_period'] != 'DNC':
        post_process_read[1] = \
            'CONV_PERIOD     %s\n' %(test_param['post_conv_period'])
    if test_param['post_conv_stf'] != 'DNC':
        post_process_read[2] = \
            'CONV_STF        %s\n' %(test_param['post_conv_stf'])
    if test_param['post_seistype'] != 'DNC':
        post_process_read[3] = \
            'SEISTYPE        %s\n' %(test_param['post_seistype'])
    if test_param['post_load_snaps'] != 'DNC':
        post_process_read[4] = \
            'LOAD_SNAPS      %s\n' %(test_param['post_load_snaps'])
    if test_param['post_data_dir'] != 'DNC':
        post_process_read[5] = \
            'DATA_DIR        %s\n' %(test_param['post_data_dir'])
    if test_param['post_negative_time'] != 'DNC':
        post_process_read[6] = \
            'NEGATIVE_TIME   %s\n' %(test_param['post_negative_time'])
    if test_param['post_3D_phi_start'] != 'DNC':
        post_process_read[7] = \
            '3D_PHI_START     %s\n' %(test_param['post_3D_phi_start'])
    if test_param['post_3D_phi_end'] != 'DNC':
        post_process_read[8] = \
            '3D_PHI_END      %s\n' %(test_param['post_3D_phi_end'])
    if test_param['post_3D_rtop'] != 'DNC':
        post_process_read[9] = \
            '3D_RTOP         %s\n' %(test_param['post_3D_rtop'])
    if test_param['post_3D_rbot'] != 'DNC':
        post_process_read[10] = \
            '3D_RBOT         %s\n' %(test_param['post_3D_rbot'])
    if test_param['post_3D_plot_top'] != 'DNC':
        post_process_read[11] = \
            '3D_PLOT_TOP     %s\n' %(test_param['post_3D_plot_top'])
    if test_param['post_3D_plot_bot'] != 'DNC':
        post_process_read[12] = \
            '3D_PLOT_BOT     %s\n' %(test_param['post_3D_plot_bot'])
    if test_param['post_3D_snap_beg'] != 'DNC':
        post_process_read[13] = \
            '3D_SNAP_BEG      %s\n' %(test_param['post_3D_snap_beg'])
    if test_param['post_3D_snap_end'] != 'DNC':
        post_process_read[14] = \
            '3D_SNAP_END      %s\n' %(test_param['post_3D_snap_end'])
    if test_param['post_3D_snap_stride'] != 'DNC':
        post_process_read[15] = \
            '3D_SNAP_STRIDE   %s\n' %(test_param['post_3D_snap_stride'])

    return post_process_read

########################################################################
########################################################################
########################################################################

if __name__ == "__main__":
    status = PyAxi()

