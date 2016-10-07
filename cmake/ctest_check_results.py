#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
    global test_param
    # ------------------Read INPUT file (Parameters)--------------------
    read_input_file()

    ##############################################################
    ############################# MSEED ##########################
    ##############################################################

    if test_param['mseed'] != 'N':
        path = os.path.join('Data_Postprocessing', 'SEISMOGRAMS')
        axisem2mseed(path = path)

    if test_param['mseed_all'] != 'N':
        path = os.path.join('Data_Postprocessing', 'SEISMOGRAMS')
        axisem2mseed_all(path = path)

    ##############################################################
    ############################# TEST ###########################
    ##############################################################

    if test_param['test'] != 'N' or test_param['plot_test'] != 'N':

        chans = eval(test_param['test_chans'])
        fmin = eval(test_param['test_fmin'])
        fmax = eval(test_param['test_fmax'])
        halfduration = eval(test_param['test_half'])
        ref_folder = sys.argv[2]

        try:
            os.makedirs(os.path.join('new_data'))
        except OSError:
            pass
        test_result_folder = os.path.join('new_data')

        subprocess.check_call(['cp',
                os.path.join('Data_Postprocessing', 'SEISMOGRAMS', 'seismograms.mseed'),
                                os.path.join(test_result_folder, 'axisem.mseed')])

        sgs = []
        st = read(os.path.join(ref_folder, 'axisem.mseed'))
        sgs.append(st)
        st = read(os.path.join(ref_folder, 'yspec.mseed'))
        sgs.append(st)
        st = read(os.path.join(test_result_folder, 'axisem.mseed'))
        sgs.append(st)

        sigma =  halfduration / np.sqrt(2.) / 3.5

        for st in sgs:
            convSTF(st, sigma=sigma)
        for sg in sgs:
            sg.filter('lowpass', freq=fmax, corners=2)
            sg.filter('highpass', freq=fmin, corners=2)

        nstat = eval(test_param['test_nstat'])
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
            fl2 = open(os.path.join(test_result_folder, 'l2misfit_%s.dat' % chan), 'w')
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

            if test_param['save_plots_test'] != 'N':
                recsec.savefig(os.path.join(test_result_folder, 'record_section_%s.%s' %
                        (chan, test_param['plot_format_test'])))

            ax = misfitplot.gca()
            ax.semilogy(np.arange(nstat) + 1, np.array(l2misfit) + 1e-12, 'o', label=chan)
            ax.set_xlabel('trace')
            ax.set_ylabel('l2 - misfit to reference data')
            ax.set_ylim(1e-12, 1.)

            ax.axhline(y=1e-8, color='k', ls='--')
            if np.max(l2misfit) > 1e-8:
                fwarn = open(os.path.join(test_result_folder, 'warning.dat'), 'a')
                fwarn.write("maximum l2 norm misfit larger then 1e-8 in chan %s trace %d\n"
                             % (chan, np.argmax(l2misfit)))
                fwarn.close()
                sys.exit("L2 norm misfit is higher than 1e-8.")
                
        ax = misfitplot.gca()
        ax.legend()
        if test_param['save_plots_test'] != 'N':
            misfitplot.savefig(os.path.join(test_result_folder, 'l2_misfit.' + test_param['plot_format_test']))
        if test_param['plot_test'] != 'N':
            plt.show()

###################### read_input_file #################################

def read_input_file():
    """
    Read inputs from inpython.cfg file.
    """

    global test_param, obspy_error
    config = ConfigParser.RawConfigParser()
    test_param = {}
    try:
        config.read(os.path.join(sys.argv[1]))
        test_param['inpython_address'] = os.path.join(sys.argv[1])
        print '\n****************************************'
        print 'Read the input file (inpython.cfg) from:'
        print os.path.join(sys.argv[1])
    except Exception, error:
        config.read(os.path.join(os.getcwd(), 'inpython.cfg'))
        test_param['inpython_address'] = os.path.join(os.getcwd(), 'inpython.cfg')
        print '\n****************************************'
        print 'Read the input file (inpython.cfg) from:'
        print os.path.join(os.getcwd(), 'inpython.cfg')

    try:
        test_param['IO_STATION'] = os.path.join(sys.argv[2])
        print '\n****************************************'
        print 'Copy the STATIONS file from:'
        print os.path.join(sys.argv[2])
        if not os.path.isabs(test_param['IO_STATION']):
            test_param['IO_STATION'] = os.path.join(os.getcwd(),
                                          test_param['IO_STATION'])
    except Exception, error:
        print 'STATIONS file is not entered as an input!'
        test_param['IO_STATION'] = False

    test_param['axi_address'] = config.get('GENERAL', 'AXISEM_DIR')
    if not os.path.isabs(test_param['axi_address']):
        test_param['axi_address'] = os.path.join(os.getcwd(), test_param['axi_address'])
    print '\nWorking Directory:'
    print test_param['axi_address']
    print '****************************************\n'

    test_param['solver_name'] = config.get('GENERAL', 'SOLVER_NAME')
    test_param['verbose'] = config.get('GENERAL', 'VERBOSE')

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

    test_param['make_axisem_mpirun'] = config.get('MAKE_AXISEM', 'MPIRUN')
    test_param['make_axisem_use_netcdf'] = config.get('MAKE_AXISEM', 'USE_NETCDF')
    test_param['make_axisem_netcdf_path'] = config.get('MAKE_AXISEM', 'NETCDF_PATH')
    test_param['make_axisem_serial'] = config.get('MAKE_AXISEM', 'SERIAL')
    test_param['make_axisem_CC'] = config.get('MAKE_AXISEM', 'CC')
    test_param['make_axisem_FC'] = config.get('MAKE_AXISEM', 'FC')
    test_param['make_axisem_FFLAGS'] = config.get('MAKE_AXISEM', 'FFLAGS')
    test_param['make_axisem_CFLAGS'] = config.get('MAKE_AXISEM', 'CFLAGS')
    test_param['make_axisem_LDFLAGS'] = config.get('MAKE_AXISEM', 'LDFLAGS')

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

    test_param['post_processing'] = 'N'

    if obspy_error != 'N':
        test_param['test'] = 'N'
        test_param['mseed'] = 'N'


########################## axisem2mseed ################################
def axisem2mseed(path):
    """
    change .dat files into MSEED format
    """
    global test_param

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
            if test_param['convSTF'] == 'Y':
                sigma =  test_param['halfduration'] / np.sqrt(2.) / 3.5
                convSTF(st, sigma=sigma)
            if test_param['filter'] == 'Y':
                st.filter('lowpass', freq=test_param['fmax'], corners=2)
                st.filter('lowpass', freq=test_param['fmax'], corners=2)
                st.filter('lowpass', freq=test_param['fmax'], corners=2)
                st.filter('lowpass', freq=test_param['fmax'], corners=2)
                st.filter('highpass', freq=test_param['fmin'], corners=2)
                st.filter('highpass', freq=test_param['fmin'], corners=2)
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
    global test_param

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
    status = PyAxi()

