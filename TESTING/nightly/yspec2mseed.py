#!/usr/bin/env python

import numpy as np
import sys
import glob
from obspy.core import Trace, Stream, UTCDateTime

t = UTCDateTime(0)

try:
    path = sys.argv[1]
except:
    sys.exit('usage: python yspec2mseed.py <path>')

traces = []

chans = ['Z', 'N', 'E']

for file in glob.iglob(path + 'yspec.out.*'):
    stationID = int(file.split('.')[-1])
    
    dat = np.loadtxt(file)
    npts = len(dat[:,0])

    for i, chan in enumerate(chans):
        stats = {'network': 'SG', 
                 'station': 'RS%02d' % stationID, 
                 'location': '',
                 'channel': chan, 
                 'npts': npts, 
                 'sampling_rate': (npts - 1.)/(dat[-1,0] - dat[0,0]),
                 'starttime': t,
                 'mseed' : {'dataquality': 'D'}}
        traces.append(Trace(data=dat[:,1+i], header=stats))
 
 
st = Stream(traces)
st.sort()

fname =  path + 'seismograms.mseed'

print fname
st.write(fname, format='MSEED')
