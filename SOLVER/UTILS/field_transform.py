#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Postprocessing for AxiSEM's Kernel Wavefields

:copyright:
    Martin van Driel (Martin@vanDriel.de), 2015
:license:
    GNU General Public License, Version 3
    (http://www.gnu.org/licenses/gpl.html)
"""
import argparse
import netCDF4
import os.path
try:
    from progressbar import Percentage, ProgressBar, Bar, ETA, \
        FileTransferSpeed
    progressbar_installed = True
except ImportError:
    progressbar_installed = False
import subprocess
import warnings

# parse command line arguments
parser = argparse.ArgumentParser(
    description="Postprocess AxiSEM's Kernel Wavefields: rechunkink for "
                "faster time series reading, optional compression.")

parser.add_argument('-r', '--resume', dest='resume', action='store_true',
                    default=False, help='If copying was interrupted, resume '
                    'at last checkpoint')

parser.add_argument('-d', '--deflate_level', dest='deflate_level', type=int,
                    default=2, help='deflate level for large variables. 0 for '
                    'no comression, 9 for maximum compression. Default = 2')

parser.add_argument('-b', '--block_size', dest='disk_block_size', type=int,
                    default=8192, help='Disk block size used for chunking')

parser.add_argument('-c', '--cache_size', dest='cache_size_mb', type=int,
                    default=10, help='Cache size for copying large variables '
                    'in MB. Default = 10')

args = parser.parse_args()

fname_in = 'axisem_output.nc4'
fname_out = 'ordered_output.nc4'

# warn for each variable when problems with resume occure
warnings.simplefilter('always', UserWarning)

# find out AxiSEM SIMULATION_TYPE and build list of paths for the wavefield
# files
with open('inparam_basic', 'r') as f:
    for line in f:
        if line.startswith('SIMULATION_TYPE'):
            simulation_type = line.split()[1]
            break

simulation_type_map = {'single': ['.'],
                       'force': ['PZ', 'PX'],
                       'moment': ['MZZ', 'MXX_P_MYY', 'MXZ_MYZ',
                                  'MXY_MXX_M_MYY']}

paths = ['%s/Data/' % st for st in simulation_type_map[simulation_type]]

# get some numbers to determine chunking from the global attributes
nc = netCDF4.Dataset(paths[0] + fname_in, "r", format="NETCDF4")

try:
    ndumps = getattr(nc, "number of strain dumps")
    npoints = getattr(nc, "npoints")
    nc.close()
except AttributeError:
    print 'Wavefield output for kernels/Instaseis was not switched for the ' + \
          'simulation! Check Option KERNEL_WAVEFIELDS in inparam_advanced!'
    raise RuntimeError


# how many gll points fit into the cache
npointsperstep = args.cache_size_mb * 1048576 / 4 / ndumps

# start putting together all options for nccopy
nc_cmd = ['nccopy']

# bufsize and chunk cache size
nc_cmd.append('-m %sM -h %sM' % (args.cache_size_mb, args.cache_size_mb))

# only include the data of these groups (Snapshots to be copied manually)
nc_cmd.append('-g Seismograms,Mesh')

# only include these groups (and not the Surface group)
nc_cmd.append('-G Seismograms,Snapshots,Mesh')

# compression
nc_cmd.append('-d %d' % (args.deflate_level,))

# chunking
chunk_gll = max(args.disk_block_size / ndumps, 1)
nc_cmd.append('-c snapshots/%d,gllpoints_all/%d' % (ndumps, chunk_gll))

# paths
nc_cmd.append('%%s%s %%s%s' % (fname_in, fname_out))

# join to a single command
cmd = ' '.join(nc_cmd)

# execute for all rundirs
for p in paths:
    print 'Transforming %s%s' % (p, fname_in)

    if not (os.path.isfile(p + fname_out) and args.resume):
        # run nccopy
        cmd_p = cmd % (p, p)
        print cmd_p
        subprocess.check_call(cmd_p.split())
    else:
        print 'Output file exists, resuming copy'

    # open files for manual copy ('r' for read, 'a' for append)
    nc_in = netCDF4.Dataset(p + fname_in, "r", format="NETCDF4")
    nc_out = netCDF4.Dataset(p + fname_out, "a", format="NETCDF4")

    # copy large fields
    for var_in in nc_in.groups['Snapshots'].variables.values():
        var_out = nc_out.groups['Snapshots'].variables[var_in.name]

        if var_in.name in ('stf_dump', 'stf_d_dump'):
            # Not really a large field, is copied in one go
            var_out[:] = var_in[:]

        else:

            if args.resume:
                # try to resume, might fail to read nstep from the attribute
                try:
                    nstep = var_out.nstep
                except:
                    warnings.warn('Restart unsuccessful, starting to copy '
                                  'from beginning')
                    nstep = 0
            else:
                nstep = 0

            if progressbar_installed:
                # start a new progressbar
                widgets = ['%s: ' % (var_out.name,), Percentage(), ' ', Bar(),
                           ' ', ETA(), ' ', FileTransferSpeed()]

                # convert to floats to avoid buffer overflow
                maxval = float(ndumps) * float(npoints)

                pbar = ProgressBar(widgets=widgets, maxval=maxval)
                pbar.start()

            # copy large fields chunkwise
            while (nstep < npoints):
                npointread = min(npointsperstep, npoints - nstep)

                var_out[:, nstep:nstep+npointread] = \
                    var_in[:, nstep:nstep+npointread]

                if progressbar_installed:
                    pbar.update(float(ndumps) * float(nstep))

                # set a checkpoint to variable attribute
                var_out.nstep = nstep + npointread

                nstep = nstep + npointread

            if progressbar_installed:
                pbar.finish()

    # close files
    nc_in.close()
    nc_out.close()

    print 'Done with %s%s' % (p, fname_in)

print 'DONE'
