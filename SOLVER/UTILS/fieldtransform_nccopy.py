#!/usr/bin/env python
import argparse
import netCDF4
import os.path
from progressbar import Percentage, ProgressBar, Bar, ETA, FileTransferSpeed
import subprocess
import warnings

parser = argparse.ArgumentParser(
    description="Postprocess AxiSEM's Kernel Wavefields.")

parser.add_argument('-r', '--resume', dest='resume', action='store_true',
                    default=True, help='If copying was interrupted, resume at '
                    'last checkpoint')

parser.add_argument('-d', '--deflate_level', dest='deflate_level', type=int,
                    default=2, help='deflate level for large variables. '
                    'Default = 2')

parser.add_argument('-b', '--block_size', dest='disk_block_size', type=int,
                    default=8192, help='Disk block size used for chunking')

parser.add_argument('-c', '--cache_size', dest='cache_size_mb', type=int,
                    default=10, help='Cache size for copying large variables '
                    'in MB. Default = 10')

args = parser.parse_args()

fname_in = 'axisem_output.nc4'
fname_out = 'ordered_output.nc4'

warnings.simplefilter('always', UserWarning)

# find out AxiSEM simulation_type and build list of paths
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
ndumps = getattr(nc, "number of strain dumps")
npoints = getattr(nc, "npoints")
nc.close()

# how many gll points fit into the cache
npointsperstep = args.cache_size_mb * 1048576 / 4 / ndumps

# start putting together all options for nccopy
nc_cmd = ['nccopy']

# bufsize and chunk cache size
nc_cmd.append('-m %sM -h %sM' % (args.cache_size_mb, args.cache_size_mb))

# only include these groups (and not the Surface group)
nc_cmd.append('-G Seismograms,Snapshots,Mesh')

# only include the data of these groups (Snapshots to be copied manually)
nc_cmd.append('-g Seismograms,Mesh')

# compression
nc_cmd.append('-d %d' % (args.deflate_level,))

# chunking
chunk_gll = max(args.disk_block_size / ndumps, 1)
nc_cmd.append('-c snapshots/%d,gllpoints_all/%d' % (ndumps, chunk_gll))

# paths
nc_cmd.append('%%s/%s %%s/%s' % (fname_in, fname_out))


# execute for all rundirs
cmd = ' '.join(nc_cmd)
for p in paths:
    print 'Transforming %s%s' % (p, fname_in)
    if not (os.path.isfile(p + fname_out) and args.resume):
        cmd_p = cmd % (p, p)
        print cmd_p
        subprocess.check_call(cmd_p.split())
    else:
        print 'Output file exists, resuming copy'

    # copy large fields:
    nc_in = netCDF4.Dataset(p + fname_in, "r", format="NETCDF4")
    nc_out = netCDF4.Dataset(p + fname_out, "a", format="NETCDF4")

    for var_in in nc_in.groups['Snapshots'].variables.values():
        var_out = nc_out.groups['Snapshots'].variables[var_in.name]
        if args.resume:
            try:
                nstep = var_out.nstep
            except:
                warnings.warn(
                    'Restart unsuccessful, starting to copy from beginning')
                nstep = 0
        else:
            nstep = 0

        # start a new progressbar
        widgets = ['%s: ' % (var_out.name,), Percentage(), ' ', Bar(), ' ',
                   ETA(), ' ', FileTransferSpeed()]
        pbar = ProgressBar(widgets=widgets, maxval=ndumps * npoints * 4)
        pbar.start()

        while (nstep < npoints):
            npointread = min(npointsperstep, npoints - nstep)

            var_out[:, nstep:nstep+npointread] = \
                var_in[:, nstep:nstep+npointread]

            pbar.update(ndumps * nstep * 4)
            var_out.nstep = nstep + npointread
            nstep = nstep + npointread
        pbar.finish()

    nc_in.close()
    nc_out.close()

    print 'Done with %s%s' % (p, fname_in)
