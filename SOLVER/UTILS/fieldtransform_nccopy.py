#!/usr/bin/env python
import netCDF4
import subprocess

DISK_BLOCK_SIZE = 8192

fname_in = 'axisem_output.nc4'
fname_out = 'ordered_output.nc4'

path_px = 'PX/Data/'
path_pz = 'PZ/Data/'

# get some numbers to determine chunking from the global attributes
nc = netCDF4.Dataset(path_px + fname_in, "r", format="NETCDF4")
ndumps = getattr(nc, "number of strain dumps")
nc.close()

# compute chunking
chunk_gll = max(DISK_BLOCK_SIZE / ndumps, 1)

cmd = "nccopy -d 2 -c snapshots/%d,gllpoints_all/%d %s/%s %s/%s"
cmd_px = cmd % (ndumps, chunk_gll, path_px, fname_in, path_px, fname_out)
cmd_pz = cmd % (ndumps, chunk_gll, path_pz, fname_in, path_pz, fname_out)

print "Transforming PX"
subprocess.call(cmd_px.split())

print "Transforming PZ"
subprocess.call(cmd_pz.split())
