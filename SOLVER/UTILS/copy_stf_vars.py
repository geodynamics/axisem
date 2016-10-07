
# coding: utf-8

import netCDF4
import argparse

helptext = 'Copy variables stf_dump and stf_d_dump from one AxiSEM file \n' + \
           'to another. Move it from the Surface group to the Snapshots \n' + \
           'group at the same time'
formatter_class = argparse.RawTextHelpFormatter
parser = argparse.ArgumentParser(description=helptext,
                                 formatter_class=formatter_class)

helptext = 'NetCDF input file name.'
parser.add_argument('input_file_name', help=helptext)

helptext = 'NetCDF output file name.'
parser.add_argument('output_file_name', help=helptext)

args = parser.parse_args()


fnam_axisem_output = args.input_file_name
fnam_ordered_output = args.output_file_name


with netCDF4.Dataset(fnam_axisem_output, 'r') as nc_ax:
    with netCDF4.Dataset(fnam_ordered_output, 'r+') as nc_or:
        nc_ax_surf = nc_ax.groups['Surface']
        nc_or_snap = nc_or.groups['Snapshots']
        stf_dump_ax = nc_ax_surf.variables['stf_dump']

        nc_or_snap.createVariable(varname='stf_dump', datatype='float32',
                                  dimensions=('snapshots'))
        stf_dump_or = nc_or_snap.variables['stf_dump']

        stf_dump_or[:] = 0.0
        stf_dump_or[0:stf_dump_ax.shape[0]] = stf_dump_ax[:]

        stf_dump_ax = nc_ax_surf.variables['stf_d_dump']

        nc_or_snap.createVariable(varname='stf_d_dump', datatype='float32',
                                  dimensions=('snapshots'))
        stf_dump_or = nc_or_snap.variables['stf_d_dump']

        stf_dump_or[:] = 0.0
        stf_dump_or[0:stf_dump_ax.shape[0]] = stf_dump_ax[:]
