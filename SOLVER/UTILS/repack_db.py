#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Repacking Instaseis databases.

Requires click, netCDF4, and numpy.

:copyright:
    Lion Krischer (krischer@geophysik.uni-muenchen.de), 2016
:license:
    GNU Lesser General Public License, Version 3 [non-commercial/academic use]
    (http://www.gnu.org/copyleft/lgpl.html)
"""
import contextlib
import math
import os
import sys

import click
import netCDF4
import numpy as np
from scipy.spatial import cKDTree


if sys.version_info.major == 2:
    str_type = (basestring, str, unicode)  # NOQA
else:
    str_type = (bytes, str)


__netcdf_version = tuple(int(i) for i in netCDF4.__version__.split("."))


@contextlib.contextmanager
def dummy_progressbar(iterator, *args, **kwargs):
    yield iterator


def repack_file(input_filename, output_filename, contiguous,
                compression_level, transpose, quiet=False):
    """
    Transposes all data in the "/Snapshots" group.

    :param input_filename: The input filename.
    :param output_filename: The output filename.
    """
    assert os.path.exists(input_filename)
    assert not os.path.exists(output_filename)

    with netCDF4.Dataset(input_filename, "r", format="NETCDF4") as f_in, \
            netCDF4.Dataset(output_filename, "w", format="NETCDF4") as f_out:
        recursive_copy(src=f_in, dst=f_out, contiguous=contiguous,
                       compression_level=compression_level, quiet=quiet,
                       transpose=transpose)


def recursive_copy(src, dst, contiguous, compression_level, transpose, quiet):
    """
    Recursively copy the whole file and transpose the all /Snapshots
    variables while at it..
    """
    if src.path == "/Seismograms":
        return

    for attr in src.ncattrs():
        _s = getattr(src, attr)
        if isinstance(_s, str_type):
            # The setncattr_string() was added in version 1.2.3. Before that
            # it was the default behavior.
            if __netcdf_version >= (1, 2, 3):
                dst.setncattr_string(attr, _s)
            else:
                dst.setncattr(attr, str(_s))
        else:
            setattr(dst, attr, _s)

    # We will only transpose the Snapshots group.
    if src.path == "/Snapshots":
        is_snap = True
    else:
        is_snap = False

    # The dimensions will be reversed for the snapshots.
    items = list(src.dimensions.items())
    if is_snap and transpose:
        items = list(reversed(items))

    for name, dimension in items:
        dst.createDimension(name, len(
            dimension) if not dimension.isunlimited() else None)

    _j = 0
    for name, variable in src.variables.items():
        _j += 1
        shape = variable.shape

        # Determine chunking - only for the snapshots.
        if is_snap and name.startswith("disp_"):
            npts = min(shape)
            num_elems = max(shape)
            time_axis = np.argmin(shape)
            # Arbitrary limit.
            _c = max(int(round(32768 / (npts * 4))), 1)

            if time_axis == 0:
                chunksizes = (npts, _c)
            else:
                chunksizes = (_c, npts)

            if transpose:
                chunksizes = list(reversed(chunksizes))
        else:
            # For non-snapshots, just use the existing chunking.
            chunksizes = variable.chunking()
            # We could infer the chunking here but I'm not sure its worth it.
            if isinstance(chunksizes, str_type) and chunksizes == "contiguous":
                chunksizes = None

        # For a contiguous output, compression and chunking has to be turned
        # off.
        if contiguous:
            zlib = False
            chunksizes = None
        else:
            zlib = True

        dimensions = variable.dimensions
        if is_snap and transpose:
            dimensions = list(reversed(dimensions))

        x = dst.createVariable(name, variable.datatype, dimensions,
                               chunksizes=chunksizes, contiguous=contiguous,
                               zlib=zlib, complevel=compression_level)
        # Non-snapshots variables are just copied in a single go.
        if not is_snap or not name.startswith("disp_"):
            if not quiet:
                click.echo(click.style("\tCopying group '%s'..." % name,
                                       fg="blue"))
            dst.variables[x.name][:] = src.variables[x.name][:]
        # The snapshots variables are incrementally copied and transposed.
        else:
            if not quiet:
                click.echo(click.style(
                    "\tCopying 'Snapshots/%s' (%i of %i)..." % (
                        name, _j,
                        len([_i for _i in src.variables
                             if _i.startswith("disp_")])),
                    fg="blue"))

            # Copy around 8 Megabytes at a time. This seems to be the
            # sweet spot at least on my laptop.
            factor = int((8 * 1024 * 1024 / 4) / npts)
            s = int(math.ceil(num_elems / float(factor)))

            if quiet:
                pbar = dummy_progressbar
            else:
                pbar = click.progressbar

            with pbar(range(s), length=s, label="\t  ") as idx:
                for _i in idx:
                    _s = slice(_i * factor, _i * factor + factor)
                    if transpose:
                        if time_axis == 0:
                            dst.variables[x.name][_s, :] = \
                                src.variables[x.name][:, _s].T
                        else:
                            dst.variables[x.name][:, _s] = \
                                src.variables[x.name][_s, :].T
                    else:
                        if time_axis == 0:
                            dst.variables[x.name][:, _s] = \
                                src.variables[x.name][:, _s]
                        else:
                            dst.variables[x.name][_s, :] = \
                                src.variables[x.name][_s, :]

    for src_group in src.groups.values():
        dst_group = dst.createGroup(src_group.name)
        recursive_copy(src=src_group, dst=dst_group, contiguous=contiguous,
                       compression_level=compression_level, quiet=quiet,
                       transpose=transpose)


def recursive_copy_no_snapshots_no_seismograms_no_surface(
        src, dst, quiet, contiguous, compression_level):
    """
    A bit of a copy of the recursive_copy function but it does not copy the
    Snapshots, Seismograms, or Surface group.
    """
    for attr in src.ncattrs():
        _s = getattr(src, attr)
        if isinstance(_s, str_type):
            # The setncattr_string() was added in version 1.2.3. Before that
            # it was the default behavior.
            if __netcdf_version >= (1, 2, 3):
                dst.setncattr_string(attr, _s)
            else:
                dst.setncattr(attr, str(_s))
        else:
            setattr(dst, attr, _s)

    items = list(src.dimensions.items())

    for name, dimension in items:
        dst.createDimension(name, len(
            dimension) if not dimension.isunlimited() else None)

    for name, variable in src.variables.items():
        if name in ["Snapshots", "Seismograms", "Surface"]:
            continue

        # Use the existing chunking.
        chunksizes = variable.chunking()
        # We could infer the chunking here but I'm not sure its worth it.
        if isinstance(chunksizes, str_type) and chunksizes == "contiguous":
            chunksizes = None

        # For a contiguous output, compression and chunking has to be turned
        # off.
        if contiguous:
            zlib = False
            chunksizes = None
        else:
            zlib = True

        dimensions = variable.dimensions

        x = dst.createVariable(name, variable.datatype, dimensions,
                               chunksizes=chunksizes, contiguous=contiguous,
                               zlib=zlib, complevel=compression_level)
        if not quiet:
            click.echo(click.style("\tCopying group '%s'..." % name,
                                   fg="blue"))
        dst.variables[x.name][:] = src.variables[x.name][:]

    for src_group in src.groups.values():
        if src_group.name in ["Snapshots", "Seismograms", "Surface"]:
            continue
        dst_group = dst.createGroup(src_group.name)
        recursive_copy_no_snapshots_no_seismograms_no_surface(
            src=src_group, dst=dst_group, contiguous=contiguous,
            compression_level=compression_level, quiet=quiet)


def merge_files(filenames, output_folder, contiguous, compression_level,
                quiet):
    """
    Completely unroll and merge both files to a single database.
    """
    assert len(filenames) in (1, 2, 4)

    files = {}
    for file in filenames:
        for k in ("PX", "PZ", "MXX_P_MYY", "MXY_MXX_M_MYY", "MXZ_MYZ", "MZZ"):
            if k in file:
                assert k not in files
                assert os.path.exists(file)
                files[k] = file

    # Only a couple of combinations are valid.
    keys = sorted(files.keys())
    assert (keys == ["PX"]) or (keys == ["PZ"]) or (keys == ["PX", "PZ"]) or \
        (keys == ["MXX_P_MYY", "MXY_MXX_M_MYY", "MXZ_MYZ", "MZZ"])

    output = os.path.join(output_folder, "merged_output.nc4")
    assert not os.path.exists(output)

    input_files = {}
    try:
        for key, value in files.items():
            input_files[key] = netCDF4.Dataset(value, "r", format="NETCDF4")
        out = netCDF4.Dataset(output, "w", format="NETCDF4")
        _merge_files(input=input_files, out=out, contiguous=contiguous,
                     compression_level=compression_level, quiet=quiet)
    finally:
        for filename in input_files.values():
            try:
                filename.close()
            except:
                pass
        try:
            out.close()
        except:
            pass


def _merge_files(input, out, contiguous, compression_level, quiet):
    # First copy everything non-snapshot related.
    c_db = list(input.values())[0]
    recursive_copy_no_snapshots_no_seismograms_no_surface(
        src=c_db, dst=out, quiet=quiet, contiguous=contiguous,
        compression_level=compression_level)

    if contiguous:
        zlib = False
    else:
        zlib = True

    # We need the stf_dump and stf_d_dump datasets. They are either in the
    # "Snapshots" group or in the "Surface" group.
    for g in ("Snapshots", "Surface"):
        if g not in c_db.groups:
            continue
        if "stf_dump" not in c_db[g].variables:
            continue
        break
    else:
        raise Exception("Could not find `stf_dump` array.")

    stf_dump = c_db[g]["stf_dump"]
    stf_d_dump = c_db[g]["stf_d_dump"]

    for data in [stf_dump, stf_d_dump]:
        chunksizes = data.shape
        if contiguous:
            chunksizes = None
        d = out.createVariable(
            varname=data.name,
            dimensions=["snapshots"],
            contiguous=contiguous,
            zlib=zlib,
            chunksizes=chunksizes,
            datatype=data.dtype)
        d[:] = data[:]

    # Get all the snapshots from the other databases.
    if "PX" in input and "PZ" in input:
        meshes = [
            input["PX"]["Snapshots"]["disp_s"],
            input["PX"]["Snapshots"]["disp_p"],
            input["PX"]["Snapshots"]["disp_z"],
            input["PZ"]["Snapshots"]["disp_s"],
            input["PZ"]["Snapshots"]["disp_z"]]
    elif "PX" in input and "PZ" not in input:
        meshes = [
            input["PX"]["Snapshots"]["disp_s"],
            input["PX"]["Snapshots"]["disp_p"],
            input["PX"]["Snapshots"]["disp_z"]]
    elif "PZ" in input and "PX" not in input:
        meshes = [
            input["PZ"]["Snapshots"]["disp_s"],
            input["PZ"]["Snapshots"]["disp_z"]]
    elif "MXX_P_MYY" in input and "MXY_MXX_M_MYY" in input and \
            "MXZ_MYZ" in input and "MZZ" in input:
        meshes = [
            input["MZZ"]["Snapshots"]["disp_s"],
            input["MZZ"]["Snapshots"]["disp_z"],
            input["MXX_P_MYY"]["Snapshots"]["disp_s"],
            input["MXX_P_MYY"]["Snapshots"]["disp_z"],
            input["MXZ_MYZ"]["Snapshots"]["disp_s"],
            input["MXZ_MYZ"]["Snapshots"]["disp_p"],
            input["MXZ_MYZ"]["Snapshots"]["disp_z"],
            input["MXY_MXX_M_MYY"]["Snapshots"]["disp_s"],
            input["MXY_MXX_M_MYY"]["Snapshots"]["disp_p"],
            input["MXY_MXX_M_MYY"]["Snapshots"]["disp_z"]]
    else:  # pragma: no cover
        raise NotImplementedError

    time_axis = np.argmin(meshes[0].shape)

    dtype = meshes[0].dtype

    # Create new dimensions.
    dim_ipol = out.createDimension("ipol", 5)
    dim_jpol = out.createDimension("jpol", 5)
    dim_nvars = out.createDimension("nvars", len(meshes))
    nelem = out.getncattr("nelem_kwf_global")
    dim_elements = out.createDimension("elements", nelem)

    # New dimensions for the 5D Array.
    dims = (dim_elements, dim_nvars, dim_jpol, dim_ipol,
            out.dimensions["snapshots"])
    dimensions = [_i.name for _i in dims]

    if contiguous:
        chunksizes = None
    else:
        # Each chunk is exactly the data from one element.
        chunksizes = [_i.size for _i in dims]
        chunksizes[0] = 1

    # We'll called it MergedSnapshots
    x = out.createVariable(
        varname="MergedSnapshots",
        dimensions=dimensions,
        contiguous=contiguous,
        zlib=zlib,
        chunksizes=chunksizes,
        datatype=dtype)

    utemp = np.zeros([_i.size for _i in dims[1:]], dtype=dtype, order="C")

    # We also re-sort the elements to follow the traversal of a kd-tree in
    # the same fashion instaseis uses it - this should allow for even faster
    # I/O for spatially adjacent elements.

    # Get the midpoints for each element.
    s_mp = c_db["Mesh"]["mp_mesh_S"][:]
    z_mp = c_db["Mesh"]["mp_mesh_Z"][:]

    # Fill kd-tree.
    midpoints = np.empty((s_mp.shape[0], 2), dtype=s_mp.dtype)
    midpoints[:, 0] = s_mp[:]
    midpoints[:, 1] = z_mp[:]
    kdtree = cKDTree(data=midpoints)

    # This is now the order in which we will write the indices.
    inds = kdtree.indices

    # Make sure all indices are available.
    assert list(range(nelem)) == sorted(inds)

    sem_mesh = c_db["Mesh"]["sem_mesh"][:].copy()

    # Resort and write the new order to the file.
    out["Mesh"]["sem_mesh"][:] = sem_mesh[inds]
    out["Mesh"]["fem_mesh"][:] = out["Mesh"]["fem_mesh"][:][inds]

    # We'll also have to resort the midpoints.
    out["Mesh"]["mp_mesh_S"][:] = c_db["Mesh"]["mp_mesh_S"][:][inds]
    out["Mesh"]["mp_mesh_Z"][:] = c_db["Mesh"]["mp_mesh_Z"][:][inds]
    # And a couple of other things.
    out["Mesh"]["eltype"][:] = out["Mesh"]["eltype"][:][inds]
    out["Mesh"]["axis"][:] = out["Mesh"]["axis"][:][inds]

    if not quiet:
        click.echo(click.style("\tCreating '/MergedSnapshots'...", fg="blue"))
        pbar = click.progressbar
    else:
        pbar = dummy_progressbar

    with pbar(range(nelem), length=nelem, label="\t  ") \
            as indices:
        for elem_id in indices:
            # Get the old and new indices.
            old_index = elem_id
            new_index = np.where(inds == elem_id)[0][0]

            gll_point_ids = sem_mesh[old_index]

            # Load displacement from all GLL points.
            for i, var in enumerate(meshes):
                # The list of ids we have is unique but not sorted.
                ids = gll_point_ids.flatten()
                s_ids = np.sort(ids)
                if time_axis == 0:
                    temp = var[:, s_ids]
                    for jpol in range(dim_jpol.size):
                        for ipol in range(dim_ipol.size):
                            idx = ipol * 5 + jpol
                            utemp[i, jpol, ipol, :] = \
                                temp[:, np.argwhere(s_ids == ids[idx])[0][0]]
                else:
                    temp = var[s_ids, :]
                    for jpol in range(dim_jpol.size):
                        for ipol in range(dim_ipol.size):
                            idx = ipol * 5 + jpol
                            utemp[i, jpol, ipol, :] = \
                                temp[np.argwhere(s_ids == ids[idx])[0][0], :]
            x[new_index] = utemp


@click.command()
@click.argument("input_folder", type=click.Path(exists=True, file_okay=False,
                                                dir_okay=True))
@click.argument("output_folder", type=click.Path(exists=False))
@click.option("--contiguous", is_flag=True,
              help="Write a contiguous array - will turn off chunking and "
                   "compression")
@click.option("--compression_level",
              type=click.IntRange(1, 9), default=2,
              help="Compression level from 1 (fast) to 9 (slow).")
@click.option('--method', type=click.Choice(["transpose", "repack", "merge"]),
              required=True,
              help="`transpose` will transpose the data arrays which "
                   "oftentimes results in faster extraction times. `repack` "
                   "will just repack the data and solve some compatibility "
                   "issues. `merge` will create a single much larger file "
                   "which is much quicker to read but will take more space.")
def repack_database(input_folder, output_folder, contiguous,
                    compression_level, method):
    found_filenames = []
    for root, _, filenames in os.walk(input_folder, followlinks=True):
        for filename in sorted(filenames, reverse=True):
            if filename not in ["ordered_output.nc4", "axisem_output.nc4"]:
                continue
            found_filenames.append(os.path.join(root, filename))
            break

    assert found_filenames, "No files named `ordered_output.nc4` found."

    input_folder = os.path.normpath(os.path.realpath(input_folder))
    output_folder = os.path.normpath(os.path.realpath(output_folder))

    if input_folder == output_folder:
        if "ordered_output.nc4" in [os.path.basename(_i) for _i in
                                    found_filenames]:
            raise ValueError("ordered_output.nc4 already exists.")
    else:
        os.makedirs(output_folder)

    if method in ["transpose", "repack"]:
        for _i, filename in enumerate(found_filenames):
            click.echo(click.style(
                "--> Processing file %i of %i: %s" %
                (_i + 1, len(found_filenames), filename), fg="green"))

            output_filename = os.path.join(
                output_folder,
                os.path.relpath(filename, input_folder))

            output_filename = output_filename.replace(
                "axisem_output.nc4", "ordered_output.nc4")

            if not input_folder == output_folder:
                os.makedirs(os.path.dirname(output_filename))

            if method == "transpose":
                transpose = True
            else:
                transpose = False

            repack_file(input_filename=filename,
                        output_filename=output_filename,
                        contiguous=contiguous,
                        transpose=transpose,
                        compression_level=compression_level)
    elif method == "merge":
        merge_files(filenames=found_filenames, output_folder=output_folder,
                    contiguous=contiguous, compression_level=compression_level,
                    quiet=False)
    else:
        raise NotImplementedError


if __name__ == "__main__":
    repack_database()
