#!/usr/bin/env python
# coding: utf-8

import os
import subprocess as sp
import shutil
import stat
import argparse
import glob


def create_inparam_mesh(mesh_file, mesh_period, nrad, ncl, ntheta=0,
                        max_depth=None, max_colat=None):
    with open('inparam_mesh', 'w') as fid:
        if mesh_file in \
             ['prem_iso', 'prem_iso_solid', 'prem_iso_onecrust',
              'prem_iso_light', 'prem_iso_solid_light',
              'prem_ani', 'prem_ani_onecrust', 'prem_ani_light',
              'ak135', 'ak135f', 'iasp91']:

            fid.write('BACKGROUND_MODEL %s\n' % mesh_file)
            fid.write('EXT_MODEL none \n')

        else:
            fid.write('BACKGROUND_MODEL external \n')
            fid.write('EXT_MODEL "%s" \n' % os.path.abspath(mesh_file))

        fid.write('DOMINANT_PERIOD %f \n' % mesh_period)
        fid.write('NRADIAL_SLICES %d \n' % nrad)

        if max_depth:
            fid.write('MAX_DEPTH %f \n' % max_depth)

        if ncl:
            fid.write('COARSENING_LAYERS %d \n' % ncl)

        if max_colat:
            fid.write('LOCAL_MAX_COLAT %f \n' % max_colat)
        
        if ntheta == 0:
            fid.write('NTHETA_SLICES %d \n' % 1)
            fid.write('ONLY_SUGGEST_NTHETA TRUE \n')
        else:
            fid.write('NTHETA_SLICES %d \n' % ntheta)


def progress():
    pct = 0
    print('[---------------------------------------------------]')
    while True:
        ts = glob.glob('timestamp*.txt')
        for i in range(pct, len(ts)):
            print('=', end='')
        pct = len(ts)
        if pct == 100:
            return


def get_ntheta(mesh_file, mesh_period, ncl,
               max_depth=None, max_colat=None):
    create_inparam_mesh(mesh_file, mesh_period, ntheta=0, nrad=1,
                        max_depth=max_depth, max_colat=max_colat, ncl=ncl)
    output = sp.check_output('./xmesh')
    output_lines = output.split(sep=b'\n')
    for iline in range(0, len(output_lines)):
        line = output_lines[iline]
        if line[0:8] == b'    sugg':
            return int(output_lines[iline + 1])


def hour2hms(hours):
    hour = int(hours)
    hours -= hour
    hours *= 60.
    minute = int(hours)
    hours -= minute
    hours *= 60.
    sec = int(hours)
    return '%02d:%02d:%02d' % (hour, minute, sec)


def define_arguments():
    helptext = 'Create AxiSEM run and submit job.'
    formatter_class = argparse.RawTextHelpFormatter
    parser = argparse.ArgumentParser(description=helptext,
                                     formatter_class=formatter_class)

    helptext = "Job directory name. \n"
    parser.add_argument('job_name', help=helptext)

    helptext = 'Mesh file \n'
    parser.add_argument('mesh_file', help=helptext)

    helptext = 'Mesh period \n'
    parser.add_argument('mesh_period', type=float, help=helptext)

    helptext = 'Number of radial slices \n'
    parser.add_argument('--nrad', type=int, help=helptext,
                        default=1)

    helptext = 'Number of theta slices \n'
    parser.add_argument('--ntheta', type=int, help=helptext,
                        default=2)

    helptext = 'Number of coarsening layers\n'
    parser.add_argument('--ncl', type=int, help=helptext,
                        default=1)

    helptext = 'Maximum colatitude\n'
    parser.add_argument('--max_colat', type=float,
                        help=helptext)

    helptext = 'Maximum depth in kilometer\n'
    parser.add_argument('--max_depth', type=float,
                        help=helptext)

    helptext = 'Source depth in kilometer\n'
    parser.add_argument('--src_depth', type=float, default=0.0,
                        help=helptext)

    helptext = 'Number of theta slices \n'
    parser.add_argument('--ntheta', type=int, help=helptext)

    helptext = 'Number of coarsening layers\n'
    parser.add_argument('--ncl', type=int, default=1, help=helptext)

    helptext = 'Maximum colatitude\n'
    parser.add_argument('--max_colat', type=float, 
                        help=helptext)
    helptext = 'Maximum depth in kilometer\n'
    parser.add_argument('--max_depth', type=float, 
                        help=helptext)

    helptext = 'Source depth in kilometer\n'
    parser.add_argument('--src_depth', type=float, default=0.0,
                        help=helptext)
    
    
    helptext = 'Wall time for the solver in hours\n'
    parser.add_argument('-w', '--walltime', type=float, default=1.0,
                        help=helptext)

    helptext = 'Mail adress for notifications\n'
    parser.add_argument('-m', '--mail_adress', type=str,
                        default='info@nosuchserver.com',
                        help=helptext)

    helptext = 'Daint project account\n'
    parser.add_argument('-a', '--account', type=str,
                        default='ACCOUNT',
                        help=helptext)

    helptext = 'Job type (local, CSCS Daint)\n'
    parser.add_argument('-j', '--jobtype', type=str,
                        default='local',
                        help=helptext)

    return parser

if __name__ == "__main__":
    parser = define_arguments()
    args = parser.parse_args()

    jobname = args.job_name
    base_dir = os.getcwd()
    os.chdir(base_dir)

    # Create rundir
    rundir = os.path.abspath(os.path.join('runs', jobname))
    os.mkdir(rundir)

    # Create directory for the mesh
    print('[MESHER]')
    os.chdir(base_dir)
    print('  Creating mesh directory')
    meshdir = os.path.abspath(os.path.join(rundir, 'Mesh'))
    os.mkdir(meshdir)
    os.mkdir(os.path.join(meshdir, 'Diags'))

    mesh_exe_path = os.path.join(meshdir, 'xmesh')
    shutil.copyfile(src=os.path.join('MESHER', 'xmesh'),
                    dst=mesh_exe_path)
    st = os.stat(mesh_exe_path)
    os.chmod(mesh_exe_path, st.st_mode | stat.S_IEXEC)

    os.chdir(meshdir)

    if args.mesh_file in [
            'prem_iso', 'prem_iso_solid', 'prem_iso_onecrust',
            'prem_iso_light', 'prem_iso_solid_light',
            'prem_ani', 'prem_ani_onecrust', 'prem_ani_light',
            'ak135', 'ak135f', 'iasp91']:
        print('  Using internal model %s' % args.mesh_file)
        int_model = True
    else:
        int_model = False
        fnam_mesh_file = os.path.split(args.mesh_file)[-1]
        shutil.copyfile(src=args.mesh_file,
                        dst=os.path.join(meshdir, fnam_mesh_file))
        print('  Using external model %s' % fnam_mesh_file)

    if args.ntheta:
        ntheta = args.ntheta
    else:
        ntheta = get_ntheta(fnam_mesh_file, args.mesh_period,
                            ncl=args.ncl,
                            max_depth=args.max_depth,
                            max_colat=args.max_colat)
        print('  Optimal number of theta slices used: %d' % ntheta)

    nrad = args.nrad
    create_inparam_mesh(args.mesh_file,
                        args.mesh_period,
                        ntheta=ntheta,
                        nrad=nrad,
                        ncl=args.ncl,
                        max_depth=args.max_depth,
                        max_colat=args.max_colat)

    ncpu = ntheta * nrad
    print('  Number of cores used: %4d' % ncpu)

    if args.jobtype == 'local':  # local run, consecutive
        os.chdir(os.path.abspath(meshdir))
        print('  Starting Mesher')
        output = sp.check_output('./xmesh > OUTPUT_MESH', shell=True)

    elif args.jobtype == 'Daint':
        batch_mesher_fmt =                                              \
            '#!/bin/bash -l \n' +                                       \
            '#SBATCH --ntasks=1 \n' +                                   \
            '#SBATCH --ntasks-per-node=1 \n' +                          \
            '#SBATCH --ntasks-per-core=1 \n' +                          \
            '#SBATCH --cpus-per-task=1 \n' +                            \
            '#SBATCH --time=00:30:00 \n' +                              \
            '#SBATCH --account=%s \n' % args.account +                  \
            '#SBATCH --mail-type=BEGIN,FAIL \n' +                       \
            '#SBATCH --mail-user=%s \n' % args.mail_adress +            \
            '#SBATCH --constraint=mc \n' +                              \
            '#SBATCH --partition=prepost \n' +                          \
            '#SBATCH --workdir=%s \n' % meshdir +                       \
            'export OMP_NUM_THREADS=8 \n' +                             \
            'module load slurm \n' +                                    \
            'echo "The current job ID is $SLURM_JOB_ID" \n' +           \
            'echo "Running on $SLURM_JOB_NUM_NODES nodes" \n' +         \
            'echo "Using $SLURM_NTASKS_PER_NODE tasks per node" \n' +   \
            'echo "A total of $SLURM_NTASKS tasks is used" \n' +        \
            'echo "using $SLURM_CPUS_PER_TASK omp threads" \n' +        \
            './xmesh > OUTPUT_MESHER'

        path_sbatch_mesher = os.path.join(rundir, 'job_%s_mesh.sh' % (jobname))
        with open(path_sbatch_mesher, 'w') as fid:
            fid.write(batch_mesher_fmt)

    # SOLVER
    print('[SOLVER]')
    inparam_source = {'PX':
                      'SOURCE_TYPE thetaforce  \n' +
                      'SOURCE_DEPTH %f  \n' +
                      'SOURCE_LAT 90.0  \n' +
                      'SOURCE_LON 0.0  \n' +
                      'SOURCE_AMPLITUDE  1.E20',
                      'PZ':
                      'SOURCE_TYPE vertforce  \n' +
                      'SOURCE_DEPTH %f  \n' +
                      'SOURCE_LAT 90.0  \n' +
                      'SOURCE_LON 0.0  \n' +
                      'SOURCE_AMPLITUDE  1.E20'}

    path_sbatch_solver = dict()

    desc = {'PX': 'horizontal', 'PZ': 'vertical'}

    for part_run in ['PX', 'PZ']:
        print('  Creating solver dir for %s force source' % desc[part_run])
        os.chdir(base_dir)

        # Create directory for the solver (PZ)
        solverdir = os.path.abspath(os.path.join(rundir, part_run))
        os.mkdir(solverdir)

        # Copy Solver executable and make it executable
        solver_exe_path = os.path.join(solverdir, 'axisem')
        shutil.copyfile(src=os.path.join('SOLVER', 'axisem'),
                        dst=solver_exe_path)
        st = os.stat(solver_exe_path)
        os.chmod(solver_exe_path, st.st_mode | stat.S_IEXEC)

        # Copy inparam files
        for file in ['inparam_basic', 'inparam_advanced']:
            shutil.copyfile(src=file,
                            dst=os.path.join(solverdir, file))

        if not int_model:
            shutil.copyfile(src=args.mesh_file,
                            dst=os.path.join(solverdir,
                                             'external_model.bm'))

        with open(os.path.join(solverdir, 'inparam_source'), 'w') as fid:
            fid.write(inparam_source[part_run] % args.src_depth)

        # Create output directories
        os.mkdir(os.path.join(solverdir, 'Data'))
        os.mkdir(os.path.join(solverdir, 'Info'))

        # Create symlink to mesh directory
        os.symlink(os.path.abspath(meshdir), os.path.join(solverdir, 'Mesh'))

        if args.jobtype == 'local':  # local run, consecutive
            os.chdir(os.path.abspath(solverdir))
            print('  Starting %s AxiSEM simulation' % desc[part_run])
            output = sp.check_output('mpirun -n %d ./axisem > OUTPUT' % ncpu,
                                     shell=True)

        elif args.jobtype == 'Daint':
            batch_solver_fmt = \
                '#!/bin/bash -l \n' +                                         \
                '#SBATCH --ntasks=%d \n' % ncpu +                             \
                '#SBATCH --ntasks-per-node=36 \n' +                           \
                '#SBATCH --ntasks-per-core=1 \n' +                            \
                '#SBATCH --cpus-per-task=1 \n' +                              \
                '#SBATCH --time=%s \n' % hour2hms(args.walltime) +            \
                '#SBATCH --account=%s \n' % args.account +                    \
                '#SBATCH --mail-type=FAIL \n' +                               \
                '#SBATCH --mail-user=%s \n' % args.mail_adress +              \
                '#SBATCH --constraint=mc \n' +                                \
                '#SBATCH --partition=normal \n' +                             \
                '#SBATCH --workdir=%s \n' % os.path.abspath(solverdir) +      \
                'export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK \n' +            \
                'module load slurm \n' +                                      \
                'echo "The current job ID is $SLURM_JOB_ID" \n' +             \
                'echo "Running on $SLURM_JOB_NUM_NODES nodes" \n' +           \
                'echo "Using $SLURM_NTASKS_PER_NODE tasks per node" \n' +     \
                'echo "A total of $SLURM_NTASKS tasks is used" \n' +          \
                'echo "using $SLURM_CPUS_PER_TASK omp threads" \n' +          \
                'srun --ntasks-per-node $SLURM_NTASKS_PER_NODE ' +            \
                '-n $SLURM_NTASKS ./axisem >& OUTPUT_%s' % part_run

            path_sbatch_solver[part_run] = \
                os.path.join(rundir, 'job_%s_%s.sh' % (jobname, part_run))
            with open(path_sbatch_solver[part_run], 'w') as fid:
                fid.write(batch_solver_fmt)

    # Run field_transform (the Instaseis version)
    repack_path = os.path.join(base_dir, 'SOLVER', 'UTILS', 'repack_db.py')
    repack_call = repack_path + ' --method transpose ' + '. ' \
        + jobname+'_packed ' + '> OUTPUT_FT'

    print('[REPACK]')
    if args.jobtype == 'local':  # local run, consecutive
        os.chdir(os.path.abspath(rundir))
        print('  Starting database repack for Instaseis')
        output = sp.check_output(repack_call, shell=True)

    elif args.jobtype == 'Daint':
        batch_FT_fmt =                                                  \
            '#!/bin/bash -l \n' +                                       \
            '#SBATCH --ntasks=1 \n' +                                   \
            '#SBATCH --ntasks-per-node=1 \n' +                          \
            '#SBATCH --ntasks-per-core=1 \n' +                          \
            '#SBATCH --cpus-per-task=1 \n' +                            \
            '#SBATCH --time=24:00:00 \n' +                              \
            '#SBATCH --mem=120GB \n' +                                  \
            '#SBATCH --account=%s \n' % args.account +                  \
            '#SBATCH --mail-type=END,FAIL \n' +                         \
            '#SBATCH --mail-user=%s \n' % args.mail_adress +            \
            '#SBATCH --constraint=mc \n' +                              \
            '#SBATCH --partition=normal \n' +                           \
            '#SBATCH --workdir=%s \n' % rundir +                        \
            'export OMP_NUM_THREADS=1 \n' +                             \
            'module load slurm \n' +                                    \
            'echo "The current job ID is $SLURM_JOB_ID" \n' +           \
            'echo "Running on $SLURM_JOB_NUM_NODES nodes" \n' +         \
            'echo "Using $SLURM_NTASKS_PER_NODE tasks per node" \n' +   \
            'echo "A total of $SLURM_NTASKS tasks is used" \n' +        \
            'echo "using $SLURM_CPUS_PER_TASK omp threads" \n' +        \
            repack_call

    # Submit the jobs

    if args.jobtype == 'Daint':  # Daint (CSCS), Slurm-based
        path_sbatch_FT = os.path.join(rundir, 'job_%s_FT.sh' % (jobname))
        with open(path_sbatch_FT, 'w') as fid:
            fid.write(batch_FT_fmt)

        res_submit_mesh = sp.check_output('sbatch ' + path_sbatch_mesher,
                                          shell=True)
        jobid_mesher = int(res_submit_mesh.split()[3])
        print('Mesher JOBID: ', jobid_mesher)

        res_submit_solver = dict()
        jobid_solver = []
        for part_run in ['PX', 'PZ']:
            res_submit_solver = sp.check_output(
                     'sbatch --dependency=afterok:%d ' % jobid_mesher +
                     path_sbatch_solver[part_run],
                     shell=True)
            jobid_solver.append(int(res_submit_solver.split()[3]))
            print('Solver %s JOBID: ' % part_run, jobid_solver[-1])

        res_submit_FT = sp.check_output(
                'sbatch --dependency=afterok:%d:%d %s' %
                (jobid_solver[0], jobid_solver[1], path_sbatch_FT),
                shell=True)
        jobid_FT = int(res_submit_FT.split()[3])
        print('FT JOBID:     ', jobid_FT)
