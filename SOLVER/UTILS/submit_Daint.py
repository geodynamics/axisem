# coding: utf-8

import os
import subprocess as sp
import shutil
import stat
import argparse
import numpy as np


def create_inparam_mesh(mesh_file, mesh_period, ntheta=0, nrad=1):
    with open('inparam_mesh', 'w') as fid:
        fid.write('BACKGROUND_MODEL external \n')
        fid.write('EXT_MODEL "%s" \n' % os.path.abspath(mesh_file))
        fid.write('DOMINANT_PERIOD %f \n' % mesh_period)
        
        if ntheta==0:
            fid.write('NTHETA_SLICES %d \n' % 1)
            fid.write('NRADIAL_SLICES %d \n' % nrad)
            fid.write('ONLY_SUGGEST_NTHETA TRUE \n')
        else:
            fid.write('NTHETA_SLICES %d \n' % ntheta)
            fid.write('NRADIAL_SLICES %d \n' % nrad)


def get_ntheta(mesh_file, mesh_period):
    create_inparam_mesh(mesh_file, mesh_period, ntheta=0, nrad=1)
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
    parser.add_argument('--nrad', type=int, help=helptext)

    helptext = 'Wall time for the solver in hours\n'
    parser.add_argument('-w', '--walltime', type=float, default=1.0, 
                        help=helptext)

    helptext = 'Mail adress for notifications\n'
    parser.add_argument('-m', '--mail_adress', type=str, 
                        default='info@nosuchserver.com', 
                        help=helptext)

    # This does not work, since you cannot SSH out from compute nodes
    # helptext = 'Transfer adress, where DB should be moved to\n'
    # parser.add_argument('-t', '--transfer_adress', type=str, 
    #                     default='info@nosuchserver.com', 
    #                     help=helptext)

    helptext = 'Daint project account\n'
    parser.add_argument('-a', '--account', type=str, 
                        default='ACCOUNT', 
                        help=helptext)

    return parser

parser = define_arguments()
args = parser.parse_args()

jobname = args.job_name
base_dir = os.getcwd()
os.chdir(base_dir)
# Create rundir
rundir = os.path.abspath(os.path.join('runs', jobname))
os.mkdir(rundir)


# Create directory for the mesh
os.chdir(base_dir)
print('Creating mesh directory')
meshdir = os.path.abspath(os.path.join(rundir, 'Mesh'))
os.mkdir(meshdir)
os.mkdir(os.path.join(meshdir, 'Diags'))

mesh_exe_path = os.path.join(meshdir, 'xmesh')
shutil.copyfile(src=os.path.join('MESHER', 'xmesh'), 
                dst=mesh_exe_path)
st = os.stat(mesh_exe_path)
os.chmod(mesh_exe_path, st.st_mode | stat.S_IEXEC)

os.chdir(meshdir)

fnam_mesh_file = os.path.split(args.mesh_file)[-1]

shutil.copyfile(src=args.mesh_file,
                dst=os.path.join(meshdir, fnam_mesh_file))
ntheta = get_ntheta(fnam_mesh_file, args.mesh_period)
print('  Optimal number of theta slices: %d' % ntheta)
nrad = args.nrad
create_inparam_mesh(args.mesh_file, 
                    args.mesh_period, 
                    ntheta=ntheta, 
                    nrad=nrad)
         
ncpu = ntheta * nrad
print('  Number of cores used:           %d' % ncpu)

batch_mesher_fmt =                                                  \
        '#!/bin/bash -l \n' +                                       \
        '#SBATCH --ntasks=1 \n' +                                   \
        '#SBATCH --ntasks-per-node=1 \n' +                          \
        '#SBATCH --ntasks-per-core=1 \n' +                          \
        '#SBATCH --cpus-per-task=1 \n' +                            \
        '#SBATCH --time=00:30:00 \n' +                              \
        '#SBATCH --account=%s \n' % args.account +                  \
        '#SBATCH --mail-type=ALL \n' +                              \
        '#SBATCH --mail-user=%s \n' % args.mail_adress +            \
        '#SBATCH --constraint=mc \n' +                              \
        '#SBATCH --partition=prepost \n' +                          \
        '#SBATCH --workdir=%s \n' % meshdir +          	            \
        'export OMP_NUM_THREADS=8 \n' +                             \
        'module load slurm \n' +                                    \
        'echo "The current job ID is $SLURM_JOB_ID" \n' +           \
        'echo "Running on $SLURM_JOB_NUM_NODES nodes" \n' +         \
        'echo "Using $SLURM_NTASKS_PER_NODE tasks per node" \n' +   \
        'echo "A total of $SLURM_NTASKS tasks is used" \n' +        \
        'echo "using $SLURM_CPUS_PER_TASK omp threads" \n' +        \
        './xmesh > OUTPUT_MESHER'

path_sbatch_mesher = os.path.join(rundir, 'sbatch_mesher.sh')
with open(path_sbatch_mesher, 'w') as fid:
    fid.write(batch_mesher_fmt)


inparam_source = {'PX': 
                    'SOURCE_TYPE thetaforce  \n' + 
                    'SOURCE_DEPTH 0.0  \n' +
                    'SOURCE_LAT 90.0  \n' +
                    'SOURCE_LON 0.0  \n' +
                    'SOURCE_AMPLITUDE  1.E20', 
                  'PZ': 
                    'SOURCE_TYPE vertforce  \n' +
                    'SOURCE_DEPTH 0.0  \n' +
                    'SOURCE_LAT 90.0  \n' +
                    'SOURCE_LON 0.0  \n' +
                    'SOURCE_AMPLITUDE  1.E20'}

path_sbatch_solver = dict()

for part_run in ['PX', 'PZ']:
    print('Creating solver dir for run %s' % part_run)
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

    shutil.copyfile(src=args.mesh_file,
                    dst=os.path.join(solverdir, 'external_model.bm'))
        
    with open(os.path.join(solverdir, 'inparam_source'), 'w') as fid:
        fid.write(inparam_source[part_run])

    # Create output directories
    os.mkdir(os.path.join(solverdir, 'Data'))
    os.mkdir(os.path.join(solverdir, 'Info'))

    # Create symlink to mesh directory
    os.symlink(os.path.abspath(meshdir), os.path.join(solverdir, 'Mesh'))
        
    batch_solver_fmt = \
            '#!/bin/bash -l \n' +                                         \
            '#SBATCH --ntasks=%d \n' % ncpu +                             \
            '#SBATCH --ntasks-per-node=36 \n' +                           \
            '#SBATCH --ntasks-per-core=1 \n' +                            \
            '#SBATCH --cpus-per-task=1 \n' +                              \
            '#SBATCH --time=%s \n' % hour2hms(args.walltime) +            \
            '#SBATCH --account=%s \n' % args.account +                    \
            '#SBATCH --mail-type=ALL \n' +                                \
            '#SBATCH --mail-user=%s \n' %args.mail_adress +               \
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
            'srun --ntasks-per-node $SLURM_NTASKS_PER_NODE -n $SLURM_NTASKS ./axisem >& OUTPUT_%s' % part_run
    
    path_sbatch_solver[part_run] = os.path.join(rundir, 'sbatch_%s.sh' % part_run)
    with open(path_sbatch_solver[part_run] , 'w') as fid:
        fid.write(batch_solver_fmt)


# Run field_transform (the Instaseis version)

repack_path = os.path.join(base_dir, 'SOLVER', 'UTILS', 'repack_db.py')
repack_call = repack_path + ' --method transpose ' + '. ' + jobname+'_packed ' + '> OUTPUT_FT'

# This does not really work, since you cannot SSH out from compute nodes
#if args.transfer_adress:
#  repack_call += '\n scp -r %s_packed %s' % (jobname, args.transfer_adress)

batch_FT_fmt =                                                      \
        '#!/bin/bash -l \n' +                                       \
        '#SBATCH --ntasks=1 \n' +                                   \
        '#SBATCH --ntasks-per-node=1 \n' +                          \
        '#SBATCH --ntasks-per-core=1 \n' +                          \
        '#SBATCH --cpus-per-task=1 \n' +                            \
        '#SBATCH --time=24:00:00 \n' +                              \
        '#SBATCH --mem=120GB \n' +                                  \
        '#SBATCH --account=%s \n' % args.account +                  \
        '#SBATCH --mail-type=ALL \n' +                              \
        '#SBATCH --mail-user=%s \n' %args.mail_adress +             \
        '#SBATCH --constraint=mc \n' +                              \
        '#SBATCH --partition=normal \n' +                           \
        '#SBATCH --workdir=%s \n' % rundir +          	            \
        'export OMP_NUM_THREADS=1 \n' +                             \
        'module load slurm \n' +                                    \
        'echo "The current job ID is $SLURM_JOB_ID" \n' +           \
        'echo "Running on $SLURM_JOB_NUM_NODES nodes" \n' +         \
        'echo "Using $SLURM_NTASKS_PER_NODE tasks per node" \n' +   \
        'echo "A total of $SLURM_NTASKS tasks is used" \n' +        \
        'echo "using $SLURM_CPUS_PER_TASK omp threads" \n' +        \
        repack_call


# Submit the jobs

path_sbatch_FT = os.path.join(rundir, 'sbatch_FT.sh')
with open(path_sbatch_FT, 'w') as fid:
        fid.write(batch_FT_fmt)    

res_submit_mesh = sp.check_output('sbatch ' + path_sbatch_mesher, shell=True)
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


res_submit_FT = sp.check_output(\
        'sbatch --dependency=afterok:%d:%d ' % (jobid_solver[0], jobid_solver[1]) +
         path_sbatch_FT,
         shell=True)
