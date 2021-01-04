#!/usr/bin/env python
# coding: utf-8
"""
Test all internal meshes of AxiSEM
Author: Simon Stähler, ETH Zürich
"""

import os
import glob
from submit import int_models, run_axisem

def test_mesher_internal(modelname):
    job_name = 'testmesher_' + modelname

    run_axisem(job_name=job_name,
               mesh_file=modelname,
               mesh_period=20.,
               mesher_only=True,
               ntheta=16)

def test_mesher_external(filename):
    job_name = 'testmesher_' + os.path.split(filename)[-1][:-3]

    run_axisem(job_name=job_name,
               mesh_file=filename,
               mesh_period=20.,
               mesher_only=True,
               ntheta=16)



if __name__=='__main__':
    ext_models = glob.glob('MESHER/Models/*.bm')
    ext_models.sort()
    for model in ext_models:
        test_mesher_external(filename=os.path.abspath(model))

    for model in int_models:
        test_mesher_internal(model)

