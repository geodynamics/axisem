#!/usr/bin/env python
# coding: utf-8
"""
Test all internal meshes of AxiSEM
Author: Simon Stähler, ETH Zürich
"""


from submit import int_models, run_axisem

def test_mesher(modelname):
    job_name = 'testmesher' + modelname

    run_axisem(job_name=job_name,
               mesh_file=modelname,
               mesh_period=20.,
               mesher_only=True,
               ntheta=16)



if __name__=='__main__':
    for model in int_models:
        test_mesher(model)