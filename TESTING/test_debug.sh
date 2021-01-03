#!/bin/bash
set -ev

cp $TEST_DIR/TEST01_elastic_isotropic/inparam_basic .
./submit.py TEST01 $TEST_DIR/TEST01_elastic_isotropic/model.bm 100.0 --jobtype local --ntheta 4 

cp $TEST_DIR/inparam_basic_el .
./submit.py test_prem_iso_el prem_iso 100.0 --jobtype local --ntheta 4
