#!/bin/bash
set -ev

cp $TEST_DIR/TEST01_elastic_isotropic/inparam_basic .
./submit.py TEST01 $TEST_DIR/TEST01_elastic_isotropic/model.bm 80.0 --job_type local --ntheta 4 

cp $TEST_DIR/TEST02_elastic_anisotropic/inparam_basic .
./submit.py TEST02 $TEST_DIR/TEST02_elastic_anisotropic/model.bm 80.0 --job_type local --ntheta 4 

cp $TEST_DIR/TEST03_anelastic_isotropic/inparam_basic .
./submit.py TEST03 $TEST_DIR/TEST03_anelastic_isotropic/model.bm 80.0 --job_type local --ntheta 4 

cp $TEST_DIR/TEST04_anelastic_anisotropic/inparam_basic .
./submit.py TEST04 $TEST_DIR/TEST04_anelastic_anisotropic/model.bm 80.0 --job_type local --ntheta 4 
