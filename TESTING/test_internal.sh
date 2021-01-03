#!/bin/bash
set -ev

./TESTING/test_mesher.py

cp $TEST_DIR/inparam_basic_anel .
./submit.py test_prem_iso_anel prem_iso 80.0 --job_type local --ntheta 4

cp $TEST_DIR/inparam_basic_anel .
./submit.py test_prem_ani_anel prem_ani 80.0 --job_type local --ntheta 4

cp $TEST_DIR/inparam_basic_el .
./submit.py test_prem_iso_el prem_iso 80.0 --job_type local --ntheta 4

cp $TEST_DIR/inparam_basic_el .
./submit.py test_prem_ani_el prem_ani 80.0 --job_type local --ntheta 4

cp $TEST_DIR/inparam_basic_el .
./submit.py test_iasp91 iasp91 80.0 --job_type local --ntheta 4

cp $TEST_DIR/inparam_basic_anel .
./submit.py test_ak135f_anel ak135f 80.0 --job_type local --ntheta 4

cp $TEST_DIR/inparam_basic_el .
./submit.py test_ak135f_el ak135f 80.0 --job_type local --ntheta 4


