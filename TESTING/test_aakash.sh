#!/bin/bash
set -ev

cp $TEST_DIR/inparam_basic_anel ./inparam_basic
./submit.py TESTAA01 $TEST_DIR/TESTAA01/prem_ani_e.txt 80.0 --job_type local --ntheta 4 

cp $TEST_DIR/inparam_basic_anel ./inparam_basic
./submit.py TESTAA02 $TEST_DIR/TESTAA02/prem_iso_e.txt 80.0 --job_type local --ntheta 4 

cp $TEST_DIR/inparam_basic_anel ./inparam_basic
./submit.py TESTAA03 $TEST_DIR/TESTAA03/ak135f_scak_1.txt 80.0 --job_type local --ntheta 4 

cp $TEST_DIR/inparam_basic_el ./inparam_basic
./submit.py TESTAA04 $TEST_DIR/TESTAA04/ak135f_scak_2.txt 80.0 --job_type local --ntheta 4 

cp $TEST_DIR/inparam_basic_anel ./inparam_basic
./submit.py TESTAA05 $TEST_DIR/TESTAA05/ak135f_scak_ani_1.txt 80.0 --job_type local --ntheta 4 

cp $TEST_DIR/inparam_basic_el ./inparam_basic
./submit.py TESTAA06 $TEST_DIR/TESTAA06/ak135f_scak_ani_2.txt 80.0 --job_type local --ntheta 4 
