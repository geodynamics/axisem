#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  test_axisem.py
#   Purpose:   run the AXISEM tests automatically
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#-------------------------------------------------------------------

import subprocess
import sys
import os
import argparse

parser = argparse.ArgumentParser(description=
    'Run some automated tests for the axisem solver.')

parser.add_argument('-a', help='run all tests', action='store_const',
                    const=True, default=False)

args = parser.parse_args()

if not args.a:
    print '====================================================='
    test_no = \
        raw_input('Please provide the number of desired tests: \n' + \
            '1. test01: explosion\n' + \
            '2. test02: dipole (mxz)\n' + \
            '3. test03: quadpole (mxz)\n' + \
            '4. test04: CMT (source: north pole)\n' + \
            '5. test05: CMT (source: 70-50)\n' + \
            '\n(format = 01,02 OR 1,2,3)' + \
            '\n')
else:
    test_no = '1,2,3,4,5'

print '====================================================='
print 'Requested Test numbers:'

for i in range(0, len(test_no.split(','))):
    num = test_no.split(',')[i]
    if len(num) == 1:
        num = '0' + num
    print 'test' + num
print '====================================================='

for i in range(0, len(test_no.split(','))):
    num = test_no.split(',')[i]
    if len(num) == 1:
        num = '0' + num
    
    print '======================='
    print 'test' + num + ' starts!'
    print '======================='
    
    address = os.path.join('.', 'automated', 'test_' + num)
    output = subprocess.check_call(['python', 'PyAxi.py', address])
    if output != 0: print output_print
    print "=============================================="
