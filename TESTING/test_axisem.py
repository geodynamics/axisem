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

print '====================================================='
test_no = \
    raw_input('Please provide the number of desired tests: (format = 01,02)\n' + \
        '1. test01: explosion\n' + \
        '2. test02: dipole (mxz)\n' + \
        '3. test03: quadpole (mxz)\n' + \
        '4. test04: CMT (source: noth pole)\n' + \
        '5. test05: CMT (source: 70-50)\n' + \
        '\n')
print '====================================================='

m = -1

for i in range(0, len(test_no.split(','))):
    num = test_no.split(',')[i]
    
    inpython_open = open('./inpython.cfg', 'r')
    inpython_read = inpython_open.readlines()
    
    search = 'test_folder = '
    for j in range(0, len(inpython_read)):
        if inpython_read[j].find(search) != -1:
            inpython_read[j] = 'test_folder = ./automated/test_' + num + '\n'
            m = j
    
    if m == -1:
        print "Something is wrong with the inpython.cfg file!"
        print "=============================================="
        sys.exit()
    
    inpython_open.close()
    inpython_open = open('./inpython.cfg', 'w')

    for k in range(0, len(inpython_read)):
        inpython_open.write(inpython_read[k])

    inpython_open.close()
    
    output = subprocess.check_call(['python', 'PyAxi.py'])
    if output != 0: print output_print
    print "=============================================="
