#!/bin/bash

gfortran -c -I $HOME/local/include/ field_transform.f90
gfortran field_transform.o -o field_transform -L $HOME/local/lib -lnetcdff -Wl,-rpath,$HOME/local/lib
