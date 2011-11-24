#!/bin/csh 
surface $1 -R0/1/-1/1 -I0.01/0.01 -V -Gtoto.grd
grdcontour -A0.2 -C0.05 toto.grd -JX10/20 -P -R toto.grd -K >test.ps
grdimage -Cplt.cpt toto.grd -JX -R -O >>test.ps
gv test.ps 
