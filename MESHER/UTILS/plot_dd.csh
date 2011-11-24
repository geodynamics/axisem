#/bin/csh -f
#dd plot
triangulate $1fort.55 -R0./1./-1./1.  -I0.001/0.001 -V -Gdd.grd
#surface $1fort.55 -R0./1.5/-1.5/1.5  -I0.001/0.001 -V -Gdd.grd
grdimage -Cpltdd.cpt -R0./1./-1./1. -P  -JX10/20 dd.grd  -K > toto.ps
#grdcontour -C0.2 -A1.  -R0./1./0./1   -JX15/15 dd.grd  -O -K >> toto.ps
psxy -M -L $1fort.56 -R -JX -O -V -W1  >> toto.ps
gv toto.ps
#rm toto.ps
