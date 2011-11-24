#!/bin/csh

# convert from  radian to degrees
set fname = $1
echo "filename is"
echo $fname
set pi = 3.1415927
set dpr = ` echo $pi | awk ' { print 180./$1 } ' `
set fact = 0.01
echo "1 radian is"
echo $dpr
echo "degrees"
cat $fname | awk ' {print 90.-$2*c,$1,$3}' c=$dpr d=$fact |cat > testu1.dat
blockmedian testu1.dat -R-90/90/0./1.  -I1./0.01  -V > nico1.dat
triangulate nico1.dat  -R-90/90/0./1. -I1./.01 -Gu1.grd
#/usr/local/GMT3.4/bin/surface testu1.dat -R-90./90/.5/1.5  -I1./0.01  -V -Gu1.grd
#surface testu1.dat -R-90./90/.0/1.   -I3./0.02  -V -Gu1.grd
#surface testu1.dat -R-90./90/0./1.   -I.2/0.01  -V -Gu1.grd
grdimage u1.grd -R -B:'':/:''::.'':WeSn -Cplt.cpt -JP8 -P -V -Y+4 -X+4 -K > test.ps
grdcontour u1.grd -R  -C.05 -JP8 -W3 -O -L0.01/1. -K >> test.ps
grdcontour u1.grd -R  -C.05 -JP8 -W3ta3 -O -L-1/-.01  >> test.ps
#grdcontour u1.grd -R  -C.05 -JP8 -W3ta3 -O -L-1/-.01 -U/0/-1/ALX"`pwd`"/"Field:$fname " >> test.ps
#grdcontour u1.grd -R  -C.004 -JP8 -W3 -O -U/0/-1/ALX"`pwd`"/"Field:$fname " >> test.ps
gv   test.ps
rm -f *.grd 
