#/bin/csh -f
#dd plot
psxy -L -M $1fort.56 -R0/1.5/-1.5/1.5  -JX10/20 -V -P > toto.ps
gv toto.ps 
rm toto.ps
