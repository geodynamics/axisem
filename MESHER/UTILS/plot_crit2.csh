psxy -JX14/24 -Cplt.cpt -R-0.2/1.2/-1.2/1.2 -Sc0.04 fort.344 -P \
-B::/:::.'Crit ':WeSn -U/0/-1/ALX"`pwd`" -K -V >test.ps
psscale -Cplt.cpt -D7/2/6/.5h.5h.5h.5h.5h -O -V >>test.ps
gv test.ps
