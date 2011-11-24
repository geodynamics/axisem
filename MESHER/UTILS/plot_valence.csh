psxy -JX14/24 -Cplt-valence.cpt -R-0.2/1.2/-1.2/1.2 -Sc0.05 fort.69 -P \
-B::/:::.'Valence ':WeSn -U/0/-1/ALX"`pwd`" -K -V >valence.ps
psscale -Cplt-valence.cpt -D7/2/6/.5h.5h.5h.5h.5h -O -V >>valence.ps
gs valence.ps
