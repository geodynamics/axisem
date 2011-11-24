psxy -JX14/24 -Cplt-valence.cpt -R-0.2/1.2/-1.2/1.2 -Sc0.05 valence_glob.dat -P \
-B::/:::.'Valence ':WeSn -U/0/-1/ALX"`pwd`" -K -V >valence_glob.ps
psscale -Cplt-valence.cpt -D7/2/6/.5h.5h.5h.5h.5h -O -V >>valence_glob.ps

psxy -JX14/24 -Cplt-valence.cpt -R-0.2/1.2/-1.2/1.2 -Sc0.05 valence_flob.dat -P \
-B::/:::.'Valence ':WeSn -U/0/-1/ALX"`pwd`" -K -V >valence_flob.ps
psscale -Cplt-valence.cpt -D7/2/6/.5h.5h.5h.5h.5h -O -V >>valence_flob.ps

psxy -JX14/24 -Cplt-valence.cpt -R-0.2/1.2/-1.2/1.2 -Sc0.05 valence_slob.dat -P \
-B::/:::.'Valence ':WeSn -U/0/-1/ALX"`pwd`" -K -V >valence_slob.ps
psscale -Cplt-valence.cpt -D7/2/6/.5h.5h.5h.5h.5h -O -V >>valence_slob.ps

