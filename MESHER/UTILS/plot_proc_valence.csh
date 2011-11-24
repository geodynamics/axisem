#!/bin/csh -f

set list_grid_files = `ls valence_glob_per_proc.dat0*`

foreach gridfile (${list_grid_files})

set procnum = `ls $gridfile | sed 's/dat/ /' | awk '{print $2}' `

psxy -JX14/24 -Cplt-valence.cpt -R-0.2/1.2/-1.2/1.2 -Sc0.05 $gridfile -P \
-B::/:::.'Valence':WeSn -U/0/-1/ALX"`pwd`" -K -V >valence_proc$procnum.ps
psscale -Cplt-valence.cpt -D7/2/6/.5h.5h.5h.5h.5h -O -V >> valence_proc$procnum.ps

echo wrote processor valence into valence_proc$procnum.ps

end

########### ONLY PLOT VALENCE FOR THE CENTRAL REGION R<0.2
set list_grid_files2 = `ls valence_glob_per_proc_central.dat0*`

foreach gridfile (${list_grid_files2})

set procnum = `ls $gridfile | sed 's/dat/ /' | awk '{print $2}' `

psxy -JX14/24 -Cplt-valence.cpt -R-0.2/1.2/-1.2/1.2 -Sc0.05 $gridfile -P \
-B::/:::.'valence r<0.2':WeSn -U/0/-1/ALX"`pwd`" -K -V >valence_central_proc$procnum.ps
psscale -Cplt-valence.cpt -D7/2/6/.5h.5h.5h.5h.5h -O -V >> valence_central_proc$procnum.ps

echo wrote processor valence into valence_central_proc$procnum.ps

end
