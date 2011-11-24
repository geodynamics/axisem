#!/bin/csh -f

################ SEND ################
set list_grid_files = `ls proc_sendtoproc.dat0*`

foreach gridfile (${list_grid_files})

set procnum = `ls $gridfile | sed 's/dat/ /' | awk '{print $2}' `

psxy -JX14/24 -Cplt-valence.cpt -R-0.2/1.2/-1.2/1.2 -Sc0.05 $gridfile -P \
-B::/:::.'sending to':WeSn -U/0/-1/ALX"`pwd`" -K -V >message_send_proc$procnum.ps
psscale -Cplt-valence.cpt -D7/2/6/.5h.5h.5h.5h.5h -O -V >> message_send_proc$procnum.ps

echo wrote processor message sending into message_send_proc$procnum.ps

end

################ RECEIVE ################
set list_grid_files = `ls proc_recvfromproc.dat0*`

foreach gridfile (${list_grid_files})

set procnum = `ls $gridfile | sed 's/dat/ /' | awk '{print $2}' `

psxy -JX14/24 -Cplt-valence.cpt -R-0.2/1.2/-1.2/1.2 -Sc0.05 $gridfile -P \
-B::/:::.'receiving from':WeSn -U/0/-1/ALX"`pwd`" -K -V >message_recv_proc$procnum.ps
psscale -Cplt-valence.cpt -D7/2/6/.5h.5h.5h.5h.5h -O -V >> message_recv_proc$procnum.ps

echo wrote processor message sending into message_recv_proc$procnum.ps

end
