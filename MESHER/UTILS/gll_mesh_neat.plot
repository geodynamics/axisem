set size ratio -1 
#set term x11
set term postscript portrait solid "Helvetica" 24
set output "grid.ps"
set noborder 
set noxtics ; set noytics
#set title "Computational Mesh"
plot "colloc_grid.dat" t'' with l , "mesh.dat" t'' with  lines lw 2
#pause -1 "Hit any key to exit..."
