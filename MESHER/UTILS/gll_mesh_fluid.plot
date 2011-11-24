set size ratio -1 
#set term x11
set term postscript portrait solid "Helvetica" 24
set output "grid_fluid.ps"
set noborder 
set noxtics ; set noytics
#set title "Elementary Conforming Mesh"
plot "colloc_grid_fluid.dat" t'' with lines lw .5, "mesh_fluid.dat" t'' with lines lw 1.
#pause -1 "Hit any key to exit..."
