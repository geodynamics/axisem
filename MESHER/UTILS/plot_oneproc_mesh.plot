set size ratio -1
set term postscript portrait solid "Helvetica" 24
set output "grid_proc0000.ps"
set noborder
set noxtics ; set noytics
plot "serend_coords_per_proc.dat0000" t'' with l lw .4
