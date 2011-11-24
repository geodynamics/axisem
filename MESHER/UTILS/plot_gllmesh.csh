#/bin/csh -f
cp $1mesh.dat .
cp $1colloc_grid.dat .
gnuplot gll_mesh.plot
ps2pdf -sPAPERSIZE=a4 grid.ps
evince grid.pdf
