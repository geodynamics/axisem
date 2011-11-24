#/bin/csh -f
gnuplot plot_mesh2.plot
gnuplot plot_mesh_central.plot
gnuplot plot_fluid.plot
gnuplot plot_solid.plot
gnuplot plot_uppermantle.plot

echo "see grid.ps, grid_central.ps, grid_fluid.ps, grid_solid.ps, uppermantle_grid.ps"
