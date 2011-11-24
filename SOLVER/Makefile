PROG =xsem

SRCS =	analyt_homog_radiation.f90 analytic_mapping.f90 \
	analytic_semi_mapping.f90 analytic_spheroid_mapping.f90 \
	apply_masks.f90 background_models.f90 clocks.f90 commpi.f90 \
	commun.f90 data_comm.f90 data_dipole.f90 data_heterogeneous.f90 \
	data_io.f90 data_matr.f90 data_mesh.f90 data_mesh_preloop.f90 \
	data_monopole.f90 data_numbering.f90 data_pointwise.f90 data_proc.f90 \
	data_quadrupole.f90 data_source.f90 data_spec.f90 data_time.f90 \
	def_grid.f90 def_precomp_terms.f90 geom_transf.f90 get_mesh.f90 \
	get_model.f90 global_parameters.f90 lateral_heterogeneities.f90 \
	main.f90 meshes_io.f90 parameters.f90 pointwise_derivatives.f90 \
	rotations.f90 seismograms.f90 source.f90 splib.f90 stiffness.f90 \
	subpar_mapping.f90 time_evol_wave.f90 unit_stride_colloc.f90 \
	unrolled_loops.f90 utlity.f90 wavefields_io.f90 dfour1.f

OBJS =	analyt_homog_radiation.o analytic_mapping.o analytic_semi_mapping.o \
	analytic_spheroid_mapping.o apply_masks.o background_models.o \
	clocks.o commpi.o commun.o data_comm.o data_dipole.o \
	data_heterogeneous.o data_io.o data_matr.o data_mesh.o \
	data_mesh_preloop.o data_monopole.o data_numbering.o data_pointwise.o \
	data_proc.o data_quadrupole.o data_source.o data_spec.o data_time.o \
	def_grid.o def_precomp_terms.o geom_transf.o get_mesh.o get_model.o \
	global_parameters.o lateral_heterogeneities.o main.o meshes_io.o \
	parameters.o pointwise_derivatives.o rotations.o seismograms.o \
	source.o splib.o stiffness.o subpar_mapping.o time_evol_wave.o \
	unit_stride_colloc.o unrolled_loops.o utlity.o wavefields_io.o \
	dfour1.o

LIBS = -lm 
CC = gcc
CFLAGS = -O3 -DF_UNDERSCORE
F90 = mpif90 -Warray-temporaries -fcheck-array-temporaries -fbounds-check -frange-check -pedantic 
FC = gfortran -Warray-temporaries -fcheck-array-temporaries -fbounds-check -frange-check -pedantic 


# cancel m2c implicit rule 
%.o : %.mod 
 

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.M *.mod *.d *.il core 

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

analyt_homog_radiation.o: data_io.o data_mesh.o data_mesh_preloop.o \
	data_proc.o data_source.o data_time.o global_parameters.o utlity.o \
	mesh_params.h
analytic_mapping.o: analytic_semi_mapping.o analytic_spheroid_mapping.o \
	data_mesh.o data_mesh_preloop.o global_parameters.o subpar_mapping.o
analytic_semi_mapping.o: global_parameters.o
analytic_spheroid_mapping.o: data_mesh.o global_parameters.o
apply_masks.o: global_parameters.o mesh_params.h
commpi.o: data_comm.o data_mesh.o data_numbering.o data_proc.o \
	global_parameters.o mesh_params.h mpif.h
commun.o: clocks.o commpi.o data_io.o data_matr.o data_mesh.o \
	data_numbering.o data_proc.o data_time.o global_parameters.o \
	meshes_io.o mesh_params.h
data_comm.o: global_parameters.o
data_dipole.o: global_parameters.o mesh_params.h
data_heterogeneous.o: global_parameters.o mesh_params.h
data_matr.o: global_parameters.o mesh_params.h
data_mesh.o: global_parameters.o mesh_params.h
data_mesh_preloop.o: global_parameters.o mesh_params.h
data_monopole.o: global_parameters.o mesh_params.h
data_numbering.o: mesh_params.h
data_pointwise.o: global_parameters.o mesh_params.h
data_quadrupole.o: global_parameters.o mesh_params.h
data_source.o: global_parameters.o mesh_params.h
data_spec.o: global_parameters.o mesh_params.h
data_time.o: global_parameters.o
def_grid.o: apply_masks.o commun.o data_comm.o data_io.o data_mesh.o \
	data_mesh_preloop.o data_numbering.o data_pointwise.o data_proc.o \
	data_spec.o geom_transf.o get_mesh.o global_parameters.o meshes_io.o \
	splib.o utlity.o mesh_params.h
def_precomp_terms.o: analytic_mapping.o background_models.o commun.o \
	data_dipole.o data_io.o data_matr.o data_mesh.o data_mesh_preloop.o \
	data_monopole.o data_pointwise.o data_proc.o data_quadrupole.o \
	data_source.o data_spec.o def_grid.o geom_transf.o get_mesh.o \
	get_model.o global_parameters.o pointwise_derivatives.o splib.o \
	utlity.o mesh_params.h
geom_transf.o: analytic_mapping.o data_mesh.o data_mesh_preloop.o \
	subpar_mapping.o
get_mesh.o: commun.o data_comm.o data_mesh.o data_mesh_preloop.o \
	data_numbering.o data_proc.o data_time.o global_parameters.o
get_model.o: background_models.o commun.o data_io.o data_mesh.o \
	data_mesh_preloop.o data_proc.o data_time.o def_grid.o \
	global_parameters.o lateral_heterogeneities.o utlity.o mesh_params.h
lateral_heterogeneities.o: data_heterogeneous.o data_io.o data_proc.o \
	data_source.o global_parameters.o utlity.o mesh_params.h
main.o: clocks.o commun.o data_io.o data_mesh.o data_proc.o data_source.o \
	data_time.o def_grid.o get_mesh.o global_parameters.o parameters.o \
	time_evol_wave.o
meshes_io.o: data_io.o data_mesh.o data_mesh_preloop.o data_pointwise.o \
	data_proc.o data_spec.o global_parameters.o utlity.o mesh_params.h
parameters.o: commun.o data_comm.o data_io.o data_mesh.o data_mesh_preloop.o \
	data_numbering.o data_proc.o data_source.o data_time.o \
	global_parameters.o utlity.o mesh_params.h
pointwise_derivatives.o: data_mesh.o data_pointwise.o data_spec.o \
	global_parameters.o unit_stride_colloc.o unrolled_loops.o \
	mesh_params.h
rotations.o: data_io.o data_mesh.o data_proc.o data_source.o \
	global_parameters.o
seismograms.o: commun.o data_io.o data_mesh.o data_mesh_preloop.o \
	data_pointwise.o data_proc.o data_source.o data_time.o \
	global_parameters.o pointwise_derivatives.o rotations.o \
	unit_stride_colloc.o utlity.o wavefields_io.o mesh_params.h
source.o: apply_masks.o commun.o data_io.o data_mesh.o data_mesh_preloop.o \
	data_proc.o data_source.o data_spec.o data_time.o global_parameters.o \
	pointwise_derivatives.o utlity.o
splib.o: global_parameters.o
stiffness.o: data_dipole.o data_matr.o data_mesh.o data_monopole.o \
	data_quadrupole.o data_source.o data_spec.o global_parameters.o \
	unit_stride_colloc.o unrolled_loops.o mesh_params.h
subpar_mapping.o: global_parameters.o
time_evol_wave.o: apply_masks.o clocks.o commun.o data_io.o data_matr.o \
	data_mesh.o data_pointwise.o data_proc.o data_source.o data_time.o \
	def_precomp_terms.o global_parameters.o meshes_io.o parameters.o \
	pointwise_derivatives.o rotations.o seismograms.o source.o \
	stiffness.o unit_stride_colloc.o wavefields_io.o mesh_params.h
unit_stride_colloc.o: global_parameters.o
unrolled_loops.o: global_parameters.o mesh_params.h
utlity.o: data_mesh.o data_mesh_preloop.o data_spec.o geom_transf.o \
	global_parameters.o
wavefields_io.o: data_io.o data_mesh.o data_pointwise.o data_proc.o \
	data_source.o global_parameters.o pointwise_derivatives.o \
	unit_stride_colloc.o mesh_params.h
