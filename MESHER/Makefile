PROG =xmesh

SRCS =	analytic_mapping.f90 analytic_semi_mapping.f90 \
	analytic_spheroid_mapping.f90 background_models.f90 \
	data_bkgrdmodel.f90 data_coarse.f90 data_diag.f90 data_gllmesh.f90 \
	data_grid.f90 data_mesh.f90 data_numbering.f90 data_pdb.f90 \
	data_spec.f90 discont_meshing.f90 gllmeshgen.f90 \
	global_parameters.f90 input.f90 main.f90 mapping_spheroid.f90 \
	mesh_info.f90 meshgen.f90 model_discontinuities.f90 numbering.f90 \
	parallelization.f90 pdb.f90 splib.f90 subpar_mapping.f90 \
	test_bkgrdmodel.f90

OBJS =	analytic_mapping.o analytic_semi_mapping.o \
	analytic_spheroid_mapping.o background_models.o data_bkgrdmodel.o \
	data_coarse.o data_diag.o data_gllmesh.o data_grid.o data_mesh.o \
	data_numbering.o data_pdb.o data_spec.o discont_meshing.o \
	gllmeshgen.o global_parameters.o input.o main.o mapping_spheroid.o \
	mesh_info.o meshgen.o model_discontinuities.o numbering.o \
	parallelization.o pdb.o splib.o subpar_mapping.o test_bkgrdmodel.o

LIBS = -lm 
CC = gcc
CFLAGS = -O3 -DF_UNDERSCORE
F90 = mpif90 -O3 
FC = gfortran  -O3 


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

analytic_mapping.o: analytic_semi_mapping.o analytic_spheroid_mapping.o \
	data_mesh.o global_parameters.o subpar_mapping.o
analytic_semi_mapping.o: data_mesh.o global_parameters.o
analytic_spheroid_mapping.o: data_mesh.o global_parameters.o
data_mesh.o: global_parameters.o
discont_meshing.o: background_models.o data_bkgrdmodel.o data_diag.o \
	data_spec.o global_parameters.o model_discontinuities.o splib.o
gllmeshgen.o: analytic_mapping.o data_diag.o data_gllmesh.o data_mesh.o \
	data_spec.o splib.o
input.o: data_bkgrdmodel.o data_coarse.o data_diag.o data_grid.o data_mesh.o \
	data_spec.o global_parameters.o
main.o: data_bkgrdmodel.o data_grid.o data_mesh.o discont_meshing.o \
	gllmeshgen.o input.o mesh_info.o meshgen.o numbering.o \
	parallelization.o pdb.o test_bkgrdmodel.o
mesh_info.o: data_bkgrdmodel.o data_diag.o data_gllmesh.o data_grid.o \
	data_mesh.o data_numbering.o data_spec.o global_parameters.o
meshgen.o: analytic_spheroid_mapping.o data_bkgrdmodel.o data_coarse.o \
	data_diag.o data_grid.o data_mesh.o data_numbering.o data_spec.o \
	global_parameters.o numbering.o
model_discontinuities.o: data_bkgrdmodel.o data_grid.o
numbering.o: data_diag.o data_gllmesh.o data_grid.o data_mesh.o \
	data_numbering.o data_spec.o
parallelization.o: data_bkgrdmodel.o data_coarse.o data_diag.o data_gllmesh.o \
	data_grid.o data_mesh.o data_pdb.o data_spec.o test_bkgrdmodel.o
pdb.o: data_bkgrdmodel.o data_diag.o data_gllmesh.o data_grid.o data_mesh.o \
	data_numbering.o data_pdb.o data_spec.o global_parameters.o \
	numbering.o
splib.o: global_parameters.o
subpar_mapping.o: global_parameters.o
test_bkgrdmodel.o: background_models.o data_bkgrdmodel.o data_diag.o \
	data_gllmesh.o data_grid.o data_mesh.o data_spec.o \
	global_parameters.o
