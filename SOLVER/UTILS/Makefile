include ../../make_axisem.macros
ifeq ($(strip $(USE_NETCDF)),true)
   FFLAGS += -Dunc
   LIBS = -L $(strip $(NETCDF_PATH))/lib -lnetcdff -Wl,-rpath,$(strip $(NETCDF_PATH))/lib
   INCLUDE = -I $(strip $(NETCDF_PATH))/include
else
   LIBS = 
   INCLUDE =
endif

all: xpost_processing 

xpost_processing: post_processing.o nc_postroutines.o
	$(FC) $(LDFLAGS) post_processing.o nc_postroutines.o -o xpost_processing $(LIBS) 

post_processing.o: post_processing.f90 nc_postroutines.o ../../make_axisem.macros ../mesh_params.h
	$(FC) -c $(FFLAGS) $(INCLUDE) post_processing.f90 

nc_postroutines.o: nc_postroutines.F90 ../../make_axisem.macros
	$(FC) -c $(FFLAGS) $(INCLUDE) nc_postroutines.F90


clean: 
	rm -rf *.o *.mod xpost_processing xfind_receiver

