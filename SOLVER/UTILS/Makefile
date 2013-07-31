include ../../make_axisem.macros
ifeq ($(strip $(USE_NETCDF)),true)
   FFLAGS += -Dunc
   LIBS = -L $(NETCDF_PATH)/lib -lnetcdff -Wl,-rpath,$(NETCDF_PATH)/lib
   INCLUDE = -I $(NETCDF_PATH)/include
else
   LIBS = 
   INCLUDE =
endif

all: xpost_processing 

xpost_processing: post_processing.o nc_postroutines.o
	$(FC) $(LDFLAGS) post_processing.o nc_postroutines.o -o xpost_processing $(LIBS) 

xfind_receiver: find_receiver.o nc_postroutines_dead.o
	$(FC) $(FFLAGS) $(LDFLAGS) find_receiver.o nc_postroutines_.o -o xfind_receiver

find_receiver.o: find_receiver.f90
	$(FC) -c $(FFLAGS) find_receiver.f90

post_processing.o: post_processing.f90 nc_postroutines.o
	$(FC) -c $(FFLAGS) $(INCLUDE) post_processing.f90

nc_postroutines.o: nc_postroutines.F90
	$(FC) -c $(FFLAGS) $(INCLUDE) nc_postroutines.F90


clean: 
	rm -rf *.o *.mod xpost_processing xfind_receiver

