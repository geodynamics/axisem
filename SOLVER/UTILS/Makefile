F90 = gfortran
F90FLAGS = -g -fbacktrace -fbounds-check
LDFLAGS_NETCDF = -lm -lnetcdf -lnetcdff -L/usr/lib
LDFLAGS = -lm 
INCLUDE = -I /usr/include

all: xpost_processing xfind_receiver
all_netcdf: xpost_processing_netcdf xfind_receiver_netcdf

xpost_processing: post_processing.o nc_postroutines_dead.o
	$(F90) $(F90FLAGS) $(LDFLAGS) post_processing.o nc_postroutines_dead.o -o xpost_processing

xfind_receiver: find_receiver.o nc_postroutines_dead.o
	$(F90) $(F90FLAGS) $(LDFLAGS) find_receiver.o nc_postroutines_dead.o -o xfind_receiver

xpost_processing_netcdf: post_processing.o nc_postroutines.o
	$(F90) $(F90FLAGS) $(LDFLAGS_NETCDF) post_processing.o nc_postroutines.o -o xpost_processing

xfind_receiver_netcdf: find_receiver.o
	$(F90) $(F90FLAGS) $(LDFLAGS_NETCDF) find_receiver.o nc_postroutines.o -o xfind_receiver

find_receiver.o: find_receiver.f90
	$(F90) -c $(F90FLAGS) find_receiver.f90

post_processing.o: post_processing.f90
	$(F90) -c $(F90FLAGS) post_processing.f90

nc_postroutines.o: 
	$(F90) -c $(F90FLAGS) $(INCLUDE) nc_postroutines.f90

nc_postroutines_dead.o: 
	$(F90) -c $(F90FLAGS) nc_postroutines_dead.f90

clean: 
	rm -rf *.o xpost_processing xfind_receiver

