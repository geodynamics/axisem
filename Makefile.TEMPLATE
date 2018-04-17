include make_axisem.macros

all: MESHER/xmesh SOLVER/axisem

MESHER/xmesh: MESHER/gitversion.h ./make_axisem.macros
	@cd MESHER; $(MAKE) $(MFLAGS)

SOLVER/axisem: SOLVER/gitversion.h ./make_axisem.macros
	@cd SOLVER; $(MAKE) $(MFLAGS)

clean:
	cd MESHER; rm -f *.o *.M *.mod *.d *.il core *.html *.gcda *.gcno *.h; cd ..
	cd SOLVER; rm -f *.o *.M *.mod *.d *.il core *.html *.gcda *.gcno *.h; cd ..
	rm -f bin/axisem bin/xmesh

SOLVER/gitversion.h: .git/HEAD .git/index
	echo "character(len=*), parameter :: gitversion = \"$(shell git rev-parse HEAD)\"" > $@

MESHER/gitversion.h: .git/HEAD .git/index
	echo "character(len=*), parameter :: gitversion = \"$(shell git rev-parse HEAD)\"" > $@


