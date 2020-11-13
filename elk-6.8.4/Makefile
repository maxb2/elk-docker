
MAKE = make

include make.inc

all:
	cd src; $(MAKE) all
	cd src/eos; $(MAKE)
	cd src/spacegroup; $(MAKE)

clean:
	cd src; $(MAKE) cleanall
	cd src/eos; $(MAKE) clean
	cd src/spacegroup; $(MAKE) clean
	rm -f *.o *.mod *~ fort.* ifc* *.gcno *.exe

test:
	cd tests; ./test.sh

test-libxc:
	cd tests-libxc; ./test.sh

test-mpi:
	cd tests; ./test-mpi.sh

test-all:
	$(MAKE) test
	$(MAKE) test-libxc
	$(MAKE) test-mpi

vim:
	cd src; ./vimelk

