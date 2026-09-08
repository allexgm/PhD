FC      = gfortran
FFLAGS  = -O2 -Wall
MPIFC   = mpif90

all: parallel

#serial: serial.f90
#	$(FC) $(FFLAGS) -o serial serial.f90

# parallel target added later
parallel: parallel.f90
 	$(MPIFC) $(FFLAGS) -o parallel parallel.f90

clean:
	rm -f serial parallel *.dat

.PHONY: all parallel clean
