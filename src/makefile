FC = gfortran
FFLAGS = -O3 -fdefault-real-8

TARGETS = euler

ALL: $(TARGETS)

SRC = $(wildcard *.f95)
OBJ = $(patsubst %.f95,%.o,$(SRC))

euler: $(OBJ)
	$(FC) -o euler $(OBJ)

comvar.mod: comvar.f95
	$(FC) -c $(FFLAGS) comvar.f95

%.o: %.f95 comvar.mod
	$(FC) -c $(FFLAGS) $*.f95

clean:
	rm -f $(OBJ) $(TARGETS) *.mod
