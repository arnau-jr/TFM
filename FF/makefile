objects = analysis.o tools.o ff_utils.o
comp = gfortran
OPT = -O3 -fbounds-check

analysis.x: $(objects)
	$(comp) -o analysis.x $(OPT) $(objects)

tools.f90: ../tools.f90
	cp ../tools.f90 .
	
ff_utils.f90: ../ff_utils.f90
	cp ../ff_utils.f90 .

tools.o: tools.f90
	$(comp) -c $(OPT) tools.f90

tools.mod: tools.f90
	$(comp) -c $(OPT) tools.f90

ff_utils.o: ff_utils.f90 tools.f90
	$(comp) -c $(OPT) ff_utils.f90 tools.f90

ff_utils.mod: ff_utils.f90 tools.f90 tools.mod
	$(comp) -c $(OPT) ff_utils.f90 tools.f90

analysis.o: analysis.f90 tools.f90 ff_utils.f90 tools.mod ff_utils.mod
	$(comp) -c $(OPT) analysis.f90 tools.f90 ff_utils.f90

clean:
	rm -f tools.f90
	rm -f tools.mod
	rm -f ff_utils.mod
	rm -f $(objects)
