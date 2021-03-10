objects = water.o tools.o ff_utils.o
objectsbis = MD.o tools.o ff_utils.o
comp = gfortran
OPT = 

water.x: $(objects)
	$(comp) -o water.x $(OPT) $(objects) -llapack

MD.x: $(objectsbis)
	$(comp) -o MD.x $(OPT) $(objectsbis)

tools.o: ../tools.f90
	$(comp) -c $(OPT) ../tools.f90

ff_utils.o: ../ff_utils.f90 ../tools.f90
	$(comp) -c $(OPT) ../ff_utils.f90 ../tools.f90

water.o: water.f90 ../tools.f90 ../ff_utils.f90 tools.o ff_utils.o
	$(comp) -c $(OPT) water.f90 ../tools.f90 ../ff_utils.f90 -llapack

MD.o: MD.f90 ../tools.f90 ../ff_utils.f90 tools.o ff_utils.o
	$(comp) -c $(OPT) MD.f90 ../tools.f90 ../ff_utils.f90
clean:
	rm -f tools.mod
	rm -f ff_utils.mod
	rm -f $(objects)
	rm -f $(objectsbis)
