objects = nitromethane.o tools.o ff_utils.o
objectsbis = MD.o tools.o ff_utils.o md_utils.o
comp = gfortran
OPT = -O3 -fbounds-check

nitromethane.x: $(objects)
	$(comp) -o nitromethane.x $(OPT) $(objects) -llapack

MD.x: $(objectsbis)
	$(comp) -o MD.x $(OPT) $(objectsbis) -llapack

tools.o: ../tools.f90
	$(comp) -c $(OPT) ../tools.f90

ff_utils.o: ../ff_utils.f90 ../tools.f90
	$(comp) -c $(OPT) ../ff_utils.f90 ../tools.f90

md_utils.o: ../md_utils.f90 ../ff_utils.f90 ../tools.f90
	$(comp) -c $(OPT) ../md_utils.f90 ../ff_utils.f90 ../tools.f90

nitromethane.o: nitromethane.f90 ../tools.f90 ../ff_utils.f90 tools.o ff_utils.o
	$(comp) -c $(OPT) nitromethane.f90 ../tools.f90 ../ff_utils.f90 -llapack

MD.o: MD.f90 tools.o ff_utils.o md_utils.o
	$(comp) -c $(OPT) MD.f90 ../tools.f90 ../ff_utils.f90 ../md_utils.f90 -llapack
clean:
	rm -f tools.mod
	rm -f ff_utils.mod
	rm -f md_utils.mod
	rm -f $(objects)
	rm -f $(objectsbis)
