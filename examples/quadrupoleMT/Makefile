
include ../Makefile.common

CXXFLAGS += $(DEBUGFLAG)

OBJ = Driver.o Model_QuadrupoleLD.o RLSolver_2DTableMT.o 

driver : $(OBJ)
	$(CXX) -o $@ $(OBJ) $(LDFLAG) 
	
	
clean:
	rm -f *.o *~
