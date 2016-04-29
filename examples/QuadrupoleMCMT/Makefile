
include ../Makefile.common

CXXFLAGS += $(DEBUGFLAG) -I../quadrupoleMT
OBJ= RLSolver_2DTableMT.o Model_QuadrupoleMC.o mainC.o
#please check your compiler for these options

LDFLAG += -lgfortran

test:$(OBJ)
	$(CXX) -o $@ $^ $(LDFLAG)
	
test_static:$(OBJ)
	$(CXX) -static -o $@ $^ $(LDFLAG) -lquadmath -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
# here the linking flag is from https://gcc.gnu.org/ml/gcc-help/2010-05/msg00029.html

%.o:%.cpp
	$(CXX) -c $(CXXFLAGS) $^

RLSolver_2DTableMT.o: RLSolver_2DTableMT.cpp
	$(CXX) -c $(CXXFLAGS) $^

clean:   
	rm -f core* test test_static *.o *.il
