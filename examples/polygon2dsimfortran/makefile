# depends on FORTRAN compilers, command to invloke the fortran compiler
FC = gfortran

#please check your compiler for these options
FFLAGS = -g3

2d  =   main.o ranf.o testoverlap.o commondataset.o calorder.o\
        $(libs)



bdpd: $(2d)	  
	$(FC) $(FFLAGS) -o $@ $^

main.o: main.f commondataset.o
	$(FC) -c $(FFLAGS) main.f
ranf.o: ranf.f
	$(FC) -c $(FFLAGS) ranf.f
testoverlap.o:testoverlap.f commondataset.o
	$(FC) -c $(FFLAGS) testoverlap.f
calorder.o:calorder.f calorder.o
	$(FC) -c $(FFLAGS) calorder.f

commondataset.o: commondataset.f
	$(FC) -c $(FFLAGS) commondataset.f

clean:   
	rm -f core* bdpd *.o *.il *.mod




