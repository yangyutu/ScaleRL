# depends on FORTRAN compilers, command to invloke the fortran compiler
FC = gfortran

#please check your compiler for these options
FFLAGS = -O3

2d  =   main.o ranf.o testoverlap3D.o commondataset.o calorder.o getrotate.o rotate3D.o\
        $(libs)



bdpd: $(2d)	  
	$(FC) $(FFLAGS) -o $@ $^

main.o: main.f commondataset.o
	$(FC) -c main.f
ranf.o: ranf.f
	$(FC) -c ranf.f
testoverlap3D.o:testoverlap3D.f commondataset.o
	$(FC) -c testoverlap3D.f
calorder.o:calorder.f commondataset.o
	$(FC) -c calorder.f
rotate3D.o:rotate3D.f commondataset.o
	$(FC) -c rotate3D.f
getrotate.o:getrotate.f commondataset.o
	$(FC) -c getrotate.f
commondataset.o: commondataset.f
	$(FC) -c commondataset.f

clean:   
	rm -f core* bdpd *.o *.il *.mod




