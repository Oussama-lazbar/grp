F90 = gfortran
OFLAGS = -O3 
DFLAGS = -g -fcheck=all -Wall -fbounds-check

CFLAGS = $(OFLAGS)

grp :  quad.o init.o func.o grp.o
	$(F90) $(CFLAGS) -o $@ $^
	
test : quad.o init.o func.o test.o
	$(F90) $(CFLAGS) -o $@ $^

init.o : init.f90
	$(F90) -c $(CFLAGS) $^

quad.o : quad.f90 
	$(F90) -c $(CFLAGS) $^

func.o : func.f90 init.f90 quad.f90 
	$(F90) -c $(CFLAGS) $^

grp.o : grp.f90 func.f90
	$(F90) -c $(CFLAGS) $^

test.o : test.f90 func.f90
	$(F90) -c $(CFLAGS) $^


clean :
	rm -f *.o sol.txt grp test *.png n*.txt *.dat

plot :
	mkdir -p plots/ordre2/
	gnuplot -persist "plot.gn"
