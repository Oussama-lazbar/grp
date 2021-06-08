gfortran -c quad.f90
gfortran -c init.f90
gfortran -c quad.f90 init.f90 func.f90
gfortran -c func.f90 grp.f90
gfortran -o exe quad.o init.o func.o grp.o
./exe
#rm *.o n*.txt exe

