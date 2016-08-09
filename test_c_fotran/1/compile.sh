gfortran -c *.f90
g++ -c *.cpp
gfortran *.o -lstdc++ -o a.out
