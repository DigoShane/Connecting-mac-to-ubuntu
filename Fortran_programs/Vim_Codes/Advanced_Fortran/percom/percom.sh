clear
ls
rm *o *.exe *.mod 
ls
gfortran -c *.f90 
ls
gfortran *.o -o percom.exe
ls
./percom.exe
