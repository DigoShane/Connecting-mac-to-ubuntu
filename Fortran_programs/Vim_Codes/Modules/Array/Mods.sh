! this shell script has everything
clear
ls
rm *mod *o *exe
ls
gfortran -c math_mod.f95
ls
gfortran -c Crammer.f95
ls
gfortran Crammer.o math_mod.o -o Crammer.exe
ls
./Crammer.exe
