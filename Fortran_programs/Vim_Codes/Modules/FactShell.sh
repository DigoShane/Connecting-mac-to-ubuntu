clear
ls
rm *o *exe *mod
ls
gfortran -c  fact1.f95
ls
gfortran -c fact_fun.f95
ls
gfortran fact1.f95 fact_fun.f95 -o factorial.exe
ls
./factorial.exe
