

psng.exe: PssngArray.o facts.o
	gfortran -o psng.exe PssngArray.o facts.o

PssngArray.o: PssngArray.f95 facts.o
	gfortran -c PssngArray.f95

facts.o: facts.f95
	gfortran -c facts.f95 
