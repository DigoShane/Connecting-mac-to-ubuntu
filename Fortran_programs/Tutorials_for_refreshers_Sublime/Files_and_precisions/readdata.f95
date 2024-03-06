program readdata
implicit none
!reads data from a file called mydata.txt
real :: x,y,z
open(2,file='mydata.txt')
read(2,*) x,y,z
print *,x,y,z
end program  readdata 
!A good read http://www.math.hawaii.edu/~hile/fortran/fort7.htm