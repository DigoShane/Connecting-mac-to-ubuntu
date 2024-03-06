program arithmetic
implicit none
real :: x,y,z,a,b,answer1,answer2
print  *, 'Enter 5 numbers'                        
read  *, x                                           
read  *, y
read  *, z
read  *, a
read  *, b
! ./. is division
! .+. is addition
! .*. is multiplication
! .**. is exponentation
! brackets () help ordering the operations
answer1=(x/a+y*b)**z
answer2=x**y**z
print  *, 'The total is ', answer1
print  *, 'The answer for x^(y^z) is ', answer2
! For a list of all the functions in fortran, http://www.ldeo.columbia.edu/~mspieg/mmm/Fortran.pdf
end program arithmetic