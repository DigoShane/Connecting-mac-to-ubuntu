!This can be found at https://www.fortrantutorial.com/basics/index.php
!In mathematics, “x = 2” means that the variable x is equal to 2.
!In FORTRAN it means “store the value 2 in the memory location that we have given the name x”.
!Thus z=x+y means we add the values at memory location x and memory location y and store
! the result in memory location z. Thus x+y=z, doesn't make sense 
program sum                                          !a: name of program
!an example of program  structure                    !b: a  comment   
real :: answer,x,y                                   !c: declarations  a good read is availabel at "https://en.wikibooks.org/wiki/Fortran/Fortran_variables"
print  *, 'Enter two numbers'                        !d: output
read  *, x                                           !e: input
read  *, y                                           !e: input
answer=x+y                                           !f: arithmetic
print  *, 'The total is ', answer                    !g: output
end  program sum                                     !h: end of program