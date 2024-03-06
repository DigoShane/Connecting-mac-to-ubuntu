program convert
implicit none
!Old FORTRAN compilers used an implicit convention that integers have names starting with
!the letters in the range i – n, all the others being real.
!FORTRAN still allows you to do this if we don't include the line, implicit none
integer :: pounds,pence,total
!character :: name*10
!The *10 tells the program that you want to input 10 characters. 
! An alternatate way is the following
CHARACTER(LEN=20) :: name
print *,'what is your name?'
read *,name
print *, 'Hi ',name,'! Enter number of pounds and  pence'
read *, pounds,pence
total =100 * pounds + pence
print *,'the total money in pence is ',total
end program convert