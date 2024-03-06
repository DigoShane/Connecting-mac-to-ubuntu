program vol
implicit none
real ::rad1,rad2,vol1,vol2
character :: response
do
print *, 'Enter the two radii'
read *, rad1,rad2
call volume(rad1,vol1)
call volume(rad2,vol2)
write(*,10) 'the difference in volume is, ', abs(vol1-vol2)
10 format(a,2f10.3)
print *, 'Any more?- hit Y for yes or else any other key'
read *, response
if (response /= 'Y' .and. response /= 'y') stop
end do
end program vol
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine volume(rad,vol)
implicit none
real :: rad,vol,pi
pi=4.0*atan(1.0)
vol=4./3. *pi*rad*rad*rad
end subroutine volume
!++++++++++++++++++++++++++_Notes_++++++++++++++++++++++++++++++++
! 1. There may be many subroutines in your program, ideally a 
! subroutine should do one task reflected by its name
! 2. All the variable sin the subroutine apart form the ones passed
! as argumens, are hidden form the main program. Thus we cna use the
! name in the subroutine as in the main program (the value stored is
! unaffected-unless the variable is passed as an argument)
! 3.Always use 'implicit none' in the subroutine
! 4. All variables in the subroutine must be declared 
! 5. The positioning of the arguments is important. It is vital that
! the arguments agree in both position and type. If an argument to 
! the subroutine is real in the main program, then it must be declared
! real in the subroutine as well.