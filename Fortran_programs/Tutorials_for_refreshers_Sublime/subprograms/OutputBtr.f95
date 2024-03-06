program OutputBtr
implicit none
real,dimension(3) :: a,b,c
a=1.5 ! sets all the rows of the ector as 1.5
b=2.5
c=3.5
write(*,1) 'a', a
call prompt()
write(*,1) 'b', b
call prompt()
write(*,1) 'c', c
call prompt()
write (*,1) 'a*b*c', a*b*c
1 format(a,3f8.3)
! Check Fortran Lesson 3, a specifies that its the format for a string
end program OutputBtr
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine prompt()
! prompts for key press
implicit none
character :: answer*1
print *, 'type y to continue or any other key to finish'
read *, answer
if (answer /='y') stop
end subroutine prompt