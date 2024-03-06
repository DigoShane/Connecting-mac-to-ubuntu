! program for basic subroutine

! a funciton cna return only one value immaterial of the no. of arguments (need not be 1)
! a subroutine can take in multiple vlaues and change multiple values at the same time. A subroutine cannot return a value, it just modifies the values when called inside the program
! Important for a large no. of programs.

! A function needs parameters passed into them, while subroutines need not have parameters passed into them

	program sub_ex

	implicit none

	print *, "WE are inside main program!"
	print *, "calling Subroutine A!"

	call A

	print *, "We are back to the main program and exiting!"


	end program sub_ex

!! If we wanted to include the subroutine as part of the main function, then we would have to use the keyword 'contains'
!
!	subroutine B
!
!	implicit none
!
!
!	print *, "WE are inside Subroutine B!"
!	print *, "calling Subroutine C!"
!
!	call C
!
!	print *, "We are back to subroutine B"
!
!
!	end subroutine B
!
!! Including a subroutine in the main function makes it local to that function and thus it can't
!! call other subroutines outside its scope. Thus Unforutnately as the above subroutine calls a
!! subroutine outsie its scope, the above will not get executed. However, if the above was not
!! calling outside its scope it would run perfectly.

	

	subroutine A

	implicit none


	print *, "WE are inside Subroutine A!"
	print *, "calling Subroutine B!"

	call B

	print *, "We are back to subroutine A"


	end subroutine A
