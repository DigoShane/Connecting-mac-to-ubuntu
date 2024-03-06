! Set of subroutines for demo



	subroutine B

	implicit none


	print *, "WE are inside Subroutine B!"
	print *, "calling Subroutine C!"

	call C

	print *, "We are back to subroutine B"


	end subroutine B
	

	subroutine C

	implicit none


	print *, "WE are inside C!"
	print *, "calling Subroutine D!"

	call D

	print *, "We are back to subroutine C"


	end subroutine C



	subroutine D

	implicit none


	print *, "WE are inside Subroutine D!"

	print *, "We are going back to the called subroutine!"


	end subroutine D
