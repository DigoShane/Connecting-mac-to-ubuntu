!This module is to be used by pssngArray

	module facts

	contains

	subroutine initialize(temp,r,c)

	implicit none
	integer, intent(in) :: r,c
	double precision, dimension(1:r,1:c),intent(inout) :: temp

	temp=0.0

	end subroutine initialize

	end module facts
