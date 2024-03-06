! This is a Program to test passing Arrays through subroutines

	program pssngArray

	use facts

	implicit none
	integer,parameter :: r=2,c=3
	integer :: i,j,a,b
	double precision, dimension(1:r,1:c) :: temp

	a=r
	b=c
	call initialize(temp,r,c)
!	temp=0.0

	do i=1,2
	do j=1,3
	print *, temp(i,j)
	enddo
	enddo

!	print *, r
!	print *, c


	end program
