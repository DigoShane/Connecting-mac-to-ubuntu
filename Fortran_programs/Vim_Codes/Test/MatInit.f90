! This is a program to check if we can initialize the matrix to zero

	program matinit

	implicit none
	integer, parameter :: r=2,c=3
	integer :: i,j
	double precision, dimension(1:r,1:c) :: temp

	temp = 0.0
	print *, temp
	print *, r,c
	print *, sizeof(temp)

	do i=1,r
	do j=1,c
	print *, temp(i,j)
	enddo
	enddo

	end program matinit
