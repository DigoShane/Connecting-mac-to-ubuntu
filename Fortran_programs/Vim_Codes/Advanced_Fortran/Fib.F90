! This is a Fibionacci seires using positive integers

	program fib

	implicit none

!	integer, dimension(10) :: fib_vals !For this one fib values are bounded by defining the array to have size 40bytes (integer has 4 bytes)
	integer(kind=16), allocatable, dimension(:) :: fib_val,c_vals! does dynamics allocation based on the process at hand
	integer, allocatable, dimension(:) :: nvals
	integer :: i,n


	read *, n
	print 5,n
	5 format("Enter the max index no. of the fibonacci series to be &
	& displayed", i3)

	allocate(fib_val(n), c_vals(n), nvals(n)) ! Allocates n integers to the array at run time

	open(unit=1, file="fib_Val.txt") !opens a file called fib_Val (if not already present, then it is created) and it assigns a no. 1 to this file

	15 format(i3, 5x, i25)

	fib_val(1)=1
	fib_val(2)=1

	do i=3,n

	fib_val(i)=fib_val(i-1)+fib_val(i-2)

	end do


	do i=1,n

	print 10, i, fib_val(i)
	10 format("The fibonacci value ", i3, " is:", i25)
	write(unit=1, fmt=15), i, fib_val(i) ! writing values to unit=1(previously identified), fmt refers to the format in which it will be extracted
! equivalenty we can write write(1,15), i, fib_vals(i)
	end do

	close(unit=1) ! Think of this as Good programming practise
	open(unit=2, file= "fib_Val.txt")

	do i=1,n

	read (unit=2, fmt =15), nvals(i), c_vals(i) ! the read command aslo has an internal formatting. 15 was the formatting we used to write the file
	print *, "i = ", i, "i-n = ", i-nvals(i), "fib-fib_check = ", &
			&  fib_val(i)-c_vals(i)

	end do

	deallocate(fib_val, c_vals, nvals) ! Deallocates the array
	close(unit=2)

	end program fib
