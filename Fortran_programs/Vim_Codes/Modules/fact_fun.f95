! this is a program to caculate the factorial of a no.

	program fact_fun

	use fact1 ! Calling the module fact1

	implicit none
	integer, parameter :: n=10
	integer*8 :: i
	integer*8, dimension(n) :: num, numfact

	num(1:n)=[(i, i=1,n)] ! allocating the 1 to n values fo the array num with the the corresponding no.
	numfact = fac_num(num,n)

	do i=1,n

	print *, "The",i,"! is:",numfact(i)

	end do


	end program fact_fun
