! program to find permutaitons and combinations via functions

	program percom

	implicit none

	integer*8, external :: per, comm, com, comb ! *8 is short kand for writing (kind=8)
! the key word external is used to define eternal functions to the program, thus in this case
! per comm and com are programs which return integers
	integer*8 :: n,r,per1,com1,com2,com3

	print *, "Enter the value of n:"
	read *, n

	print *, "Enter Values of r"
	read *, r

		per1=per(n,r)
		com1=com(per1,r)
		com2=comm(n,r)
		com3=comb(n,r)

	print *, "The permutation of ", r," in", n, "is :", per1
	print *, "The combination of ", r," in", n, "is :", com1
	print *, "The combination of ", r," in", n, "is :", com2
	print *, "The combination of ", r," in", n, "is :", com3

	end program percom


