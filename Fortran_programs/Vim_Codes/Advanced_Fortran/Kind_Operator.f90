	program kind
	implicit none
	integer(kind =8) :: a,r,i,n ! kind is used ot indicate that these are integers using 4 bytes
! with kind=4, and a=4, r=8, the series starts giving 0 vlaues, set kind=8 and that error goes away, more data cna be stored


	print *, "Enter the vlaue of the initial term, "
	read *, a,r

	print *, "Enter the no. of terms,"
	read *, n

	do i=0,n,1

	print *, "The value ", i, " of the series is:", a*(r**i)


	end do




	end program kind
