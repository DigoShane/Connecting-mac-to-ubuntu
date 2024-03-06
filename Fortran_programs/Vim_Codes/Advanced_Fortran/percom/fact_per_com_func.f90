! This funciton is called by the program percom, this program caluclates the factorial, perumation and combination which is passed to the "main" program which prints the vlaues of nPr and nCr

! the key word recursive is to allow for the function to be a recursive type
	recursive integer*8 function fact(n) result(fact1) 
! the key word return tells the compiler that the function reutrns the variable fact1
! there is no need to specify the return type of the variable a

	implicit none

	integer*8, intent(in) ::n!INTENT(IN) in function's discussion. It indicates that an argument will receives some input from outside of the function and its value will not, actually cannot, be changed within the function. Since a subroutine cannot return a value through its name, it must return the computation results, if any, through its argument.


	if(n==1 .or. n==0) then

	fact1=1
	return

	else

	fact1=fact(n-1)
	fact1=fact1*n
	return

	end if

	end function fact



	integer*8 function per(n,r)

	implicit none
	integer*8, intent(in) :: n,r
	integer*8, external :: fact ! this is becuase the fact function is external to this per function


	per=fact(n)/fact(n-r)

	end function per



	integer*8 function comm(n,r)

	implicit none
	integer*8, intent(in) :: n,r
	integer*8, external :: fact ! this is becuase the fact function is external to this per function


	comm=fact(n)/(fact(n-r)*fact(r))

	end function comm



	integer*8 function com(per,r)

	implicit none
	integer*8, intent(in) :: per,r
	integer*8, external :: fact ! this is becuase the fact function is external to this per function


	com=per/fact(r)

	end function com
	


	integer*8 function comb(n,r)

	implicit none
	integer*8, intent(in) :: n,r
	integer*8, external :: fact, per


	comb=per(n,r)/fact(r) ! this is a multi funciton

	end function comb
