! This is a program to test function handles in fortran

	program fnHndl

	implicit none

	real*8, parameter :: piby2=2*atan(1.0)
	double precision :: f,a
	integer :: i,n

	f(a)=sin(a)

	print *, 'pi by 2 = ', piby2

	print *, 'enter the no. of divisions of pi'
	read (*,*) n

	 DO 10, I = 1,n
         print *, 'sin(', I,'/',n,'pi/2)=', f(I*piby2/n)
10       CONTINUE


	end program fnHndl
