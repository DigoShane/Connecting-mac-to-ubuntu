!Calcualte the area of a triangle using a function

	program area_qt

	implicit none

	real*8 :: area
	character :: fig
	real*8, external :: tri, sq !tells the compiler that a function will be defined externalls form whom it can collect values to be used in this program
	real*8 :: a,b,c,s

	print *, "This program calculates the area od a square/triangle"
	print *, "Type s for a square and tyoe t for a triangle"

	read *, fig

	if (fig =='t' .or. fig == 'T') then

	print *, "Enter the sides of the triangle"
	read *, a,b,c
	area = tri(a,b,c)

	print *, "The area of the triangle is: ", area

	else if (fig == 's' .or. fig =='S') then

	print *, "Enter the side of the square"
	read *, s
	area = sq(s)

	print *, "The area of the triangle is: ", area

	else

	print *, "Invalid Option"

	end if

	end program area_qt

! we cna define the program in side the same script, inside the program itself or as a separate script
! we define one function here, the tri function

	real*8 function tri(a,b,c)

	implicit none

	real*8, intent(in) :: a,b,c !these variables are external and being imported to this function. The in ensures that whatever values came into this function will remain unchanged

	real*8 :: s

	s = (a+b+c)*0.5
	tri = (s*(s-a)*(s-b)*(s-c))**0.5 !The function name and the variable passing the answer have to be the same
! functions can take multiple values inside but returns only one answer. For returning multiple variables we need subroutines

	print *, "The sides are:"
	print *, a,b,c

	end function tri
