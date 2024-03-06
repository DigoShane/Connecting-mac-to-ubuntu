	program String_array

	implicit none

! We write down the different ways of defining a string
	character :: a 
!thefollowing 4 definitions are equivalnet
	character (len=5) :: b !string that can store 5 character
	character(5) :: be
	character*5 :: bf
	

	character :: c*6, d*8, e*3 ! *6 is the size of the character array


! Array of Strings
	character (len=6), dimension(2) :: Str1 !a 1by2 array with each unit being a string made up of 6 characrters
	character , dimension(2) :: Str2


	a= 'X'
	b= 'white'
	c= '456321'
	d= "asdghft6"
	e= "_ _"

	print *, "a = ",a
	print *, "b = ",b
	print *, "c = ",c
	print *, "d = ",d
	print *, "e = ",e

	Str1(1) = "Carrom"
	Str1(2) = "Codng"

	print *, "First String is ", Str1(1), "Second String is", Str1(2)


	end program String_array
