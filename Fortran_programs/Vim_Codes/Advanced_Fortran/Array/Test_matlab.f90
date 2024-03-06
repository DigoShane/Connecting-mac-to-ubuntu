	program test_mtlab

	implicit none
	integer, dimension(2,3) :: A
	integer, dimension(2,1) :: b
	integer, dimension(3,1) :: c
	integer, allocatable, dimension(:) :: D
	integer :: sizeD
	integer :: i

	A=reshape((/1,4,7,2,5,8/),(/2,3/))
	b=reshape((/6,9/),(/2,1/))

	c(1:2,1)=A(2,2:3)
	c(3,1)=b(2,1)

	do i=1,3
	print *, c(i,1)
	end do

! now we want to set b to zero and print it
	b=0
	do i=1,2
	print *, b(i,1)
	end do

	print *, 'for i=1,n,2 and n is even, then it should stop at n-1 right?'
	do i=1,6,2
	print *, i
	end do

	print *, 'This is to test out the Shape command,b(:,1)=shape(A)'
	b(:,1)=shape(A)
	print *, 'The dim of A is', b(1,1),'x', b(2,1) 

	print *, 'This is a test to see when we should deallocate memory'
	print *, 'Enter the no. of elements of the 2D array D'
	read (*,*) sizeD
	allocate(D(sizeD))
	print *, 'Please enter the elements of D'
	do i=1,sizeD
	read *, D(i)
	end do
	print *,'Now we print each element' 
	do i=1,sizeD
	print *, D(i)
	end do
	print *, 'We now deallocate the memory'
	deallocate(D)
!Had we taken this route, when we would have received a segementation fault, thisproves that we must deallocate the memory at the very end, just as all the programs do. :)
!	print *, 'We now deallocate the memory'
!	deallocate(D)
!	print *,'Now we print each element' 
!	do i=1,sizeD
!	print *, D(i)
!	end do

	end program test_mtlab
