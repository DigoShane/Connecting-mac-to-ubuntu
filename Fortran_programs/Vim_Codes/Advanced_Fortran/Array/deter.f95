
! the determinant of a matrix

	program deter

	implicit none
	integer, parameter :: n=3 !Size of the matrix
	real*8, dimension(n,n) :: mat
	real*8, external :: detf! This is a function to be called, more specifically , this is a recursive function.
	real*8 :: det
	integer :: i,j


!	mat= reshape((/1.0,2.0,0.0,2.0,-4.0,0.0,3.0,1.0,-1.0/),(/3,3/)) ! det=1
	mat= reshape((/3,2,0,2,-4,0,3,1,-1/),(/3,3/)) ! det=4
!	mat= reshape((/2,2,3,4,5,6,7,8,9/),(/3,3/)) ! det=0
!	mat= reshape((/-2,-1,2,2,1,0,-3,3,-1/),(/3,3/)) ! det=18
!	mat= reshape((/6,4,2,1,-2,8,1,5,7/),(/3,3/)) ! det=306
!	mat= reshape((/5,0,2,-2,3,0,1,-1,7/),(/3,3/)) ! det=103

	     det=detf(mat,n)

	print *, "The matrix is"

	do i=1,n

		do j=1,n

		write(*,'(f8.3,t3)', advance='no'),mat(i,j) ! This is used to supress linebreaks???

		end do

	write(*,*)

	end do


	print *, "The determinant is = ", det

	end program deter


	recursive real*8 function detf(mat,n) result(det)

	implicit none
	integer, intent(in) :: n !to keeo the values fixed, we use intent(in)
	real*8, intent(in), dimension(n,n) :: mat 
	real*8, dimension(n-1,n-1) :: sl

	integer :: i

	det=0

	if (n == 2) then

	det=mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)

	return

	else if (n==1) then

	det = mat(1,1)
	return

	else

		do i = 1,n

			call slicef(sl,mat,n,1,i)
			det=det+((-1.0)**(i+1)*mat(1,i)*detf(sl,n-1))

		end do

	return

	end if
	end function detf


	subroutine slicef(sl, mat, n, row, column)

	implicit none
	integer, intent(in) :: n,row,column
	real*8, dimension(n,n), intent(in) :: mat
	real*8, dimension(n-1,n-1), intent(out) :: sl! this is the output value
	logical, dimension(n,n) :: mask


	mask = .true.
	mask(row,:)=.false.
	mask(:,column)=.false.

	sl = reshape(pack(mat,mask),(/n-1,n-1/))

	end subroutine slicef
