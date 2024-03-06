! Subroutine for the program crammer.f95

	subroutine cram(x,a,b,n)

	implicit none
	integer, intent(in) :: n
	real*8, dimension(n,n,n+1), intent(inout) :: a !this allows for a to be both inout output
	real*8, intent(in), dimension(n) :: b
	real*8, external :: detf
	integer :: i
	real*8, dimension(n,1), intent(out) :: x

	do i=1,n

	a(:,:,i+1)=a(:,:,1)
	a(:,i,i+1)=b
	x(i,1) = detf(a(:,:,i+1),n)/detf(a(:,:,1),n)

	end do

	end subroutine cram


	recursive real*8 function detf(mat,n) result(det)

	implicit none
	integer, intent(in) :: n !to keeo the values fixed, we use intent(in)
	real*8, intent(in), dimension(n,n) :: mat 
	real*8, dimension(n-1,n-1) :: sl
!	real*8 :: det1
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
