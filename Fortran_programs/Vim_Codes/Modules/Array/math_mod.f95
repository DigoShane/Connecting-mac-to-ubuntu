! this is the module
! Module is a collection of functions and subroutines in an orderly manner.
! the usefull ness of this is that the module need not be in the same directory
!==========================Some info about modules =======================================
!1. All functions inside the module can access each other. No "external" referencing required
!2. Module name and file name need not be the same.
!3(a). If the file names are different, then  then the .o file has the file name (math_mod.f95)
!	 and .mod file has the module name (maths1).
!3(b). We will use the .mod file name inside the call and .o filename for building and execution
!4. Module cnanot be written as part of the main program
!5. In the shell script, we must compile the module file first and then the main program. 
!   This way the compiler knows where to look for the module when it comes across the 
!   use <module_name> command in the main program.
!6. When functions are written inside module, they can return arrays for an example, check fact_fun


	module maths1 ! BEcuase of point 2. check by :echo @%

	contains


! Subroutine for the program crammer.f95

	subroutine cram(x,a,b,n)

	implicit none
	integer, intent(in) :: n
	real*8, dimension(n,n,n+1), intent(inout) :: a !this allows for a to be both inout output
	real*8, intent(in), dimension(n) :: b
!	real*8, external :: detf   ! Due to point 1. If this was not a module this would be required.
	integer :: i
	real*8, dimension(n,1), intent(out) :: x

	do i=1,n

	a(:,:,i+1)=a(:,:,1)
	a(:,i,i+1)=b
	x(i,1) = detf(a(:,:,i+1),n)/detf(a(:,:,1),n)

	end do

	end subroutine cram

! this function recursively calculates the determinant
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

! This subroutine gives us back the cofactor/ or minor, which ever you prefer , here they call it slice
	subroutine slicef(sl, mat, n, row, column)

	implicit none
	integer, intent(in) :: n,row,column
	real*8, dimension(n,n), intent(in) :: mat
	real*8, dimension(n-1,n-1), intent(out) :: sl! this is the output value
	logical, dimension(n,n) :: mask


	mask = .true.
	mask(row,:)= .false.
	mask(:,column)= .false.

	sl = reshape(pack(mat,mask),(/n-1,n-1/))

	end subroutine slicef

	end module maths1
