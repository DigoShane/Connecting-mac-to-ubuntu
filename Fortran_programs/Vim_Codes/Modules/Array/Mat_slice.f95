! Matrix slicing

	program mat_slice

	implicit none
	real*8, dimension(3,3) :: mat
	real*8, dimension(2,2) :: slice
	logical, dimension(3,3) :: mask ! logocal matrix
!	real*8, dimension(4) ::temp
	integer :: i=3,j=3

	mat = reshape((/1,4,7,2,5,8,3,6,9/),(/3,3/)) ! fragemnts the 1d array into 3 rows each with 3 cols and each will get stacked to mat

	mask =.true. ! this sets each value to be true
	mask(i,:) = .false.
	mask(:,j) = .false.
! after this step, teh 3rd row (i=3) and the 3rd colum (j=3) gets set to false


	slice = reshape(pack(mat, mask),(/2,2/)) 
! The pack function will take a matrix (mat) and compares each and every index position to another index position of a mask matrix (whihc should be a logical matrix)
! the only true arguments of the matrix is the 2by2 cofactor of the 3,3 element. Thus pack returns a 1by4 array. This is the input ot the reshape command, which reshapes it into
! a 2by2 matrix.

	do i=1,3

	do j=1,3

	print *, '("Mat(,", i1,"," ,i1,") is = ", f10.3)', i,j,mat(i,j)

	end do

	end do

	end program mat_slice

! incomplete program
