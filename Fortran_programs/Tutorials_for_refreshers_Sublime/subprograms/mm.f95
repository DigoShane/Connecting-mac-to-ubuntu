program mm
implicit none
real :: a(2,2),b(2,2),c(2,2),d
integer :: i,j
print *, 'please enter the first matrix'
do i=1,2
do j=1,2
read *, a(i,j)
end do
end do
print *, 'please enter the second matrix'
do i=1,2
do j=1,2
read *, b(i,j)
end do
end do
c=MATMUL(a,b)
write(*,1) 'The two matrices multiplies is'
call printmat(c)
write(*,1) 'The first matrix transposed is'
c=TRANSPOSE(a)
call printmat(c)
write(*,1) 'The min val of the second matrix is'
d=MINVAL(b)
print *, d
write(*,1) 'The sum of the second matrix is'
d=SUM(b)
print *, d
1 format(a,2i2)
end program mm
!+++++++++++++++++++Subroutine+++++++++++++++++
subroutine printmat(a)
implicit none
real :: a(2,2)
integer :: i,j
do i=1,2
do j=1,2
print *, a(i,j)
end do
end do
end subroutine printmat