program avg
implicit none
integer, parameter :: imax=10
real,dimension(imax) :: x,y
real :: aveg
integer :: i
print *, 'enter' ,imax, 'numbers'
i=1
aveg=0
do i=1,imax,1
read *, x(i)
aveg=aveg+x(i)
end do
aveg=aveg/imax
print *, 'enter' ,imax, 'numbers'
i=1
aveg=0
do i=1,imax,1
read *, y(i)
end do
print *, aveg
print *, 3*y+2*x
end program avg