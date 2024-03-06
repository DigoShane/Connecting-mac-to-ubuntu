program avg
implicit none
real :: x(10),avg
integer :: i
print *, 'enter 10 nos.'
i=1
avg=0
do i=1,10,1
read *, x(i)
avg=avg+x(i)
end do
avg=avg/10
print *, avg
end program avg