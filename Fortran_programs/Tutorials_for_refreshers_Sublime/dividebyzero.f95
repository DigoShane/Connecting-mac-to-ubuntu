program dividebyzero
implicit none
real :: ans
integer :: i,x 
print *, 'the value of x'
read *, x
do i=-x,x,1
if (i == 0) then
print *, 'cannot divide by zero'
else
ans=1.0/i
print *, ans
end if
end do
end program dividebyzero