program factorial
implicit none
integer :: i,n,ans
ans=1;
print *, 'enter the integer whose factorial you wanted evaluated'
read *, n
! In general the do statement if of the form
! do i=<start value>,<end value>,<increment>
!      statements
! end do
do i=1,n
ans=ans*i
end do
print *, 'the value of',n,'! is', ans
end program factorial