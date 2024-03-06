program ImpliesDo
implicit none
real :: ra(3,3)
integer :: row,col
print *, 'please enter the 3x3 matrix'
do row=1,3
do col=1,3
read *, ra(row,col)
end do
end do
do row=1,3
write (*,2) (ra(row,col),col=1,3)
end do
2 format(10f5.1)
end program ImpliesDo