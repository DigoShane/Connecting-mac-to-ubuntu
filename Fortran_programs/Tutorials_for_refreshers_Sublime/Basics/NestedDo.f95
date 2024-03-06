program NestedDo
implicit none
real :: x,y,ans
integer :: i,j
x=1.0;
y=1.0;
do i=0,2
do j=0,2
ans=(x+i*0.5)**(y+j*0.5)
print *, ans
end do
end do
end program NestedDo