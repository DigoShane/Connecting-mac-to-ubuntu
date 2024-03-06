program todara
implicit none
integer,dimension(2,3) :: a
integer :: row,col,count
count=0
do row=1,2
count=0
do col=1,3
count=count+1
a(row,col)=count
end do
end do
do row=1,2
do col=1,3
print *, a(row,col)
end do
end do
end program todara