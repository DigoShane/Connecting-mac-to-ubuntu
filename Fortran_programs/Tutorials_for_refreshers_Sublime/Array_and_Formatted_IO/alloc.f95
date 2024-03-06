program alloc
implicit none
integer, allocatable,dimension(:):: A 
!note syntax - dimension(:) 
integer  :: elements,i
print *,'enter the number of  elements in the A'
read *,elements    
allocate(A(elements))  
!allocates the correct amount  of memory 
print *,' your vector is of  size ',elements,'. Now enter each element'
do i=1,elements
read *,A(i)
end do
print *,'This is your vector'
do i=1,elements
print *,A(i)
end do            
deallocate(A) 
!tidies up the memory
end program alloc