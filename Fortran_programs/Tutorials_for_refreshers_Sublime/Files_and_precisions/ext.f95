program ext
implicit none
! we set the precision we want by the argument 'p' below. In this case
! we are looking for 15 digit precision
integer, parameter :: ikind=selected_real_kind(p=17)
! 'ikind' is a new kind of data type, a parameter
! Fortran returns a value to the parameter ikind that will be adequate
! to return 15 digit precision
real (kind=ikind) :: sum,x
integer :: i
sum=0.0
do i=1,100
x=i
sum=sum+1.0/(x**6)
end do
print *, sum
end program ext