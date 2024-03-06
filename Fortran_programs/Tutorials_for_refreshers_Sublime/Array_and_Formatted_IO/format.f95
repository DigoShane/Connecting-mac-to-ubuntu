program format
implicit none
integer, parameter :: ikind=selected_real_kind(p=15)
real, dimension(4) :: x
integer, dimension(4) :: nums
integer :: i
real (kind=ikind),dimension(4) :: computed
do i=1,4
nums(i)=i*10
computed(i)=cos(0.1*i)
x(i)=computed(i)
! as x is stored as a real value, it doesn't have the same amount of precision as computed
end do
print *, 'nums-integer'
write(*,1) nums
1 format(2i10)
! 2 tells the comouter to print 2 variables in a row
print *, 'x - real'
write(*,2) x
2 format(f6.2)
print *, 'computed - double precision'
write(*,3) computed
3 format(f20.7)
end program format
! the usual syntax of the write statement is
! write(outputdevice,label) variable(s)
! label format(specification)
!  please check https://ww.fortrantutorials.com/documents/IntroductionToFTN95.pdf
! for details on exponential,character,etc specifications