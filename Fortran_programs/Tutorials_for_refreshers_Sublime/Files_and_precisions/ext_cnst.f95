program extended_constants
! this prgram demonstrates the use of extended precision
implicit none
integer, parameter :: ikind=selected_real_kind(p=18)
real (kind=ikind) :: val,x,y
val=10/3
print *, val
print *, 'calculates an integer becuase both 10 and 3 are integer and division gives 3'
x=10.0
y=3.0
val=x/y
print *, '10/3 is', val
print *, 'x/y assigned to extended precision - right!'
val=10.0_ikind/3.0
print *, val
print *, 'extended precision const - right!'
val=10.0/3.0
print *, val
print *, 'real constant - wrong!'
val=0.12345678901234567890
print *, val
print *, 'real constant - wrong!'
val=0.12345678901234567890_ikind
print *, val
print *, 'extended precision constant - right!'
end program extended_constants