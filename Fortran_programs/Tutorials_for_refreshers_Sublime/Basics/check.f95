program check
!Integer and real arithmetic
implicit none
real :: x,y
integer i
x=2.0
i=2
y=x*((2**i)/3) ! here fortran calculates (2**i) as n integer thus 4 then it does 4/3 whose integer part is 1 This *2 is just 2.
print *,y
y=x*((2.0**i)/3) ! here the 2.0 signals that everything is being calculated using real nos. thu you get a much better answer.
print *,y
end program check