program test
implicit none
real :: x,y,choice,answer
print *, 'enter the value for x'
read *, x
print *, 'enter the value for y'
read *, y
print *, 'enter one of the following options :-'
print *, '1. Multiply x.y'
print *, '2. Divide x/y'
print *, '3. Add x+y'
print *, '4. Subtract x-y'
print *, '5. Exponentiate x^y'
print *, '6. Sine sin(x)*y'
print *, '7. Cosine cos(x)*cos(y)'
print *, '8. tan and inverse tan tan(x)*atan(y)'
read *, choice;
if (choice == 1) then
answer=x*y
!print *, 'choice was' choice
!print *,'result = ',answer 
else if (choice == 2) then 
answer=x/y
!print *, 'choice was' choice
!print *,'result = ',answer 
else if (choice == 3) then
answer=x+y
!print *, 'choice was' choice
!print *,'result = ',answer 
else if (choice == 4) then
answer=x-y
!print *, 'choice was' choice
!print *,'result = ',answer 
else if (choice == 5) then
answer=x**y
!print *, 'choice was' choice
!print *,'result = ',answer 
else if (choice == 6) then
answer=sin(x)*y
!print *, 'choice was' choice
!print *,'result = ',answer 
else if (choice == 7) then
answer=cos(x)*cos(y)
!print *, 'choice was' choice
!print *,'result = ',answer 
else if (choice == 8) then
answer=tan(x)*atan(y)
!print *, 'choice was' choice
!print *,'result = ',answer 
else 
print *, 'Please enter a choice form the above list'
end if
!print *, 'Your choice was', choice , 'and thus the computed result is' answer
! Note tha above commented print option gives an error the one below works perfectly, thus the spacing is important
print *, 'You selected' ,int(choice),' The result is' ,answer
end program test