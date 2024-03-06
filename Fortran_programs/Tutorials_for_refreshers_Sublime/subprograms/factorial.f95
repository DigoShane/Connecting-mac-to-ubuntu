program factorial
integer k
k=8
print *, 'the value of',k,'! is', fact(k)
end program factorial
! This is a subfunction for the factorial function
function fact(n)
integer p,n,fact
p=1
do i=1,n
p=p*i
end do
fact=p
end