program divisors
! 	This program finds the divisors of an integer input by the user.
! 	The divisors are printed to a file.
integer n, k, d(10)
open (unit = 1, file = "divisors")
print *, "Enter a positive integer :"
read *, n
write (1,*) "Here are the divisors of ", n, " :"
k = 0
do i = 1, n
if (mod(n,i) .eq. 0) then
! .eq. is the logical expression for equality, .eq.\equiv ==
k = k + 1
d(k) = i
end if
if (k .eq. 10) then
write (1,5) (d(j), j = 1, 10)! check Fortran Lesson 3. pdf, this input is listed ther, basically it prints the array in a single line in the output file '1'.
k = 0
end if
end do
write (1,5) (d(j), j = 1, k)
! To write in fortran we use the format write (*,20).
! * tells fortran to write to the screen
! 20 refers to the label of the format statement fir the write command
! An example of such a format will be 
! 20 format (3f10.4)
! 3 refers to the fact that 3 entries will be printed
! f denotes the fact that they will be floating point real numbers
! .4 mandates that there will be 4 digits after the decimal place
! thus examples include 12345.6789, -1234.5678, 10002.3400
! the letter 'f' is referred to as the format code letter
! some common examples cna be dounf in the Fortran Lesson 3.pdf file in this folder.
5 format (10i7)
close (1)
! Thus weite (1,5) does the same thing as write (*,5) but the output is in the file numbered '1'.
print *, "The divisors are listed in the file 'divisors'. Bye."
end