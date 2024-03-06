 program circle
      !Notes on running fortran:-
      ! This code can be found at https://web.stanford.edu/class/me200c/tutorial_77/03_basics.html
      ! The lines beginning with '!' are commented, 'c' for some reason doesn't seem to work
      ! in the command terminal enter:
      ! gfortran -std=f95 <programname>.f95 -0 <executable file name>
      ! For this program, the above is
      ! gfortran -std=f95 programcircle.f95 -0 programme
      ! The above creates an executable file by the name <executable file name> (programme)
      ! To run the program, you have to execute the executable file which can be done by the follwoing command
      ! ./<executable file name> 
      ! which for this case would be
      ! ./programme
      ! The above runs the executable file and then you have your 
      real r, area
      ! This program reads a real number r and prints
      ! the area of a circle with radius r.
      write (*,*) 'Give radius r:'
      read  (*,*) r
      area = 3.14159*r*r
      write (*,*) 'Area = ', area
 	  ! I wanna get this done asap
      stop
      end