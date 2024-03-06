


      program formats
      implicit none
      integer :: a
      real :: ar
      !character :: ac
      
      read *, a,ar
      5 format(i6,f8.3)
      
      10 format("With formating", i5, 4X,f5.2)

      print 10, a,ar
      print 5, a,ar
      print *, "Without formatting", a,ar
      end program formats












