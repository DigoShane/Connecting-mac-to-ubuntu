!=====================================================================
! To find the length and angle of a line joining a point and the origin


      program quadrant
      implicit none
      real x,y,l,a
      real, parameter :: pi=4*atan(1.0)
      !	To keep some values fized, we define them as parameters

	print *, "Enter the vlae fo coordinate, x and y :"
	read *, x,y

!      pi=6;! this will give na error 'pi has no IMPLICIT type'
      ! this is because this was defined to be a constant.
      
      l=sqrt(x**2+y**2)
      a=atan(y/x)


	print *, "l = ", l
        print *, "a = ", a


      end program quadrant
