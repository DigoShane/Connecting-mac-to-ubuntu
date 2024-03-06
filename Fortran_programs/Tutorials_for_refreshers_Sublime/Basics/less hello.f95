 $ less hello.f95
 write(*,*) 'Hello, world.'
 end
 $ f95 -0 hello hello.f95
 $./hello
 Hello, world.
 $