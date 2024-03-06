program hello_world
include '/usr/include/mpif.h'
integer ierr

call MPI_IINIT ( ierr )
print *, ''Hello world''
call MPI_FINALIZE ( ierr )

stop
end
