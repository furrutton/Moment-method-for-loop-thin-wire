PROGRAM main
  INCLUDE 'mpif.h'
  ! use mpi
  integer message(3)
  integer myrank,ierr,rc
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
  if (myrank==0) then
    message(1) = 5
    message(2) = 6
    message(3) = 7
    call MPI_SEND(message,3,MPI_INTEGER,1,99,MPI_COMM_WORLD,ierr)
  else if (myrank==1) then
    call MPI_RECV(message,3,MPI_INTEGER,0,99,MPI_COMM_WORLD,status,ierr)
    write(*,*) "received:",message
  end if
  call MPI_FINALIZE(rc)
end program main
