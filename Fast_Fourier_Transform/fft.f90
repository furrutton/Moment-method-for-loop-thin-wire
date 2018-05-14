module constants

  integer*4, parameter :: n = 2**3

end module constants


module tools

  complex*16, parameter :: j = (0.0d0, 1.0d0)
  real*8, parameter :: pi = 3.141592653589793d0

contains
  recursive subroutine fft(x)
    implicit none
    complex*16, dimension(:), intent(inout) :: x
    complex*16, allocatable :: even(:), odd(:)
    complex*16 :: t
    integer*4 :: n, i

    n = size(x)
    if ( n .le. 1 ) return

    allocate(odd((n)/2))
    allocate(even(n/2))

    odd  = x(1:n:2)
    even = x(2:n:2)

    call fft(odd)
    call fft(even)

    do i = 1, n/2, 1
      t = exp(-2.0d0*j*pi*(i-1)/n)*even(i)
      x(i)     = odd(i)+t
      x(i+n/2) = odd(i)-t
    end do

    deallocate(odd)
    deallocate(even)
    return
  end subroutine fft

  subroutine ifft(x)
    implicit none
    complex*16, dimension(:), intent(inout) :: x
    integer*4 :: n

    n = size(x)

    x = conjg(x)
    call fft(x)
    x = conjg(x)/n
    return
  end subroutine ifft

end module tools

program main
  use constants
  use tools
  implicit none
  complex*16 :: a(n)

! ------------------------------------------------------
! 判断 n 的取值是否为 2 的幂次
! ------------------------------------------------------
  if ( mod(n, 2) .ne. 0.0d0 ) then
    write(unit=*, fmt=*) "N SET ERROR!!!"
    goto 100
  end if
! ------------------------------------------------------

  data a / (0.0d0, 0.0d0), (1.0d0, 0.0d0), (2.0d0, 0.0d0), (3.0d0, 0.0d0), &
        &  (4.0d0, 0.0d0), (5.0d0, 0.0d0), (6.0d0, 0.0d0), (7.0d0, 0.0d0) /


  

  call fft(a)
  call ifft(a)
  write(unit=*, fmt=*) a

  write(unit=*, fmt=*) a*a

  write(unit=*, fmt=*) "SUCCED!!"
100  stop
end program main
