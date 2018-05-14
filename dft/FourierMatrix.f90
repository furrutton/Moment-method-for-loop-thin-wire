module tools
    real*8, parameter :: pi = 3.141592653589793d0
    contains
    function genFM(n)    ! 用来产生 N 阶傅里叶矩阵
        integer*4 :: n, k, j
        complex*16 :: genFM(n, n)
        complex*16, parameter :: i = (0.0d0, 1.0d0)
        do j = 1, n
            do k = 1, n
                genFM(k, j) = exp(-2*pi*i*(j-1)*(k-1)/n)
            enddo
        enddo
    end function genFM

    function genIFM(n)  ! 产生 n 阶傅里叶反变换矩阵
        integer*4 :: n
        complex*16 :: genIFM(n, n)
        genIFM = conjg(genFM(n))/n
    end function genIFM

    function Mat_Vec_FFT(A, b, N)
        integer*4 :: N
        complex*16 :: A(N, N), b(N), Mat_Vec_FFT(N)
        complex*16 :: z(2*N), v(2*N), temp(2*N)
        v(1:N) = b
        v(N+1: 2*N) = (0.0d0, 0.0d0)
        z(1:N) = A(:,1)
        z(N+1) = (0.0d0, 0.0d0)
        z(N+2:2*N) = A(1, N:2:-1)
        ! write(unit=*, fmt=*) z
        temp = matmul(genIFM(2*N), matmul(genFM(2*N), z)*matmul(genFM(2*N), v))
        Mat_Vec_FFT = temp(1:N)
    end function Mat_Vec_FFT
end module tools



program main
    use tools
    integer*4, parameter :: n = 2
    complex*16 :: a(n, n), b(n), temp(n)
    real :: start, finish
    ! 输入数据按列输入
    data a / (2.d0, 0.0d0), (1.0d0, 0.0d0), (3.d0, 0.0d0), (2.0d0, 0.0d0) /
    data b / (3.d0, 0.0d0), (7.0d0, 0.0d0) /
    call cpu_time(start) ! timer start

    temp = Mat_Vec_FFT(a, b, n)
    write(unit=*, fmt=*) matmul(a, b)
    write(unit=*, fmt=*) temp

! ////////////////////////////////////////////////////////////////////////    
    call cpu_time(finish)  ! timer stop                             ! ////
    write(unit=*, fmt=*) "cpu time: ",(finish-start)*1000, "s"      ! ////
! ////////////////////////////////////////////////////////////////////////
    stop
end program main