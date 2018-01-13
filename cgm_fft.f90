module tools
    real*8, parameter :: pi = 3.141592653589793238d0
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

    function rcg_fft(ca, cy, cx, n)
        integer*4, parameter :: mxiter = 2000   ! 允许的最大迭代次数
        real*8 ,parameter :: err = 1.0d-16       ! 设定的容许误差
        integer*4 :: iter, n
        complex*16 :: ca(n, n), cx(n), cy(n)
        real*8 :: ak, ay, ek, qk, sk, sk2
        complex*16 :: cp(n), cprod(n), cr(n), rcg_fft(N)
        iter = 0
        ay = dble(dot_product(cy, cy))
        cprod = Mat_Vec_FFT(ca, cx, n)
        cr = cy - cprod
        cp = Mat_Vec_FFT(transpose(conjg(ca)), cr, n)
        sk = dble(dot_product(cp, cp))
        ! open(11, file='text.csv')
        do while(iter .lt. mxiter)
            cprod = Mat_Vec_FFT(ca, cp, n)
            ak = dble(dot_product(cprod, cprod))
            ak = sk/ak
            cx = cx + ak*cp
            cr = cr - ak*cprod
            ek = dble(dot_product(cr, cr))
            ! write(unit=11, fmt=*) dsqrt(ek/ay), cx
            if(dsqrt(ek/ay) .le. err) exit
            cprod = Mat_Vec_FFT(transpose(conjg(ca)), cr, n)
            sk2 = dble(dot_product(cprod, cprod))
            qk = sk2/sk
            cp = cprod + qk*cp
            sk = sk2
            iter = iter+1
            if(iter .ge. mxiter)then
                write(unit=*, fmt=*) "beyound limitation, not converged!"
                exit
            endif
        enddo
        ! close(11)
        write(unit=*, fmt=*) iter
        rcg_fft = cx
        return
    end function rcg_fft

    function rcg(ca, cy, cx, n)
        integer*4, parameter :: mxiter = 2000   ! 允许的最大迭代次数
        real*8 ,parameter :: err = 1.0d-16       ! 设定的容许误差
        integer*4 :: iter, n
        complex*16 :: ca(n, n), cx(n), cy(n)
        real*8 :: ak, ay, ek, qk, sk, sk2
        complex*16 :: cp(n), cprod(n), cr(n), rcg(N)
        iter = 0
        ay = dble(dot_product(cy, cy))
        cprod = matmul(ca, cx)
        cr = cy - cprod
        cp = matmul(transpose(conjg(ca)), cr)
        sk = dble(dot_product(cp, cp))
        ! open(11, file='text.csv')
        do while(iter .lt. mxiter)
            cprod = matmul(ca, cp)
            ak = dble(dot_product(cprod, cprod))
            ak = sk/ak
            cx = cx + ak*cp
            cr = cr - ak*cprod
            ek = dble(dot_product(cr, cr))
            ! write(unit=11, fmt=*) dsqrt(ek/ay), cx
            if(dsqrt(ek/ay) .le. err) exit
            cprod = matmul(transpose(conjg(ca)), cr)
            sk2 = dble(dot_product(cprod, cprod))
            qk = sk2/sk
            cp = cprod + qk*cp
            sk = sk2
            iter = iter+1
            if(iter .ge. mxiter)then
                write(unit=*, fmt=*) "beyound limitation, not converged!"
                exit
            endif
        enddo
        ! close(11)
        write(unit=*, fmt=*) iter
        rcg = cx
        return
    end function rcg   



end module tools

program main
    use tools
    integer*4, parameter :: n = 2
    complex*16 :: A(n, n), b(n), x0(n), x1(n), x2(n)
    real*8 :: start, finish
    call cpu_time(start)

    data A / (-0.73492d0, 7.11486d0), (0.274492957d0, 3.7885537d0), (0.274492957d0, 3.7885537d0), (-0.73492d0, 7.11486d0) /
    data b / (0.289619736d0, 0.895562183d0), (-0.28475616d0,-0.892163111d0) /
    data x0 / (0.0d0, 0.0d0), (0.0d0, 0.0d0) /
    

    x1 = rcg_fft(A, b, x0, n)
    x2 = rcg(A, b, x0, n)
    ! do i = 1, n
    !     write(unit=*, fmt=*) x(i)
    ! enddo
    write(unit=*, fmt=*) x1
    write(unit=*, fmt=*) x2
    call cpu_time(finish)
    write(unit=*, fmt=*) "time elapsed: ", (finish-start), 'ms'
    
    stop
end program main
