module commoninfo
    implicit none
    real*8, parameter :: pi = 3.141592653589793238d0
    real*8, parameter :: freq = 3.0d9   ! 频率
    real*8, parameter :: omega = 2*pi*freq ! 角频率
    real*8, parameter :: vcc = 3.0d8    ! 光速
    real*8, parameter :: lambda = vcc/freq ! 波长
    real*8, parameter :: zhouchang = 100.0d0*lambda     ! 圆环周长
    real*8, parameter :: b = zhouchang/pi/2.0d0! 圆环半径
    real*8, parameter :: a = b/10000    ! 细线半径 假设远远小于细线环半径
    real*8, parameter :: ita = 120.d0*pi! 真空波阻抗
    real*8, parameter :: k = 2.0d0*pi/lambda    ! 波数
    complex*16 :: j = (0.0d0, 1.0d0)
end module commoninfo


module tools
    use commoninfo
    implicit none
    contains
    
    function func1(x)   ! zmn 的被积函数
        use commoninfo
        implicit none
        complex*16 :: func1
        real*8 :: r, x  ! x 代表 phim - phi'
        r = b*dsqrt(4.0d0*dsin(x/2.0d0)**2.0d0+(a/b)**2.0d0)
        func1 = exp(-j*k*r)*dcos(x)/(4.0d0*pi*r)
    end function func1

    function Ei(phi)    ! 右端向量被积表达式
        use commoninfo
        implicit none
        real*8 :: phi
        complex*16 :: Ei
        Ei = dcos(phi)*exp(-j*k*b*dcos(phi))
    end function Ei

    function gauss(func1, l, u)  ! gauss 12 点 积分
        use commoninfo
        implicit none
        real*8 :: l, u
        complex*16 :: gauss
        complex*16, external :: func1
        integer*4 :: i
        real*8 :: node(12), w(12), t(12)
        data node / -0.981560634246732d0, -0.904117256370452d0, -0.7699026741943177d0, -0.5873179542866143d0,&
                &   -0.3678314989981804d0, -0.12523340851114688d0, 0.12523340851114688d0, 0.3678314989981804d0, &
                &   0.5873179542866143d0, 0.7699026741943177d0, 0.904117256370452d0, 0.981560634246732d0 / 
        data w / 0.04717533638647547d0, 0.1069393259953637d0, 0.1600783285433586d0, 0.2031674267230672d0, &
                &   0.2334925365383534d0, 0.2491470458134027d0, 0.2491470458134027d0, 0.2334925365383534d0, &
                &   0.2031674267230672d0, 0.1600783285433586d0, 0.1069393259953637d0, 0.04717533638647547d0 /
        t = (u+l)/2.0d0+(u-l)*node/2.0d0
        gauss = (0.0d0, 0.0d0)
        do i = 1, 12
            gauss = gauss + func1(t(i))*w(i)
        enddo
        gauss = gauss*(u-l)/2.0d0
    end function gauss

    function genZMN(func1, m, n)
        use commoninfo
        implicit none
        integer*4 :: m, n, rows, cols
        complex*16, external :: func1! 被积函数
        real*8 :: delta, phi_n, phi_m
        complex*16 :: genZMN(m, n)
        real*8 :: r1, r2, r3, lb, ub
        complex*16 :: t1, t2, t3
        delta = 2.0d0*pi/n
        do rows = 1, m
            phi_m = 2.0d0*pi*rows/m-delta/2.0d0
            do cols = 1, n
                phi_n = 2.0d0*pi*cols/n-delta/2.0d0
                r1 = b*dsqrt(4.0d0*dsin((phi_m-phi_n)/2.0d0)**2.0d0+(a/b)**2.0d0)
                r2 = b*dsqrt(4.0d0*dsin((phi_m-phi_n-delta)/2.0d0)**2.0d0+(a/b)**2.0d0)
                r3 = b*dsqrt(4.0d0*dsin((phi_m-phi_n+delta)/2.0d0)**2.0d0+(a/b)**2.0d0)
                t1 = 2.0d0*exp(-j*k*r1)/(4.0d0*pi*r1)
                t2 = exp(-j*k*r2)/(4.0d0*pi*r2)
                t3 = exp(-j*k*r3)/(4.0d0*pi*r3)
                lb = phi_n-phi_m-delta/2.0d0
                ub = phi_n-phi_m+delta/2.0d0
                genZMN(rows, cols) = j*k*ita*gauss(func1, lb, ub)-(j*ita/k)*(t1-t2-t3)
            enddo
        enddo
        write(unit=*, fmt=*) "function genZMN(m, n) SUCCED!!"
    end function genZMN

    function genB(Ei, n)    ! 得到右端向量
        use commoninfo
        implicit none
        integer*4 :: n, i
        complex*16, external :: Ei
        complex*16 :: genB(n)
        real*8 :: lb, ub, phi_i, delta
        delta = 2.0d0*pi/n
        do i = 1, n
            phi_i = 2.0d0*pi*i/n-delta/2.0d0
            lb = phi_i-delta/2.0d0
            ub = phi_i+delta/2.0d0
            genB(i) = gauss(Ei, lb, ub)
        enddo
    end function genB

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
    use commoninfo
    use tools
    implicit none
    integer*4, parameter :: n=8
    integer*4, parameter :: m=8
    integer*4 :: i
    complex*16, allocatable :: Zmn(:,:), rB(:), An(:), x0(:)
    real*8 :: start, finish
    call cpu_time(start)

    allocate(Zmn(n, n))
    allocate(rB(n))
    allocate(An(n))
    allocate(x0(n))
    x0 = (0.0d0, 0.0d0)
    Zmn = genZMN( func1, m, n)
    rB = genB(Ei, n)
    An = rcg_fft(Zmn, rB, x0, n)

    open(11, file='Zmn.txt')
    open(12, file='RightB.txt')
    open(13, file='An.txt')

    do i = 1, n
        write(unit=11, fmt=*) Zmn(i, :)
        write(unit=12, fmt=*) rB(i)
        write(unit=13, fmt=*) An(i)
    enddo

    close(11)
    close(12)
    close(13)

    call cpu_time(finish)
    write(unit=*, fmt=*) "time elapsed: ", finish-start
    stop
end program main
