module commoninfo
    real*8, parameter :: pi = 3.141592653589793238d0
    real*8, parameter :: freq = 3.0d9   ! 频率
    real*8, parameter :: omega = 2*pi*freq
    real*8, parameter :: mu = 0.5d0
    real*8, parameter :: vcc = 3.0d8    ! 光速
    real*8, parameter :: b = mu*vcc/freq! 圆环半径
    real*8, parameter :: a = b/10000    ! 细线半径
    real*8, parameter :: ita = 120.d0*pi! 真空波阻抗
    real*8, parameter :: k = (2.0d0*pi/vcc)*freq    ! 波数
    complex*16 :: j = (0.0d0, 1.0d0)
end module commoninfo


module tools
    use commoninfo
    contains

    function func1(x)   ! zmn 的被积函数
        use commoninfo
        complex*16 :: func1
        real*8 :: x

    end function func1
    
    function gauss(func, a, b)  ! gauss 12 点 积分
        real*8 :: a, b
        complex*16 :: gauss
        complex*16, external :: func
        integer*4 :: i
        real*8 :: node(12), w(12), t(12)
        data node / -0.981560634246732d0, -0.904117256370452d0, -0.7699026741943177d0, -0.5873179542866143d0,&
                &   -0.3678314989981804d0, -0.12523340851114688d0, 0.12523340851114688d0, 0.3678314989981804d0, &
                &   0.5873179542866143d0, 0.7699026741943177d0, 0.904117256370452d0, 0.981560634246732d0 / 
        data w / 0.04717533638647547d0, 0.1069393259953637d0, 0.1600783285433586d0, 0.2031674267230672d0, &
                &   0.2334925365383534d0, 0.2491470458134027d0, 0.2491470458134027d0, 0.2334925365383534d0, &
                &   0.2031674267230672d0, 0.1600783285433586d0, 0.1069393259953637d0, 0.04717533638647547d0 /
        t = (b+a)/2.0d0+(b-a)*node/2.0d0
        gauss = (0.0d0, 0.0d0)
        do i = 1, 12
            gauss = gauss + func(t(i))*w(i)
        enddo
        gauss = gauss*(b-a)/2.0d0
    end function gauss

    function genZMN(fun1, m, n)
        complex*16, external :: fun1! 被积函数
        ! complex*16, external :: gf  ! 格林函数
        integer*4 :: m, n, totalN
        real*8 :: delta, phi_n, phi_m
        complex*16 :: genZMN(m, n)
        real*8 :: r1, r2, r3
        totalN = n
        delta = 2.0d0*pi/totalN
        do m = 1, totalN
            phi_m = 2.0d0*pi*m/totalN-delta/2.0d0
            do n = 1, totalN
                phi_n = 2.0d0*pi*n/totalN-delta/2.0d0
                r1 = b*dsqrt(4.0d0*dsin((phi_m-phi_n)/2.0d0)**2.0d0+(a/b)**2.0d0)
                r2 = b*dsqrt(4.0d0*dsin((phi_m-phi_n-delta)/2.0d0)**2.0d0+(a/b)**2.0d0)
                r3 = b*dsqrt(4.0d0*dsin((phi_m-phi_n+delta)/2.0d0)**2.0d0+(a/b)**2.0d0)
                genZMN(m, n) = j*k*ita*gauss(fun1, phi_n-delta/2.0d0, phi_n+delta/2.0d0)-(j*ita/k)*&
                        &(2.0d0*exp(-j*k*r1)/(4.0d0*pi*r1)-exp(-j*k*r2)/(4.0d0*pi*r2)-exp(-j*k*r3)/(4.0d0*pi*r3))
            enddo
        enddo
        write(unit=*, fmt=*) "function genZMN(m, n) SUCCED!!"
    end function genZMN

end module tools

program main
    use commoninfo
    use tools
    integer*4 :: n=4
    integer*4 :: m=4
    complex*16, allocatable :: x(:,:)
    complex*16, external :: func1
    allocate(x(n, n))
    x = genZMN( func1, m, n)
    open(11, file='text.txt')
    do i = 1, n
        write(unit=11, fmt=*) x(i, :)
    enddo
    close(11)
    write(unit=*, fmt=*) "SUCCED!!"
    stop
end program main

