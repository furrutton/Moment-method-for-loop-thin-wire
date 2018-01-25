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
    complex*16, parameter :: j = (0.0d0, 1.0d0) ! 虚数单位
end module commoninfo

module tools
    use commoninfo
    implicit none
    contains

! ---------------------------- zmn 的被积函数 ----------------------------
    function func1(x)   
        complex*16 :: func1
        real*8 :: r, x  ! x 代表 phim - phi'
        r = b*dsqrt(4.0d0*dsin(x/2.0d0)**2.0d0+(a/b)**2.0d0)
        func1 = exp(-j*k*r)*dcos(x)/(4.0d0*pi*r)
    end function func1

! -------------------------- 右端向量被积表达式 --------------------------
    function Ei(phi)
        real*8 :: phi
        complex*16 :: Ei
        Ei = dcos(phi)*exp(-j*k*b*dcos(phi))
    end function Ei

! ------------------------- gauss 12 点积分 -----------------------------
    function gauss(func1, l, u)
        real*8 :: l, u                  ! 积分 下, 上 限
        complex*16 :: gauss             ! 积分返回值
        complex*16, external :: func1   ! 被积函数
        integer*4 :: i
        real*8 :: node(12), w(12), t(12)
        data node / -0.981560634246732d0, -0.904117256370452d0, -0.7699026741943177d0, -0.5873179542866143d0,   &
                &   -0.3678314989981804d0, -0.12523340851114688d0, 0.12523340851114688d0, 0.3678314989981804d0, &
                &   0.5873179542866143d0, 0.7699026741943177d0, 0.904117256370452d0, 0.981560634246732d0 / 
        data w / 0.04717533638647547d0, 0.1069393259953637d0, 0.1600783285433586d0, 0.2031674267230672d0,       &
                &   0.2334925365383534d0, 0.2491470458134027d0, 0.2491470458134027d0, 0.2334925365383534d0,     &
                &   0.2031674267230672d0, 0.1600783285433586d0, 0.1069393259953637d0, 0.04717533638647547d0 /
        t = (u+l)/2.0d0+(u-l)*node/2.0d0
        gauss = (0.0d0, 0.0d0)
        do i = 1, 12
            gauss = gauss + func1(t(i))*w(i)
        enddo
        gauss = gauss*(u-l)/2.0d0
    end function gauss

! ----------------------------- 产生 n 维 Z 矩阵 ---------------------------------
    function genZMN(func1, m, n)
        integer*4 :: m, n, rows, cols
        complex*16, external :: func1   ! 被积函数
        real*8 :: delta, phi_n, phi_m   ! delta: 分段间隔
        complex*16 :: genZMN(m, n)      ! 函数返回值: 数组, m×n阶
        real*8 :: r1, r2, r3, lb, ub    ! r1, r2, r3: 计算的中间量  lb, ub: 积分下上限
        complex*16 :: t1, t2, t3        ! 中间量
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

! --------------------------- 产生 n 维右端向量 -----------------------------------
    function genB(Ei, n)
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
        write(unit=*, fmt=*) "function genB(n) SUCCED!!"
    end function genB

! ------------------------- 产生 N 阶傅里叶矩阵 ---------------------------
    function genFM(n)
        integer*4 :: n, cols, rows
        complex*16 :: genFM(n, n)
        do rows = 1, n
            do cols = 1, n
                genFM(rows, cols) = exp(-2*pi*j*(rows-1)*(cols-1)/n)
            enddo
        enddo
    end function genFM

! ------------------------ 产生 N 阶傅里叶反变换矩阵 -----------------------
    function genIFM(n)
        integer*4 :: n
        complex*16 :: genIFM(n, n)
        genIFM = conjg(genFM(n))/n
    end function genIFM

! ----------------------- 矩阵向量积 with fft & ifft --------------------------
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
        temp = matmul(genIFM(2*N), matmul(genFM(2*N), z)*matmul(genFM(2*N), v)) ! MVF = IFFT(FFT(Z)*FFT(v))
        Mat_Vec_FFT = temp(1:N)
    end function Mat_Vec_FFT

! --------------------------- cg with fft & ifft method ----------------------------
    function rcg_fft(ca, cy, cx, n)
        integer*4, parameter :: mxiter = 200    ! 允许的最大迭代次数
        real*8 ,parameter :: err = 1.0d-16      ! 设定的容许误差
        integer*4 :: iter, n                    ! 迭代次数, 秩
        complex*16 :: ca(n, n), cx(n), cy(n)    ! ca->阻抗矩阵    cx->未知向量    cy->右侧系数   Ax=y
        real*8 :: ak, ay, ek, qk, sk, sk2       ! 中间变量
        complex*16 :: cp(n), cprod(n), cr(n), rcg_fft(N)    ! 中间变量
        iter = 0    ! 迭代计数置零
        ay = dble(dot_product(cy, cy))
        cprod = Mat_Vec_FFT(ca, cx, n)          ! 使用 fft 做矩阵向量积, ca是输入Z矩阵(Toeplitz), cx是n维列向量
        cr = cy - cprod
        cp = Mat_Vec_FFT(transpose(conjg(ca)), cr, n)
        sk = dble(dot_product(cp, cp))
        ! open(11, file='text.csv')
        do while(iter .lt. mxiter)              ! 迭代开始, 条件小于最大迭代次数
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
            if(iter .eq. mxiter)then            ! 当等于最大迭代次数的时候 输出不收敛
                write(unit=*, fmt=*) "beyound limitation, not converged!"
                exit
            endif
        enddo
        ! close(11)
        write(unit=*, fmt=*) "converged step by: ", iter    ! 显示迭代次数
        rcg_fft = cx    ! 返回程序运算结果
        write(unit=*, fmt=*) "function rcg_fft(n) SUCCED!!"
    end function rcg_fft
end module tools

program main
    use commoninfo
    use tools
    implicit none
    integer*4, parameter :: n=128
    integer*4, parameter :: m=128
    integer*4 :: i
    complex*16, allocatable :: Zmn(:,:), rB(:), An(:), x0(:)
    real*8 :: start, finish
    call cpu_time(start)    ! 计时开始

! ----------------- 分配空间 -----------------
    allocate(Zmn(n, n))             ! Zmn为阻抗矩阵            
    allocate(rB(n))                 ! rB 表示右端向量
    allocate(An(n))                 ! An 是展开系数
    allocate(x0(n))                 ! 给定初始的迭代向量

! ----------------------- 计算过程 -------------------------
    x0 = (0.0d0, 0.0d0)             
    Zmn = genZMN(func1, m, n)       
    rB = genB(Ei, n)                
    An = rcg_fft(Zmn, rB, x0, n)    
    ! write(unit=*, fmt=*) matmul(Zmn, An)-rB   ! 在-15次以内

! ----------------- 存储计算结果 ------------------------   
    open(11, file='Zmn.txt')
    open(12, file='RightB.txt')
    open(13, file='An.txt')
    open(14, file='An_dB.txt')      ! out put current expanssion by dB
    do i = 1, n
        write(unit=11, fmt=*) Zmn(i, :)
        write(unit=12, fmt=*) rB(i)
        write(unit=13, fmt=*) An(i)
        write(unit=14, fmt=*) 20.0*dlog10(dsqrt(dble(An(i)*conjg(An(i)))))
    enddo
    close(11)
    close(12)
    close(13)
    close(14)
! ------------------------------------------------------

    call cpu_time(finish)       ! 计时结束
    write(unit=*, fmt=*) "time elapsed: ", finish-start ! 输出 cpu 耗时
    stop
end program main
