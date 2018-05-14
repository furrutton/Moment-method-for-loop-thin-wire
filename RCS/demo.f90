! ----------------------------------------------------------------------------------
! 常量定义
! -----------------------------------------------------------------------------------
module constants
  real*8, parameter       :: pi = 3.141592653589793238d0
  real*8, parameter       :: freq = 3.0d9                   ! 频率
  real*8, parameter       :: omega = 2*pi*freq              ! 角频率
  real*8, parameter       :: vcc = 3.0d8                    ! 光速
  real*8, parameter       :: lambda = vcc/freq              ! 波长

  real*8, parameter       :: zhouchang = lambda*100.0d0     ! 圆环周长 >>>>> >>>>>>>>>>波长的倍数

  real*8, parameter       :: b = zhouchang/pi/2.0d0         ! 圆环半径
  real*8, parameter       :: a = b/10000                    ! 细线半径 假设远远小于细线环半径
  real*8, parameter       :: ita = 120.d0*pi                ! 真空波阻抗
  real*8, parameter       :: k = 2.0d0*pi/lambda            ! 波数
  complex*16, parameter   :: j = (0.0d0, 1.0d0)             ! 虚数单位

  integer*4, parameter    :: n = 2**11                      ! 分段数 >>>>>>>>>>>>>>>>>>>>>> n

  real*8, parameter       :: delta = 2.0d0*pi/n             ! 分段间隔 delta
end module constants
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
! 计算过程中使用到的函数合集
! -----------------------------------------------------------------------------------
module tools
  use constants

contains
  ! -----------------------------------
  ! 定义计算 阻抗矩阵 时需要使用的被积函数
  ! -----------------------------------
  function f1(x) result(ans)
    implicit none
    complex*16  :: ans
    real*8      :: r, x

    r = b*dsqrt(4.0d0*dsin(x/2.0d0)**2.0d0+(a/b)**2.0d0)
    ans = exp(-j*k*r)*dcos(x)/(4.0d0*pi*r)

    return
  end function f1

  ! -----------------------------------
  ! 定义计算 右端向量 时需要使用的入射场被积函数
  ! -----------------------------------
  function Ei(phi) result(ans)
    implicit none
    real*8      :: phi
    complex*16  :: ans

    ans = dcos(phi)*exp(-j*k*b*dcos(phi))*b

    return
  end function Ei

  ! ------------------------------------------
  ! 定义 一维高斯积分 函数
  ! ------------------------------------------
  function Gauss1d(func, l, u) result(gauss)
    implicit none
    complex*16            :: gauss
    complex*16, external  :: func
    integer*4             :: i
    real*8                :: node(12), w(12), t(12)
    real*8                :: l, u

    data node / -0.981560634246732d0, -0.904117256370452d0, -0.7699026741943177d0, -0.5873179542866143d0,   &
            &   -0.3678314989981804d0, -0.12523340851114688d0, 0.12523340851114688d0, 0.3678314989981804d0, &
            &   0.5873179542866143d0, 0.7699026741943177d0, 0.904117256370452d0, 0.981560634246732d0 /
    data w / 0.04717533638647547d0, 0.1069393259953637d0, 0.1600783285433586d0, 0.2031674267230672d0,       &
            &   0.2334925365383534d0, 0.2491470458134027d0, 0.2491470458134027d0, 0.2334925365383534d0,     &
            &   0.2031674267230672d0, 0.1600783285433586d0, 0.1069393259953637d0, 0.04717533638647547d0 /
    t = (u+l)/2.0d0+(u-l)*node/2.0d0
    gauss = (0.0d0, 0.0d0)
    do i = 1, 12, 1
      gauss = gauss + func(t(i))*w(i)
    end do
    gauss = gauss*(u-l)/2.0d0

    return
  end function Gauss1d

  ! ------------------------------------------
  ! 定义 产生阻抗矩阵 的函数
  ! ------------------------------------------
  function Generate_Zmn(func) result(ans)
    implicit none
    complex*16            :: ans(n, n)
    complex*16, external  :: func
    real*8                :: phi_n, phi_m
    real*8                :: r1, r2, r3, l, u
    complex*16            :: t1, t2, t3
    integer*4             :: rows, cols

    do rows = 1, n, 1
      phi_m = 2.0d0*pi*rows/n-delta/2.0d0
      do cols = 1, n, 1
        phi_n = 2.0d0*pi*cols/n-delta/2.0d0
        r1 = b*dsqrt(4.0d0*dsin((phi_m-phi_n)/2.0d0)**2.0d0+(a/b)**2.0d0)
        r2 = b*dsqrt(4.0d0*dsin((phi_m-phi_n-delta)/2.0d0)**2.0d0+(a/b)**2.0d0)
        r3 = b*dsqrt(4.0d0*dsin((phi_m-phi_n+delta)/2.0d0)**2.0d0+(a/b)**2.0d0)
        t1 = 2.0d0*exp(-j*k*r1)/(4.0d0*pi*r1)
        t2 = exp(-j*k*r2)/(4.0d0*pi*r2)
        t3 = exp(-j*k*r3)/(4.0d0*pi*r3)
        l = phi_n-phi_m-delta/2.0d0
        u = phi_n-phi_m+delta/2.0d0
        ans(rows, cols) = -j*k*ita*b**2*delta*Gauss1d(func, l, u)-(j*ita/k)*(t1-t2-t3)
      end do
    end do

    return
  end function Generate_Zmn

  ! ------------------------------------------
  ! 定义 产生右端向量 的函数
  ! ------------------------------------------
  function Generate_B(func) result(ans)
    implicit none
    complex*16, external  :: func
    integer*4             :: i
    real*8                :: l, u, phi_i
    complex*16            :: ans(n)

    do i = 1, n, 1
      phi_i = 2.0d0*pi*i/n-delta/2.0d0
      l = phi_i-delta/2.0d0
      u = phi_i+delta/2.0d0
      ans(i) = Gauss1d(func, l, u)
    end do

    return
  end function Generate_B

  ! -----------------------------------------------------------
  ! 定义 计算快速傅里叶变换 子过程 递归调用
  !   使用 库利－图基(Cooley-Tukey) 快速傅里叶变换算法
  !   参考网页 https://en.wikipedia.org/wiki/Cooley-Tukey_FFT_algorithm
  ! -----------------------------------------------------------
  recursive subroutine fft(x)
    implicit none
    complex*16, dimension(:), intent(inout) :: x                  ! 一维, 双向数据参数
    complex*16, allocatable                 :: even(:), odd(:)    ! 一维, 可变数组
    complex*16                              :: t                  ! 临时变量, fft因子
    integer*4                               :: n, i               ! 数组长度 和 迭代器

    n = size(x)
    if ( n .le. 1 ) return

    allocate(odd((n)/2))
    allocate(even(n/2))


    odd  = x(1:n:2)
    even = x(2:n:2)

    call fft(odd)
    call fft(even)


    do i = 1, n/2, 1
      t        = exp(-2.0d0*j*pi*(i-1)/n)*even(i)
      x(i)     = odd(i)+t
      x(i+n/2) = odd(i)-t
    end do

    deallocate(odd)
    deallocate(even)

    return
  end subroutine fft

  ! ------------------------------------------
  ! 定义 逆快速傅里叶变换 子过程
  !   使用 fft 子过程求共轭
  ! ------------------------------------------
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

  ! ------------------------------------------
  ! 定义 矩阵向量乘积 函数
  !   使用 快速傅里叶变换 子过程
  ! ------------------------------------------
  function Mat_Vec_FFT(A, b) result(vector)
    implicit none
    complex*16 :: vector(n)
    complex*16 :: A(n, n), b(n)
    complex*16 :: z(n), v(n), tmp(n)

    z = A(:,1)
    v = b

    call fft(z)
    call fft(v)
    tmp = z*v
    call ifft(tmp)
    vector = tmp

    return
  end function Mat_Vec_FFT

  ! ------------------------------------------
  ! 定义 复共轭梯度 函数
  !   凡遇到 Toeplitz矩阵与任意向量做内积 均使用 fft 和 ifft 实现
  ! ------------------------------------------
  function rcg(ca, cy, cx) result(vector)
    implicit none
    complex*16    :: vector(N)                ! 返回值
    integer*4     :: mxiter = n               ! 允许的最大迭代次数
    real*8        :: err = 1.0d-7
    integer*4     :: iter
    complex*16    :: ca(n, n), cx(n), cy(n)   ! ca->A    cx->未知向量    cy->右侧系数   Ax=y
    real*8        :: ak, ay, ek, qk, sk, sk2  ! 中间变量
    complex*16    :: cp(n), cprod(n), cr(n)   ! 中间变量

    ay = dble(dot_product(cy, cy))
    cprod = Mat_Vec_FFT(ca, cx)           ! 使用 fft 做矩阵向量积, ca是输入Z矩阵(Toeplitz), cx是n维列向量
    cr = cy - cprod
    cp = Mat_Vec_FFT(transpose(conjg(ca)), cr)
    sk = dble(dot_product(cp, cp))
    do iter = 0, mxiter, 1
      cprod = Mat_Vec_FFT(ca, cp)
      ak = dble(dot_product(cprod, cprod))
      ak = sk/ak
      cx = cx + ak*cp
      cr = cr - ak*cprod
      ek = dble(dot_product(cr, cr))
      if(dsqrt(ek/ay) .le. err) exit
      cprod = Mat_Vec_FFT(transpose(conjg(ca)), cr)
      sk2 = dble(dot_product(cprod, cprod))
      qk = sk2/sk
      cp = cprod + qk*cp
      sk = sk2
      write(unit=*, fmt=*) iter
    end do
    if(iter .ge. mxiter) write(unit=*, fmt=*) "beyound limitation, not converged!"
    write(unit=*, fmt=*) "converged step by: ", iter ! 显示收敛步数
    vector = cx

    return
  end function rcg

end module tools
! -----------------------------------------------------------------------------------


! --------------------- 主程序 -----------------------
program main
  use constants
  use tools

  implicit none
  integer*4   :: i
  complex*16  :: zmn(n, n), rb(n), an(n), x0(n)
  real*8      :: start, finish

  ! ------------------------------------------------------
  ! 判断 n 的取值是否为 2 的幂次
  ! 如果不是 2 的指数次则提示错误退出程序
  ! ------------------------------------------------------
    if ( mod(n, 2) .ne. 0 ) then
      write(unit=*, fmt=*) "N SET ERROR!!!"
      goto 100
    end if
  ! ------------------------------------------------------

  call cpu_time(start)            ! 开始计时

  x0  = (0.0d0, 0.0d0)
  zmn = Generate_Zmn(f1)
  rb  = Generate_B(Ei)
  an  = rcg(zmn, rb, x0)

  ! open(11, file='Zmn100.txt')
  open(12, file='RightB100.txt')
  open(13, file='An100.txt')
  open(14, file='An_dB100.txt')      ! out put current expanssion by dB
  open(15, file='time.csv', position='APPEND')
  open(16, file='rb.txt')

  !$OMP PARALLEL DO
  do i = 1, n
    ! write(unit=11, fmt=*) Zmn(i, :)
    write(unit=12, fmt=*) rB(i)
    write(unit=13, fmt=*) An(i)
    write(unit=14, fmt=*) 10.0*dlog10(dsqrt(dble(An(i)*conjg(An(i)))))
    write(unit=16, fmt=*) real(rb(i)), imag(rb(i))
  enddo
  !$OMP END PARALLEL DO
  ! close(11)
  close(12)
  close(13)
  close(14)
  close(15)
  close(16)

  call cpu_time(finish)            ! 计时结束

  write(unit=*, fmt=*) "time elapsed:", finish-start
  write(unit=15, fmt=*) finish-start
100  stop
end program main
