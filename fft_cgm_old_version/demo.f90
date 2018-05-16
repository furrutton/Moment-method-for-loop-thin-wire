module constants

  real*8, parameter :: pi = 3.141592653589793238d0
  real*8, parameter :: freq = 3.0d9   ! 频率
  real*8, parameter :: omega = 2*pi*freq ! 角频率
  real*8, parameter :: vcc = 3.0d8    ! 光速
  real*8, parameter :: lambda = vcc/freq ! 波长
  real*8, parameter :: zhouchang = 400.0d0*lambda     ! 圆环周长
  real*8, parameter :: b = zhouchang/pi/2.0d0! 圆环半径
  real*8, parameter :: a = b/10000    ! 细线半径 假设远远小于细线环半径
  real*8, parameter :: ita = 120.d0*pi! 真空波阻抗
  real*8, parameter :: k = 2.0d0*pi/lambda    ! 波数
  complex*16, parameter :: j = (0.0d0, 1.0d0) ! 虚数单位
  integer*4, parameter :: n = 1024
  integer*4, parameter :: m = n
  real*8, parameter :: delta = 2.0d0*pi/n

end module constants

module tools
  use constants
contains
  function f1(x) result(ans)
    implicit none
    complex*16 :: ans
    real*8 :: r, x
    r = b*dsqrt(4.0d0*dsin(x/2.0d0)**2.0d0+(a/b)**2.0d0)
    ans = exp(-j*k*r)*dcos(x)/(4.0d0*pi*r)
  end function f1

  function Ei(phi) result(ans)
    implicit none
    real*8 :: phi
    complex*16 :: ans
    ans = dcos(phi)*exp(-j*k*b*dcos(phi))
  end function Ei

  function Gauss1d(func, l, u) result(gauss)
    implicit none
    complex*16 :: gauss
    complex*16, external :: func
    integer*4 :: i
    real*8 :: node(12), w(12), t(12)
    real*8 :: l, u
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
  end function Gauss1d

  function genZMN(func) result(ans)
    implicit none
    complex*16 :: ans(m, n)
    complex*16, external :: func
    real*8 :: phi_n, phi_m
    real*8 :: r1, r2, r3, l, u
    complex*16 :: t1, t2, t3
    integer*4 :: rows, cols
    do rows = 1, m, 1
      phi_m = 2.0d0*pi*rows/m-delta/2.0d0
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
        ans(rows, cols) = j*k*ita*Gauss1d(func, l, u)-(j*ita/k)*(t1-t2-t3)
      end do
    end do
  end function genZMN

  function genB(func) result(ans)
    implicit none
    complex*16, external :: func
    integer*4 :: i
    real*8 :: l, u, phi_i
    complex*16 :: ans(n)
    do i = 1, n, 1
      phi_i = 2.0d0*pi*i/n-delta/2.0d0
      l = phi_i-delta/2.0d0
      u = phi_i+delta/2.0d0
      ans(i) = Gauss1d(func, l, u)
    end do
  end function genB

  function genFM(dems) result(fft_matrix)
    implicit none
    integer*4 :: dems
    complex*16 :: fft_matrix(dems, dems)
    integer*4 :: cols, rows
    do rows = 1, dems, 1
      do cols = 1, dems, 1
        fft_matrix(rows, cols) = exp(-2*pi*j*(rows-1)*(cols-1)/dems)
      end do
    end do
  end function genFM

  function genIFM(dems) result(ifft_matrix)
    implicit none
    integer*4 :: dems
    complex*16 :: ifft_matrix(dems, dems)
    ifft_matrix = conjg(genFM(dems))/dems
  end function genIFM

  function Mat_Vec_FFT(A, b) result(vector)
    implicit none
    complex*16 :: vector(n)
    complex*16 :: A(n, n), b(n)
    complex*16 :: z(2*n), v(2*n), tmp(2*n)
    v(1: n) = b
    v(n+1: 2*n) = (0.0d0, 0.0d0)
    z(1: n)     = A(:,1)
    z(n+1)      = (0.0d0, 0.0d0)
    z(n+2: 2*n) = A(1, n:2:-1)
    tmp = matmul(genIFM(2*n), matmul(genFM(2*n), z)*matmul(genFM(2*n), v))
    vector = tmp(1: n)
  end function Mat_Vec_FFT

  function rcg(ca, cy, cx) result(vector)
    implicit none
    complex*16 :: vector(N) ! 返回值
    integer*4 :: mxiter = n ! 允许的最大迭代次数
    real*8 :: err = 1.0d-16
    integer*4 :: iter
    complex*16 :: ca(n, n), cx(n), cy(n)  ! ca->A    cx->未知向量    cy->右侧系数   Ax=y
    real*8 :: ak, ay, ek, qk, sk, sk2       ! 中间变量
    complex*16 :: cp(n), cprod(n), cr(n)    ! 中间变量
    ay = dble(dot_product(cy, cy))
    cprod = Mat_Vec_FFT(ca, cx)          ! 使用 fft 做矩阵向量积, ca是输入Z矩阵(Toeplitz), cx是n维列向量
    cr = cy - cprod
    cp = Mat_Vec_FFT(transpose(conjg(ca)), cr)
    sk = dble(dot_product(cp, cp))
    do iter = 1, mxiter, 1
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
    end do
    if(iter .eq. mxiter) write(unit=*, fmt=*) "beyound limitation, not converged!"
    write(unit=*, fmt=*) "converged step by: ", iter ! 显示收敛步数
    vector = cx
  end function rcg

end module tools




program main
  use constants
  use tools
  implicit none

  integer*4 :: i
  complex*16 :: zmn(n, n), rb(n), an(n), x0(n)
  real*8 :: start, finish
  call cpu_time(start)

  x0 = (0.0d0, 0.0d0)
  zmn = genZMN(f1)
  rb = genB(Ei)
  an = rcg(zmn, rb, x0)

  open(11, file='Zmn.txt')
  open(12, file='RightB.txt')
  open(13, file='An.txt')
  open(14, file='An_dB.txt')      ! out put current expanssion by dB

  do i = 1, n
    write(unit=11, fmt=*) Zmn(i, :)
    write(unit=12, fmt=*) rB(i)
    write(unit=13, fmt=*) An(i)
    write(unit=14, fmt=*) 10.0*dlog10(dsqrt(dble(An(i)*conjg(An(i)))))
  enddo

  close(11)
  close(12)
  close(13)
  close(14)


  call cpu_time(finish)

  write(unit=*, fmt=*) "time elapsed:", finish-start
  stop
end program main
