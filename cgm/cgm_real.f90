! 共轭梯度法的计算
module conj
!--------------------------------------
!   参数:   1: 最大允许迭代次数
!           2: tol--误差限
    integer*4, parameter :: IMAX = 200
    real*8 :: tol = 1.0d-7
    contains
    function cgm(A, b, x0, N)
        integer*4 :: N, i, j, k
        real*8 :: A(N, N), b(N), cgm(N),x0(N)
        real*8 :: r0(N), r1(N), p0(N), p1(N)
        real*8 :: x1(N), x2(N)
        real*8 :: temp1, temp2, afa, dr_s, temp3, beta
        ! write(unit=*, fmt=*) "共轭梯度"
        x1 = x0
        r0 = b- Ar(A, x1, N)
        p0 = r0
        do k = 1, IMAX
            temp1 = dr(r0, N)
            temp2 = rAr(A, p0, N)
            afa = temp1/temp2
            x2 = x1+afa*p0
            ! write(unit=*, fmt=*) k, x2 ! 迭代中间值
            ! 如果 r0 接近, 则停止迭代
            dr_s = dsqrt(dr(r0, N))
            if(dr_s .lt. tol)exit
            r1 = r0-afa*Ar(A, p0, N)
            temp3 = dr(r1, N)
            beta = temp3/temp1
            p1 = r1 + beta*p0
            r0 = r1 
            p0 = p1 
            x1 = x2 
        enddo
    cgm = x2
    end function cgm

    function dr(r, N)
    ! 计算向量长度平方
        integer*4 :: N, i
        real*8 :: r(N), dr, s
        s = 0.0d0
        do i = 1, N
            s = s+r(i)**2
        enddo
    end function dr

    function Ar(A, r, N)
        integer*4 :: N, i, j
        real*8 :: A(N, N), r(N), temp(N), Ar(N)
        temp = 0
        do i = 1, N
            do j = 1, N
                temp(i) = temp(i) + A(i, j)*r(j)
            enddo
        enddo
        Ar = temp
    end function Ar

    function rAr(A, r, N)
        integer*4 :: i, N
        real*8 :: A(N, N), r(N), temp(N)
        temp = Ar(A, r, N)
        rAr = dot_product(r, temp)
    end function rAr
end module conj

program main
    use conj
    implicit real*8(a-z)
    integer*4, parameter :: N=4
    real*8 :: A(N, N), b(N), x0(N), x(N)
    data A / 5.0d0,7.0d0,6.0d0,5.0d0,7.0d0,10.0d0,8.0d0,7.0d0,6.0d0,8.0d0,10.0d0,9.0d0,5.0d0,7.0d0,9.0d0,10.0d0 /
    x0=(/0.0d0, 0.0d0, 0.0d0, 0.0d0/)
    b = (/62.0d0, 87.0d0, 91.0d0, 90.0d0/)
    A = reshape(A,(/4,4/))
    x = cgm(A, b, x0, N)
    write(unit=*, fmt=*) x
    ! write(unit=*, fmt=*) Ar(A, x, N)-b
    
    stop
end program main
