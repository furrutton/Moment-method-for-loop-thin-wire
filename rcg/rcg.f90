module tools
    contains
    function rcg(ca, cy, cx, n)
        integer*4, parameter :: mxiter = 2000   ! 允许的最大迭代次数
        real*8 ,parameter :: err = 1.0d-7       ! 设定的容许误差
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
        rcg = cx
        return
    end function rcg    

end module tools


program main
    use tools
    integer*4, parameter :: n = 2
    real*8 :: start, finish
    complex*16 :: A(n, n), b(n), x(n)
    data A / (1.0d0, 1.0d0), (2.0d0, 2.0d0), (2.0d0, 1.0d0), (3.0d0, 4.0d0) /
    data b / (5.0d0, 3.0d0), (7.0d0, 8.0d0) /
    data x / (0.0d0, 0.0d0), (0.0d0, 0.0d0) /

    call cpu_time(start)
    x = rcg(A, b, x, n)
    call cpu_time(finish)
    write(unit=*, fmt=*) x
    write(unit=*, fmt=*) "time elapsed : ", (finish-start), 'ms'
    
    stop
end program main


