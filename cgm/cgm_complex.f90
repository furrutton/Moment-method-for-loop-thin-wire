module cgm
    contains
    function bicgm(A, b, x0, N)
        integer*4 :: N, i, j, k, maxit
        complex*16 :: A(N, N)
        complex*16 :: b(N), x0(N)
        complex*16 :: x0s(N), p(N), ps(N)
        complex*16 :: x(N), xs(N)
        complex*16 :: alpha, rho0, rho1, r_rs
        complex*16 :: r(N), rs(N), res_vec(N)
        complex*16 :: Ax(N), ATx(N)
        complex*16 :: Ap(N), Aps(N)
        complex*16 :: beta
        real*8 :: tol, res, n2b, n2r0, rel_res
        complex*16 :: bicgm(N)

        x0s = conjg(x0)
        Ax = matmul(A, x0)
        ! ATx = matmul(transpose(conjg(A)),x0s)
        ATx = matmul(transpose(conjg(A)),x0s)
        x = x0

        r = b - Ax
        rs = conjg(b) - ATx
        p = r
        ps = rs

        n2b = dsqrt(dble(dot_product(b, b)))
        res = dsqrt(dble(dot_product(r, r)))/n2b

        n2r0 = dble(dot_product(r, rs))

        if(n2r0 .eq. 0.0d0) then
            res = 1.0d-20
            ! write(unit=*, fmt=*) "(r, rs) == 0 !!!"
        endif
        ! write(unit=*, fmt=*) "n2r0 = ", n2r0

        tol = 1.0d-6
        maxit = 200
        k = 0
        
        ! main loop
        open(11, file='residul.txt')
        do while (res .gt. tol .and. k .lt. maxit)
            k = k+1
            Ap = matmul(A, p)
            Aps = matmul(conjg(A), ps)
            rho0 = dot_product(r, rs)
            rho1 = dot_product(Ap, ps)
            alpha = rho0/rho1
            x = x + alpha*p
            r = r - alpha*Ap
            rs = rs-alpha*Aps
            r_rs = dot_product(r, rs) 
            beta = r_rs/rho0
            p = r + beta*p
            ps = rs + beta*ps
            Ax = matmul(A, x)
            res_vec = b-Ax 
            res = dsqrt(dble(dot_product(res_vec, res_vec)))/n2b
            write(unit=11, fmt=*) res
        enddo
        close(11)
        
        if(k .lt. maxit)then
            write(unit=*, fmt=*) "Converged in", k, "iterations"
        else
            write(unit=*, fmt=*) "Stopped after", k, "due to over-limitation !!"
        endif     
        bicgm = x
    end function bicgm
end module cgm

program main
    use cgm
    integer*4, parameter :: n = 2
    complex*16 :: A(n, n), b(n), x0(n), x(n)
    data A / (-0.73492d0, 7.11486d0), (0.024839d0, 4.12154d0), (0.274492957d0, 3.7885537d0), (-0.632557864d0, 1.95397735d0) /
    data b / (0.289619736d0, 0.895562183d0), (-0.28475616d0,-0.892163111d0) /
    x0 = (0.0d0, 0.0d0)
    x = bicgm(A, b, x0, n)
    ! do i = 1, n
    !     write(unit=*, fmt=*) x(i)
    ! enddo
    write(unit=*, fmt=*) x
    stop
end program main
