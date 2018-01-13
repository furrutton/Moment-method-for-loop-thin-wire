module cgm
    contains
    function bicgm(A, b, x0, N)
        integer*4 :: N, i, j, k, maxit
        complex*16 :: A(N, N)
        complex*16 :: b(N), x0(N), x(N)
        complex*16 :: r0(N), r(N), d(N), r1(N)
        complex*16 :: bicgm(N)
        complex*16 :: beta(N)

        r0 = b-matmul(A, x0)
        ! d0 = r0
        d = r0
        r = r0
        k = 0
        maxit = 2000

        open(11, file='text.txt')
        do i = 1, maxit
            k = k+1
            alpha = dot_product(matmul(transpose(conjg(A)),r),matmul(transpose(conjg(A)),r))/&
                &dot_product(matmul(transpose(conjg(A)),d),matmul(transpose(conjg(A)),d))
            x = x + alpha*d
            if ( dble(dot_product(matmul(A, x), matmul(A, x))) .lt.  1.0d-6) exit
            ! write(unit=11, fmt=*) dble(dot_product(matmul(A, x), matmul(A, x)))
            r1 = r
            r = r - alpha*matmul(A, d)
            ! write(unit=11, fmt=*) r
            beta = dot_product(matmul(transpose(conjg(A)), r), matmul(transpose(conjg(A)), r))/&
                &dot_product(matmul(transpose(conjg(A)), r1), matmul(transpose(conjg(A)), r1))
            d = matmul(transpose(conjg(A)), r) + beta*d
        enddo
        close(11)
        if ( k .eq. maxit )then
            write(unit=*, fmt=*) k,"not Converged !"
        else
            write(unit=*, fmt=*) k
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
    data x0 / (0.0d0, 0.0d0), (0.0d0, 0.0d0) /
    x = bicgm(A, b, x0, n)
    write(unit=*, fmt=*) x
    stop
end program main
