module tools
    contains
    function func1(x)
        complex*16 :: func1
        real*8 :: x
        complex*16 :: j = (0.0d0, 1.0d0)
        func1 = exp(-2*j*x)
    end function func1


    function gauss(func, a, b)
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
        ! gauss = sum(dot_product(func(t), w))
        gauss = gauss*(b-a)/2.0d0
    end function gauss

end module tools

program main
    use tools
    complex*16 :: s
    real*8 :: a = 0.0d0
    complex*16 :: j = (0.0d0, 1.0d0)
    write(unit=*, fmt=*)  gauss(func1, 1.0d0, 2.0d0)
    s = (sin(4.0d0)+j*cos(4.0d0))/2.0d0 - (sin(2.0d0)+j*cos(2.0d0))/2.0d0
    write(unit=*, fmt=*) s
    
    stop
end program main

