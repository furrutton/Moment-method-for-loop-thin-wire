program main
    integer*4 :: a(2, 2), i
    a(1,1) = 11
    a(1,2) = 12
    a(2,1) = 21
    a(2,2) = 22
    do i = 1, 2
        write(unit=*, fmt=*) a(i,:) 
    enddo
    data a / 11, 12, 21, 22 /
    write(unit=*, fmt=*) "data ways: "
    
    do i = 1, 2
        write(unit=*, fmt=*) a(i,:) 
    enddo
    stop
end program main
