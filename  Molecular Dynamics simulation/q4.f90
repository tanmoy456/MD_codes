program q4
    implicit none
    integer   :: i,N,c
    real*8    :: x,fac,ex,term
    N = 100
    fac = 1.0;ex=1.0
    x =1.0 ; c = 0
    do i =1,N
        fac = fac*i
        term = x**i/fac
        c = c+1
        if (term < .000001) exit
        ex = ex + term
    end do
    print*,"The factorial:",ex
    print*,"The term:",term
    print*,"The number of iteration:",c
end program q4
