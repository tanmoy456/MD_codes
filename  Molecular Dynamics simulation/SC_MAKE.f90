!Inititalise the system in simple cubic format
program init_cubic
    implicit none
    integer,parameter  :: n_part = 108
    real*8,parameter   :: rho = 0.8442
    real*8,parameter   :: l = (n_part/rho)**(1.0/3.0) 
    real*8,parameter   :: a = (l*l*l/n_part)**(1.0/3.0)
    integer,parameter  :: n = nint(l/a)

    integer  :: i,j,k,c
    real*8   :: ran
    open(11,file='sc_jmol.xyz') !format for view in jmol
    open(22,file='sc_initial.dat') !format for view in gnuplot
    !open(33,file='init_avogadro.xyz')
    print*,"l:",l
    print*,"Lattice constant:",a
    print*,"n:",n
    write(11,*) n_part
    write(11,*) ""
    c = 1
    do i=1,n
        do j=1,n
            do k=1,n
                !write(11,*)"n_part",i*a,j*a,k*a
                write(11,*) "Ar",(i-1)*a,(j-1)*a,(k-1)*a
                write(22,*) (i-1)*a,(j-1)*a,(k-1)*a
                c = c+1
                if(c>n_part) exit
            end do
            if(c>n_part) exit
        end do
        if(c>n_part) exit 
    end do
    write(11,*) ""
    close(11)
    close(22)
end program init_cubic
    
