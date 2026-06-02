!NAME = TANMOY GHOSH
!Monte Carlo PROGRAM in 3d 
module parameters
    implicit none
    integer, parameter  :: n_part = 108 !number of particles
    real*8, parameter   :: rho = 0.8442 !density of the system
    real*8, parameter   :: sigma = 1.0  !particle diameter
    real*8, parameter   :: Temp = 0.728 !temperature of the system
    integer, parameter  :: save_time = 100 !after 1000 steps save data
    real*8, parameter   :: l = (n_part/rho)**(1.0/3.0) !system size(for cubic lx=ly=lz=l)
    integer, parameter  :: a = (l*l*l/n_part)**(1.0/3.0) !lattice constant
    integer, parameter  :: n = l/a !number of particles in any one direction
    integer, parameter  :: ncycle = 100000 !Total Monte Carlo steps
    real*8, parameter   :: del_max = 0.2  !delx 
	real*8, parameter   :: rc = 2.5, eps = 1.0 !cutoff and epsilon
    real*8, parameter   :: uc = 4.0*eps*(((sigma/rc)**12)-((sigma/rc)**6)) !cutoff potential
    real*8, parameter   :: fc = 4.0*eps*((12.0*sigma**12/rc**13)-(6.0*sigma**6/rc**7)) !cutoff force

    real*8, dimension(n_part) :: x,y,z    !array for positions at(t)
    
    integer            :: i,j,k,icycle,o
    real*8             :: ran_num,counter
    real*8             :: dx,dy,dz,t,xo,yo,zo,xn,yn,zn,xp,yp,zp
    real*8             :: r2i,r6i,r_sq
    real*8             :: PE,eno,enn,u,dE
    real*8             :: strt_time,end_time
end module parameters
program mc
    use parameters
    call cpu_time(strt_time)
    open(11,file='fcc_initial.dat') !fcc
    !open(11,file='sc_initial.dat') !sc
    open(22,file='energy_MC.dat')
    open(33,file='final_pos.xyz')
    open(44,file='enn.txt')
    
    call initialisation 
    call LJ_energy
    print*,' Initial Potential Energy: ',PE/n_part
    counter = 0.0

    do icycle =1,ncycle
        call mcmove
       if (mod(icycle,save_time) == 0) then
      ! if (icycle>1000) then
            write(22,*)icycle,PE/n_part
            call save_pos
       end if
    end do

    print*,' Final Potential Energy: ',PE/n_part
    print*,"The numbe rof accepted move",counter
    print*,' The acceptance probability: ',(counter*100)/float(ncycle)
    close(11)
    close(22)
    close(33)
    call cpu_time(end_time)
    print*, "Program run time : ", (end_time-strt_time)
end program mc
!++++++++++++++++++++
!INITIALIZE POSITION 
!++++++++++++++++++++
subroutine initialisation
    use parameters
    !POSITION
    do i=1,n_part
        read(11,*)x(i),y(i),z(i)
    end do
end subroutine initialisation
!+++++++++++++++++++++++++++++++++++++++
!L_J ENERGY CALCULATION(pairwise)
!+++++++++++++++++++++++++++++++++++++++
subroutine LJ_energy
    use parameters
    PE = 0.0
    do i = 1, n_part-1
        do j = i+1, n_part
            dx = x(i) - x(j)
            dy = y(i) - y(j)
            dz = z(i) - z(j)
            
            dx = dx - l*nint(dx/l)
            dy = dy - l*nint(dy/l)
            dz = dz - l*nint(dz/l)
            
            r_sq = dx*dx+dy*dy+dz*dz
            if(r_sq <= rc*rc) then 
                
                r2i = (1.0/r_sq)  !1/r^2
                r6i = (r2i)**3    !1/r^6

                PE = PE + 4.0*r6i*(r6i-1.0) + fc*sqrt(r_sq) - uc -fc*rc  !potential energy 
            end if
        end do
    end do
end subroutine LJ_energy
!+++++++++++++++++++++++++++++++++++++++
!ENERGY with o'th particle with the other
!+++++++++++++++++++++++++++++++++++++++
subroutine energy
    use parameters
    u = 0.0
    do i = 1, n_part
        if( i .ne. o) then
            dx = x(i) - xp
            dy = y(i) - yp
            dz = z(i) - zp
            
            dx = dx - l*nint(dx/l)
            dy = dy - l*nint(dy/l)
            dz = dz - l*nint(dz/l)

            r_sq = dx*dx+dy*dy+dz*dz
            if(r_sq <= rc*rc) then 
                
                r2i = (1.0/r_sq)  !1/r^2
                r6i = (r2i)**3    !1/r^6

                u = u + 4.0*r6i*(r6i-1.0) + fc*sqrt(r_sq) - uc -fc*rc  !potential energy 
            end if
        end if
    end do
end subroutine energy
!+++++++++++
!MCMOVE
!+++++++++++
subroutine mcmove
    use parameters
    call random_number(ran_num)
    o = int(ran_num*n_part) + 1
    xo = x(o) ; yo = y(o) ; zo = z(o)
    xp = xo ; yp = yo ; zp = zo
    call energy
    eno = u
    call random_number(ran_num) ; xn = xo + (ran_num - 0.5d0)*del_max
    call random_number(ran_num) ; yn = yo + (ran_num - 0.5d0)*del_max
    call random_number(ran_num) ; zn = zo + (ran_num - 0.5d0)*del_max
    
    if (xn > l) xn=xn-l
    if (yn > l) yn=yn-l
    if (zn > l) zn=zn-l
    
    if (xn < 0.0) xn=xn+l
    if (yn < 0.0) yn=yn+l
    if (zn < 0.0) zn=zn+l

    xp = xn ; yp = yn ; zp = zn

    call energy
    enn = u
    dE = enn - eno
    write(44,*) dE
    if (dE < 0.0) then
        PE = PE + dE
        x(o) = xn ; y(o) = yn ; z(o) = zn
        counter = counter + 1.0
    else 
        call random_number(ran_num)
        if (ran_num <= exp(-dE/(Temp))) then
             PE = PE + dE
             x(o) = xn ; y(o) = yn ; z(o) = zn
             counter = counter + 1.0
        end if
    end if
end subroutine mcmove
!++++++++++++++
!save position
!++++++++++++++
subroutine save_pos
    use parameters
  !  write(33,*) n_part
   ! write(33,*) ""
    do i=1,n_part
       ! write(33,*) 'Ar', x(i),y(i),z(i)  !writing final position 
        write(33,*) x(i),y(i),z(i)
    end do
    !write(33,*) ""
end subroutine save_pos
