!NAME = TANMOY GHOSH
!MD PROGRAM in 3d 
module parameters
    implicit none
    integer, parameter  :: n_part = 108 !number of particles
    real*8, parameter   :: rho = 0.5 !density of the system
    real*8, parameter   :: sigma = 1.0  !particle diameter
    real*8, parameter   :: Temp = 1.228 !temperature of the system
    real*8, parameter   :: dt = 0.001, tf = 200.0 !time step & total time
    real*8, parameter   :: save_time = 100.0 
    real*8, parameter   :: l = (n_part/rho)**(1.0/3.0) !system size(for cubic lx=ly=lz=l)
    integer, parameter  :: a = (l*l*l/n_part)**(1.0/3.0) !lattice constant
    integer, parameter  :: n = l/a !number of particles in any one direction

	real*8, parameter  :: rc = 2.5, eps = 1.0 !cutoff and epsilon
    real*8, parameter  :: uc = 4.0*eps*(((sigma/rc)**12)-((sigma/rc)**6)) !cutoff potential
    real*8, parameter  :: fc = 4.0*eps*((12.0*sigma**12/rc**13)-(6.0*sigma**6/rc**7)) !cutoff force

    real*8, dimension(n_part)   :: x,y,z    !array for positions at(t)
    real*8, dimension(n_part)   :: vx,vy,vz !array for velocity
    real*8, dimension(n_part)   :: fx,fy,fz !array for force
    real*8, dimension(n_part)   :: xp,yp,zp !positions at(t-dt)
    real*8, dimension(n_part)   :: xf,yf,zf !positions at(t+dt)
    
    integer            :: i,j,k
    real*8             :: ran1,ran2,ran3
    real*8             :: dx,dy,dz,t
    real*8             :: r2i,r6i,r_sq,fx_ij,fy_ij,fz_ij
    real*8             :: sum_vx,sum_vy,sum_vz,sum2_vx,sum2_vy,sum2_vz
    real*8             :: fs,PE,KE 
    real*8             :: strt_time,end_time
end module parameters
program md
    use parameters
    call cpu_time(strt_time)  !for calculation of program runtime
    open(11,file='init_cubic.txt')
    open(33,file='energy.txt')
    open(44,file='final_pos.txt')
    open(55,file='velocity.txt')
    
    call initialisation  
    call force_cal !calling force to get the initial potential energy
    !Printing Initial TEMP.,PE,KE.
    print*,"Initial Temperature:",(2.0d0*KE)/(3.0d0*n_part)
    !print*,"Initial Potential Energy:", PE/float(n_part)
    !print*,"Initial Kinetic Energy:", KE/float(n_part)
    t = 0.0
    !Time loop
    do while(t <= tf)
        call force_cal
        call verlet
        call PBC
        !write(33,*)t,KE/float(n_part),PE/float(n_part),(PE+KE)/float(n_part) 
        if (t>(tf/2.0)) then 
        	!print*,t
            if(mod(t/dt,save_time)==0.0) then
            !print*,t/dt
                call visualise
            end if
        end if
        !if(mod(t/dt,save_time)==0.0) then
        	 call velocity_save
        !end if
        t = t + dt
    end do
    !Final Position
    
    print*,"Final Temperature:",(2.0d0*KE)/(3.0d0*n_part)
    !print*,"Final Potential Energy:", PE/float(n_part)
   ! print*,"Final Kinetic Energy:", KE/float(n_part)

    close(11)
    !close(22)
    close(33)
    close(44)
    close(55)
    call cpu_time(end_time)
    print*, "Program run time : ", (end_time-strt_time)
end program md
!+++++++++++++++++++++++++++++++++
!INITIALIZE POSITION AND MOMENTUM
!+++++++++++++++++++++++++++++++++
subroutine initialisation
    use parameters
    !POSITION
    do i=1,n_part
        read(11,*)x(i),y(i),z(i)  !initial position from init_cubic.txt
    end do
    sum_vx = 0.0 ; sum2_vx = 0.0
    sum_vy = 0.0 ; sum2_vy = 0.0
    sum_vz = 0.0 ; sum2_vz = 0.0
    KE = 0.0
    !VELOCITY
    do i=1,n_part
        call random_number(ran1); vx(i) = ran1 - 0.5
        call random_number(ran2); vy(i) = ran2 - 0.5
        call random_number(ran3); vz(i) = ran3 - 0.5
        
        sum_vx = sum_vx + vx(i) ; sum2_vx = sum2_vx + vx(i)**2
        sum_vy = sum_vy + vy(i) ; sum2_vy = sum2_vy + vy(i)**2
        sum_vz = sum_vz + vz(i) ; sum2_vz = sum2_vz + vz(i)**2
    end do
    sum_vx = sum_vx/n_part ; sum2_vx = sum2_vx/n_part 
    sum_vy = sum_vy/n_part ; sum2_vy = sum2_vy/n_part
    sum_vz = sum_vz/n_part ; sum2_vz = sum2_vz/n_part
    
    fs = sqrt((3.0*Temp)/(sum2_vx+sum2_vy+sum2_vz))
    print*,'scaling factor(fs):',fs
    !C.O.M CORRECTION
    do i=1,n_part
        vx(i) = (vx(i) - sum_vx)*fs !rescaled velocities
        vy(i) = (vy(i) - sum_vy)*fs 
        vz(i) = (vz(i) - sum_vz)*fs 
        !this calculation of K.E needed to calculate initial temperature 
        KE = KE + 0.5*(vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i))
        xp(i) = x(i) - vx(i)*dt  !positions at (t-dt)
        yp(i) = y(i) - vy(i)*dt
        zp(i) = z(i) - vz(i)*dt
    end do
end subroutine initialisation
!+++++++++++++++++++++++++++++++++++++++
!FORCE AND POTENTIAL ENERGY CALCULATION
!+++++++++++++++++++++++++++++++++++++++
subroutine force_cal
    use parameters
    PE = 0.0
    do i=1,n_part
        fx(i) = 0.0
        fy(i) = 0.0
        fz(i) = 0.0
    end do
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

                fx_ij= (48.0*r2i*r6i*(r6i-0.5) - fc)*dx*sqrt(r2i)                                        
                fy_ij= (48.0*r2i*r6i*(r6i-0.5) - fc)*dy*sqrt(r2i)
                fz_ij= (48.0*r2i*r6i*(r6i-0.5) - fc)*dz*sqrt(r2i) 

                fx(i) = fx(i) + fx_ij
                fy(i) = fy(i) + fy_ij 
                fz(i) = fz(i) + fz_ij 
            
                fx(j) = fx(j) - fx_ij   !as F_{ij}=-F_{ji}
                fy(j) = fy(j) - fy_ij
                fz(j) = fz(j) - fz_ij 
    
                PE = PE + 4.0*r6i*(r6i-1.0) + fc*sqrt(r_sq) - uc -fc*rc  !potential energy 
            end if
        end do
       ! write(22,*) i,sqrt(fx(i)**2+fy(i)**2+fz(i)**2)
    end do
end subroutine force_cal
!+++++++++++++++++++
!VERLET INTEGRATION
!+++++++++++++++++++
subroutine verlet
    use parameters
    KE = 0.0
    do i=1,n_part
        xf(i) = 2.0*x(i) - xp(i) + fx(i)*dt**2
        yf(i) = 2.0*y(i) - yp(i) + fy(i)*dt**2
        zf(i) = 2.0*z(i) - zp(i) + fz(i)*dt**2
        
        vx(i) = (xf(i) - xp(i))/(2.0*dt)
        vy(i) = (yf(i) - yp(i))/(2.0*dt)
        vz(i) = (zf(i) - zp(i))/(2.0*dt)
        
        KE = KE + 0.5*(vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i)) !Kinetic energy
     
        xp(i) = x(i) ; yp(i) = y(i) ;  zp(i) = z(i)
        x(i) = xf(i) ; y(i) = yf(i) ;  z(i) = zf(i) 
    end do
end subroutine verlet
!++++++++++++++++++++++++++++++++
!periodic boundary condition(PBC)
!++++++++++++++++++++++++++++++++
subroutine PBC
    use parameters
    do i=1,n_part
        if(x(i) > l ) then 
            x(i) = x(i) - l
            xp(i) = xp(i) - l
        else if(x(i) < 0.0) then 
            x(i) = x(i) + l
            xp(i) = xp(i) + l
        end if
        if(y(i) > l ) then 
            y(i) = y(i) - l
            yp(i) = yp(i) - l
        else if(y(i) < 0.0) then 
            y(i) = y(i) + l
            yp(i) = yp(i) + l
        end if
        if(z(i) > l ) then 
            z(i) = z(i) - l
            zp(i) = zp(i) - l
        else if(z(i) < 0.0) then 
            z(i) = z(i) + l
            zp(i) = zp(i) + l
        end if 
    end do
end subroutine PBC
!+++++++++++
!Visualise
!+++++++++++
subroutine visualise
    use parameters
    !write(44,*) n_part
    !write(44,*) ""
    do i=1,n_part
        write(44,*) x(i),y(i),z(i)  !writing final position 
    end do
   ! write(44,*) ""
end subroutine visualise
!+++++++++++++++++++
!Save Velocity for each time step
!++++++++++++++++++++
subroutine velocity_save
    use parameters
    do i=1,n_part
        write(55,*) vx(i),vy(i),vz(i)
    end do
end subroutine velocity_save
