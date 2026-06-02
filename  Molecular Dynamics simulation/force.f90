!program to calculete force on each particle
program force_cal
    implicit none
    integer, parameter  :: n_part = 216 !number of particles
    real*8, parameter   :: rho = 0.8442 !density of the system
    real*8, parameter   :: sigma = 1.0  
    real*8, parameter   :: l = (n_part/rho)**(1.0/3.0) !system size(for cubic lx=ly=lz=l)
    integer, parameter  :: a = (l*l*l/n_part)**(1.0/3.0) !lattice constant
    integer, parameter  :: n = l/a !number of particles in any one direction

	real*8, parameter  :: rc = 1.5, eps = 1.0 !cutoff and epsilon
    real*8, parameter  :: uc = 4.0*eps*(((sigma/rc)**12)-((sigma/rc)**6)) !cutoff potential
    real*8, parameter  :: fc = 4.0*eps*((12.0*sigma**12/rc**13)-(6.0*sigma**6/rc**7)) !cutoff force

    real*8, dimension(n_part)   :: x,y,z !array for positions
    real*8, dimension(n_part)   :: px,py,pz !array for momentum
    real*8, dimension(n_part)   :: fx,fy,fz !array for force
    
    integer            :: i,j,k
    !real*8             :: ran1,ran2,ran3
    real*8             :: dx,dy,dz,dr
    real*8             :: r2i,r6i,r_sq,fx_ij,fy_ij,fz_ij
    real*8             :: PE !potential energy
    open(11,file='init_cubic.txt')
    do i=1,n_part
        read(11,*)x(i),y(i),z(i)  !reading initial position from init_cubic.txt
    end do 
    !print*,x
    !+++++++++++++++++++++++++++++++++++++++
    !FORCE AND POTENTIAL ENERGY CALCULATION
    !+++++++++++++++++++++++++++++++++++++++
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
            
            if(abs(dx) .ge. l)  dx = (l-abs(dx))*((-1.0d0*dx)/abs(dx))
            if(abs(dy) .ge. l)  dy = (l-abs(dy))*((-1.0d0*dy)/abs(dy))
            if(abs(dz) .ge. l)  dz = (l-abs(dz))*((-1.0d0*dz)/abs(dz))

            r_sq = dx*dx+dy*dy+dz*dz
            if(r_sq <= rc*rc) then 
                
                r2i = (1.0/r_sq)
                r6i = (r2i)**3

                fx_ij = (48.0*eps*r2i*r6i*(r6i-0.5) - fc)*dx*sqrt(r2i)
                fy_ij = (48.0*eps*r2i*r6i*(r6i-0.5) - fc)*dy*sqrt(r2i)
                fz_ij = (48.0*eps*r2i*r6i*(r6i-0.5) - fc)*dz*sqrt(r2i)

                fx(i) = fx(i) + fx_ij
                fy(i) = fy(i) + fy_ij
                fz(i) = fz(i) + fz_ij

                fx(j) = fx(j) - fx_ij   !as F{ij} = - F{ji}
                fy(j) = fy(j) - fy_ij
                fz(j) = fz(j) - fz_ij

                PE = PE + 4.0*eps*r6i*(r6i-1.0) + fc*sqrt(r_sq) - uc -fc*rc !potential energy
            end if
        end do
    end do

end program force_cal
