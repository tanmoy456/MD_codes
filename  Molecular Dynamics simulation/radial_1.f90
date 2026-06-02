program radial
    implicit none
    integer, parameter :: n_part = 108 !number of particle
    real*8,  parameter :: rho = 0.5!density of system
    real*8, parameter   :: dt = 0.001, tf = 200.0 !time step & total time
    !real*8, parameter   :: save_time = 100.0
    integer, parameter  :: time_step = 1000!dt/(tf*save_time*2)
   !integer, parameter :: save_step = 1000 !after every 100 step save the result
    real*8,  parameter :: pi = 4.0d0*atan(1.0d0) !value of pi
    real*8,  parameter :: l = (n_part/rho)**(1.0d0/3.0d0)!size of the system
    real*8,dimension(time_step,n_part)   :: x,y,z
    integer  :: i,j,k,t
    real*8   :: dx,dy,dz,r,vol,dr
    integer,parameter  :: n_bin = 800
    real*8, dimension(n_bin) :: H
    integer     :: bin_index,ngr

    open(11,file='final_pos.txt')
    open(22,file='rdf.txt')

    do t=1,time_step
        do i=1,n_part
            read(11,*) x(t,i),y(t,i),z(t,i)
        end do
    end do
    dr = l/(2*n_bin)
    !dr = 0.001
    ngr = 0
    H = 0.0d0
    do t=1,time_step
        !if(t>20) then
       ! if(mod(t,save_step)==0) then
            ngr = ngr + 1
           ! print*,t
            do i=1,n_part-1
                do j=i+1,n_part
                    dx = x(t,i) - x(t,j)
                    dy = y(t,i) - y(t,j)
                    dz = z(t,i) - z(t,j)
                    
                    dx = dx - nint(dx/l)*l
                    dy = dy - nint(dy/l)*l
                    dz = dz - nint(dz/l)*l
                
                    r = sqrt(dx**2+dy**2+dz**2)
                    if (r .lt. l/2.0) then 
                        bin_index = int(r/dr)
                        H(bin_index) = H(bin_index) + 2  
                    end if
                end do
            end do
       ! end if
    end do
    do i = 1, n_bin-1
        vol = (4.0d0/3.0d0)*pi*((i+1)**3-i**3)*dr**3*rho
        write(22,*) dr*(i+0.5d0), H(i)/(ngr*n_part*vol)
    end do
   
 close(11)
 close(22)

end program radial 
