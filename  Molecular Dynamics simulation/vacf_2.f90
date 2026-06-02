program vcf
    implicit none 
    integer, parameter  :: n_part = 108
    real*8, parameter   :: dt = 0.001, tf = 200.0 !time step & total time
    integer, parameter  :: time_step = 20000 !2000/(.001*100)
    integer, parameter  :: ndelay = 1000
    real*8,  dimension(time_step,n_part)  :: vx,vy,vz
    real*8,  dimension(0:ndelay)  :: vacf
    real*8    :: avg_vacf
    integer   :: i,j,k,t,t0,tau,m
    real*8    :: s

    open(11,file='velocity.txt')
    open(22,file='vacf.txt')
    open(33,file='avg_vacf.txt')

    do t=1,time_step
        do i=1,n_part
            read(11,*) vx(t,i),vy(t,i),vz(t,i)
        end do
    end do
    m = 0
    do tau=0,ndelay
    	m = 0
    	avg_vacf = 0.0
    	do t0 =100,18000,100
        vacf(tau) = 0.0
        do i=1,n_part
            vacf(tau) = vacf(tau) + vx(t0,i)*vx(t0+tau,i) &  
            +vy(t0,i)*vy(t0+tau,i)+vz(t0,i)*vz(t0+tau,i)
        end do
        vacf(tau) = vacf(tau)/dfloat(n_part)
        write(22,*) tau,vacf(tau)
        avg_vacf = avg_vacf + vacf(tau)
        m = m + 1
        end do
        avg_vacf = avg_vacf/m
        write(33,*) tau,avg_vacf
    end do
    close(11)
    close(22)
    close(33)
end program vcf
