!Inititalise the system randomly
program init_random 
    implicit none 
    integer     :: i,j,a,b
    real        :: r1,r2,r3
   ! real        :: x,y,z
    real*8,parameter   :: rho = 0.8442
    real*8,parameter   :: sigma = 1.0
    integer,parameter  :: n_part = 216
    real*8,parameter   :: l = (n_part/rho)**(1.0/3.0)
    
    real*8,dimension(n_part)   :: x,y,z
    real*8             :: dx,dy,dz,r
    
    open(11,file = 'init_random_jmol.xyz') !saving in xyz format to view in jmol
    open(22,file = 'init_random.txt') !saving in txt format to view in gnuplot
    a = 0
    write(11,*) n_part
    write(11,*) ""
    
    	call random_number(r1) 
        call random_number(r2) 
        call random_number(r3)
         
        x(1) = r1*(l-a) + a
        y(1) = r2*(l-a) + a
        z(1) = r3*(l-a) + a
        write(22,*) x(1),y(1),z(1)
        write(11,*) "Ar",x(1),y(1),z(1)    
    do i=2,n_part
    
44      call random_number(r1) 
        call random_number(r2) 
        call random_number(r3)
         
        x(i) = r1*(l-a) + a
        y(i) = r2*(l-a) + a
        z(i) = r3*(l-a) + a
       
        do j=1,i-1
        	dx = x(i)-x(j)
        	dy = y(i)-y(j)
        	dz = z(i)-z(j)
        	
        	r = sqrt(dx*dx+dy*dy+dz*dz)
        	if (r .le. sigma) goto 44 
        end do	
        write(11,*) "Ar",x(i),y(i),z(i)        
        write(22,*) x(i),y(i),z(i)
    end do
    write(11,*) ""

    close(11)
    close(22)
end program init_random
