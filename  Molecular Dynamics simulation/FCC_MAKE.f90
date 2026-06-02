!NAME = TANMOY GHOSH
!To make FCC LAttice
program fcc
	implicit none
	integer,parameter  :: n_part = 108 
	real*8,parameter   :: rho = 0.8442
	real*8,parameter   :: l = (n_part/rho)**(1.0/3.0)
	real*8  		   :: x(n_part),y(n_part),z(n_part)
	real*8             :: a
	integer  :: i,j,k,m,cell,ref,counter
	
	open(11,file='fcc_config.txt')
	open(22,file='fcc_initial.dat')
	
	print*,'Enter number of cell in any direction(enter 3 for 108 particle): ' !enter 3 to get 108 particle
	read(*,*) cell
	print*,'The number of cell in any direction: ',cell
	print*,'The number of particle:  ',4*cell**3  !fcc contains 4 particle in a unit cell
	a = l/float(cell)
	print*,'The lattice constant: ', a 
	
	!Sublattice 1
	x(1) = 0.0 ; y(1) = 0.0 ; z(1) = 0.0
	!Sublattice 2
	x(2) = a/2.0 ; y(2) = a/2.0 ; z(2) = 0.0
	!Sublattice 3
	x(3) = 0.0 ; y(3) = a/2.0 ; z(3) = a/2.0
	!Sublattice 4
	x(4) = a/2.0 ; y(4) = 0.0 ; z(4) = a/2.0
	
	m = 0
	counter = 0
	do i=1, cell
		do j=1,cell
			do k=1,cell
				do ref=1,4
					counter = counter + 1
					x(ref+m) = x(ref) + a*float(i-1)
					y(ref+m) = y(ref) + a*float(j-1)
					z(ref+m) = z(ref) + a*float(k-1)	
				end do 
				m = m + 4
			end do
		end do
	end do
	print*,'m:',m
	print*,'counter:',counter
	write(11,*) n_part
	write(11,*) ""
	do i=1,n_part
		write(11,*) 'Ar',x(i),y(i),z(i)
		write(22,*) x(i),y(i),z(i)
	end do 
	write(11,*) ""
	
end program fcc


