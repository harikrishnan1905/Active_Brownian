use omp_lib
implicit none
integer :: m,t,tmax,i,iseed,n,tequl
real (kind = 8) :: gasdev,ran2,vp
real (kind = 8), dimension (:) , allocatable :: xx,xy,ax,ay,the
real (kind =8) :: L,dt,temp,rho,d,pe
external :: gasdev,ran2
n = 64							
dt = 0.00001D0  
!gam = 10000.0d0						! friction coefficient
m = n*n 						! total number of particles
Pe = 150.00d0
allocate(xx(m),xy(m),the(m))				! 'the(i)' specifies the direction of the i-th particle	
allocate(ax(m),ay(m))
ax = 0.D0 ; ay = 0.D0 
xx = 0.D0 ; xy = 0.D0 ; iseed = -887  
the = 0.D0
temp = 1.0d0 
vp = Pe  						! self propulsion velocity
d = 1.0!temp/gam	      				! translational diffusion coefficient
open (400,file='finvel150_second.txt') 
open (100,file='density150_second.txt') 
tmax = 2000000
tequl = 1000000
rho = 0.80d0						! density
l = sqrt(m/rho)						! simulation box is square with side l
call initial(xx,xy,the,m,n,l,dt,iseed)			! initialzing the system
do t = 1,tmax						! time evolution
	call force(ax,ay,xx,xy,m,L)			! force calculation in each time step	
	call integrate(ax,ay,xx,xy,the,m,dt,L,iseed,d,vp)	! updating positions and angles at each time step
	if ((t.ge.tequl).and.(mod(t,10000).eq. 0))then
		do i=1,m
			write(100,*) xx(i),xy(i)
		end do	
	end if
	if (t.eq.tmax) then
		do i=1,m
			write(400,*) xx(i),xy(i)
		end do
	end if
end do
end


subroutine initial(xx,xy,the,m,n,l,dt,iseed)	
implicit none
integer ,intent(in) :: iseed,m,n
integer :: i,j,k,p,g
real (kind = 8),intent(inout) :: xx(m),xy(m),the(m)
real (kind = 8)::gasdev,a,ran2
real (kind = 8), intent (in) :: dt,l
REAL (kind = 8), PARAMETER :: Pi = 3.14159265
external :: gasdev,ran2					
! inputs are the density and the toatal number of particles
a = l/n							! a specifies the base of the triangle
g = 1							! triangular lattice but not equilateral
do i =0,m-1						
	if (g.le.m) then
		xx(g) = a*i				
		xy(g) = 0.0d0
		p = int((i)/n)
		k = mod(p,2)
		if (i.gt.(n-1)) then			
			xx(g) = xx(g-n) + a*0.5D0	
			if (k.eq. 0) then		
				xx(g) = xx(g) - a*1.0d0
			end if
			xy(g) = xy(g-n) + a		
		end if
	end if	
	g = g + 1
end do
! assigning random orientations to each particle
do i =1,m
	the(i) = 2.0d0*Pi*ran2(iseed)			! the angle is a radom number between 0 and 2 pi
	!write(10,*) xx(i),xy(i)
end do
end subroutine initial


subroutine force(ax,ay,xx,xy,m,L)
implicit none
integer , intent (in) :: m
integer :: i,j
real (kind = 8), intent (inout) :: ax(m),ay(m)
real (kind = 8), intent (in) :: xx(m),xy(m),L
real (kind = 8) :: rsqd,ff,r2i,r6i,rcut,rcut2,virij,rx,ry
rx = 0.D0 ; ry = 0.D0 
ax = 0.D0 ; ay = 0.D0 
rcut = 1.12246205 !2**(1.0d0/6.0d0)			! cut off distance for the WCA potential is 2^(1/6)
rcut2 = rcut*rcut !1.259921054
call omp_set_num_threads(10)
!$omp parallel do private(rx,ry,rsqd,r2i,r6i,virij,ff,j) schedule(dynamic) 
! calculation of pair-wise distances and forces on each particle
do i=1,m-1
	do j=i+1,m
		rsqd = 0.D0
		rx = xx(i) - xx(j)
		rx = rx-L*nint(rx/L)	! distance btw i and nearest periodic image of j
		ry = xy(i) - xy(j)
		ry = ry-L*nint(ry/L)						
		rsqd = rsqd + rx*rx +ry*ry 
	if (rsqd.le.rcut2) then;
		r2i = 1/rsqd
		r6i = r2i**3
		virij = 48*(r6i*r6i-0.5D0*r6i) 		
                ff = virij*r2i	
                if (ff .gt. 10000) then	
			print*, 'HAHAHAHAHAHA  You are stuck mate!!!!'
			EXIT
		end if
		! force on the particle
		!$omp atomic
		ax(i) = ax(i) + rx*ff		
		!$omp atomic				
		ay(i) = ay(i) + ry*ff	
		!$omp atomic									
		ax(j) = ax(j) - rx*ff
		!$omp atomic
		ay(j) = ay(j) - ry*ff
!		end if						
	end if			
	end do	
end do
!$omp end parallel do
end subroutine force



subroutine integrate(ax,ay,xx,xy,the,m,dt,L,iseed,d,vp)
implicit none
integer , intent (in) :: m,iseed
integer :: i,j
real (kind = 8), intent (in) :: ax(m),ay(m),dt,L,d,vp  ! d and dr translational and rotational diffusion coefficients
real (kind = 8), intent (inout) :: xx(m),xy(m),the(m)
real (kind = 8) :: g1,g2,g3,gasdev,dr
external :: gasdev
dr = 3.0d0*d						! in the low Reynolds number regime
do i=1,m	
	g1 = gasdev(iseed)				! Gaussian random number with 0 mean and unit variance
	xx(i) = xx(i) + dt*(vp*cos(the(i))+ ax(i))*d + sqrt(2*d*dt)*g1	
        g2 = gasdev(iseed)				
	xy(i) = xy(i) + dt*(vp*sin(the(i))+ ay(i))*d + sqrt(2*d*dt)*g2        
        g3 = gasdev(iseed)
	the(i) = the(i) + sqrt(2*dr*dt)*g3			! orientation of the particle is modified due to the the fluctuations only
	xx(i) = xx(i) - L*anint((xx(i)/L)-(0.5D0))				
	xy(i) = xy(i) - L*anint((xy(i)/L)-(0.5D0))	! if the particle is outside the simulation box, bring it back	
end do
end subroutine integrate
