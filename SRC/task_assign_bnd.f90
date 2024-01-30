
subroutine task_assign_bnd(f,dimx1,dimx2,dimy1,dimy2,dimz,buff,tag)
	
! this routine assignes the boundary info after MPI exchange

use grid
implicit none
	
integer dimx1, dimx2, dimy1, dimy2, dimz
real f(dimx1:dimx2, dimy1:dimy2, dimz)
	
real buff(*)	! buff for sending data
integer tag

integer i, j, k, n, proc
integer i1, i2, j1, j2

!       The dimensions of the fields in common com3d. Needed by MPI

!       Decode the tag:

! MO,2/26/2016: changed from base 10 to base 5 to reduce tag to below
!  MPI_TAG_UB limit on Cray MPI (e.g., edison, cori @NERSC). It increases 
!  # of allowed variables from 21 to 133 (Marat K., pers. communication)
!  Also changed task_exchange.f90 and task_dispatch.f90
!
! i= tag/100000               ! base 10
! i1 = (tag-i*100000)/10000
! i2 = (tag-i*100000-i1*10000)/1000
! j1 = (tag-i*100000-i1*10000-i2*1000)/100
! j2 = (tag-i*100000-i1*10000-i2*1000-j1*100)/10
! proc =tag-i*100000-i1*10000-i2*1000-j1*100-j2*10

!i= tag/15625                  ! base 5
!i1 = (tag-i*15625)/3125
!i2 = (tag-i*15625-i1*3125)/625
!j1 = (tag-i*15625-i1*3125-i2*625)/125
!j2 = (tag-i*15625-i1*3125-i2*625-j1*125)/25
!proc =tag-i*15625-i1*3125-i2*625-j1*125-j2*25

i= tag/4096                   ! base 4 (WON'T work with UM5 advection)
i1 = (tag-i*4096)/1024
i2 = (tag-i*4096-i1*1024)/256
j1 = (tag-i*4096-i1*1024-i2*256)/64
j2 = (tag-i*4096-i1*1024-i2*256-j1*64)/16
proc =tag-i*4096-i1*1024-i2*256-j1*64-j2*16
	
! From "North":

	  if    (proc.eq.1) then

	     n=0
	     do k=1,dimz
	       do j=nyp1,nyp1+j2
	         do i=1,nx
	           n = n+1
	           f(i,j,k) = buff(n)
	         end do
	       end do
	     end do
	  
! From "North-East":

	  elseif(proc.eq.2) then

	     n=0
	     do k=1,dimz
	       do j=nyp1,nyp1+j2
	         do i=nxp1,nxp1+i2
	           n = n+1
	           f(i,j,k) = buff(n)
	         end do
	       end do
	     end do

! From "East":

	  elseif(proc.eq.3) then
	  
	     n=0
	     do k=1,dimz
	       do j=1,ny
	         do i=nxp1,nxp1+i2
	           n = n+1
	           f(i,j,k) = buff(n)
	         end do
	       end do
	     end do
	  
! From "South-East":

	  elseif(proc.eq.4) then
	  
	     n=0
	     do k=1,dimz
	       do j=-j1,0
	         do i=nxp1,nxp1+i2
	           n = n+1
	           f(i,j,k) = buff(n)
	         end do
	       end do
	     end do
	  
! From "South":

	  elseif(proc.eq.5) then
	  
	     n=0
	     do k=1,dimz
	       do j=-j1,0
	         do i=1,nx
	           n = n+1
	           f(i,j,k) = buff(n)
	         end do
	       end do
	     end do
	     
! From "South-West":

	  elseif(proc.eq.6) then
	  
	     n=0
	     do k=1,dimz
	       do j=-j1,0
	         do i=-i1,0
	           n = n+1
	           f(i,j,k) = buff(n)
	         end do
	       end do
	     end do
	     
! From "West":

	  elseif(proc.eq.7) then
	  
	     n=0
	     do k=1,dimz
	       do j=1,ny
	         do i=-i1,0
	           n = n+1
	           f(i,j,k) = buff(n)
	         end do
	       end do
	     end do
	     	  
! From "North-West":

	  elseif(proc.eq.8) then
	  
	     n=0
	     do k=1,dimz
	       do j=nyp1,nyp1+j2
	         do i=-i1,0
	           n = n+1
	           f(i,j,k) = buff(n)
	         end do
	       end do
	     end do
	     
	  endif
	
	
end
	
