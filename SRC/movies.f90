! module for production of movie (visualization) files written in raw format 
! and postprocessed into gif files.
! Original implementation by Peter Bogenschutz, U. Uta, April 2008
! 

module movies

use grid

implicit none

  ! Variables required for production of movie files
  real cldttmp_xy(nx,ny), cldttmpl_xy(nx,ny), cwppath(nx,ny)
  real sfcflux(nx,ny), sfcpflog10dat(nx,ny)
  real dtmovie, amin, amax 
  integer irecc

CONTAINS

subroutine openmoviefiles()

use vars
use params
implicit none

character *4 rankchar

! Subroutine to open the files needed to create the movies.

write(rankchar,'(i4)') rank

open(unit=80,file='./OUT_MOVIES/'//trim(case)//'_'//trim(caseid)//'_usfc.raw_'//rankchar(5-lenstr(rankchar):4), &
                                              form='unformatted',access='direct',recl=nx*ny)

open(unit=81,file='./OUT_MOVIES/'//trim(case)//'_'//trim(caseid)//'_vsfc.raw_'//rankchar(5-lenstr(rankchar):4), &
                                              form='unformatted',access='direct',recl=nx*ny)

open(unit=82,file='./OUT_MOVIES/'//trim(case)//'_'//trim(caseid)//'_cldtop.raw_'//rankchar(5-lenstr(rankchar):4),& 
                                              form='unformatted',access='direct',recl=nx*ny)

open(unit=83,file='./OUT_MOVIES/'//trim(case)//'_'//trim(caseid)//'_thsfc.raw_'//rankchar(5-lenstr(rankchar):4),&
                                              form='unformatted',access='direct',recl=nx*ny)

open(unit=84,file='./OUT_MOVIES/'//trim(case)//'_'//trim(caseid)//'_qvsfc.raw_'//rankchar(5-lenstr(rankchar):4),&
                                              form='unformatted',access='direct',recl=nx*ny)

open(unit=85,file='./OUT_MOVIES/'//trim(case)//'_'//trim(caseid)//'_sfcprec.raw_'//rankchar(5-lenstr(rankchar):4), &
                                              form='unformatted',access='direct',recl=nx*ny)

open(unit=86,file='./OUT_MOVIES/'//trim(case)//'_'//trim(caseid)//'_cwp.raw_'//rankchar(5-lenstr(rankchar):4), &
                                              form='unformatted',access='direct',recl=nx*ny)

open(unit=87,file='./OUT_MOVIES/'//trim(case)//'_'//trim(caseid)//'_iwp.raw_'//rankchar(5-lenstr(rankchar):4), &
                                              form='unformatted',access='direct',recl=nx*ny)

open(unit=88,file='./OUT_MOVIES/'//trim(case)//'_'//trim(caseid)//'_mse.raw_'//rankchar(5-lenstr(rankchar):4), &
                                              form='unformatted',access='direct',recl=nx*ny)


return

end subroutine openmoviefiles

!----------------------------------------------------------------------

subroutine mvmovies()

use vars
use params
implicit none

call t_startf ('movies')
! Subroutine to call the subroutines required to make the movies.  
! Modify amin and max as you see fit.  

call openmoviefiles()

call cldtoptmp()

!  Surface U-Wind Component
amin=20.
amax=-20.
call plotpix(u(1:nx,1:ny,1),nx,ny,nstep,irecc,amin,amax,80)

!  Surface V-Wind Component
amin=20.
amax=-20.
call plotpix(v(1:nx,1:ny,1),nx,ny,nstep,irecc,amin,amax,81)

!  Cloud Top Temperature
amin=300.
amax=200.
call plotpix(cldttmp_xy(:,:),nx,ny,nstep,irecc,amin,amax,82)


!  Surface MSE
amin=303.
amax=293.
call plotpix(t(1:nx,1:ny,1),nx,ny,nstep,irecc,amin,amax,88)

!  Surface potential temperature
amin=303.
amax=293.
call plotpix(tabs(:,:,1),nx,ny,nstep,irecc,amin,amax,83)

!  Water vapor mixing ratio
amin=0.022
amax=0.013
call plotpix(qv(:,:,1),nx,ny,nstep,irecc,amin,amax,84)

!  Surface precip rate
call getvtrr(nx,ny,nz_gl,qpl(:,:,1),tabs(:,:,1),rho,sfcflux)
amin=alog10(0.01)
amax=alog10(0.00000001)
call convertlog10(sfcflux,nx,ny,amax,sfcpflog10dat)
call plotpix(sfcpflog10dat,nx,ny,nstep,irecc,amin,amax,85)

!  Cloud Water Path
amin=alog10(1000.)
amax=alog10(0.1)
sfcpflog10dat(:,:)=0.
call pathvar(nx,ny,nz_gl,z,qcl(:,:,:),cwppath)
call convertlog10(cwppath,nx,ny,amax,sfcpflog10dat)
call plotpix(sfcpflog10dat,nx,ny,nstep,irecc,amin,amax,86)

!  Ice Water Path
amin=alog10(1000.)
amax=alog10(5.)
sfcpflog10dat(:,:)=0.
cwppath(:,:)=0.
call pathvar(nx,ny,nz_gl,z,qci(:,:,:)+qpi(:,:,:),cwppath)
call convertlog10(cwppath,nx,ny,amax,sfcpflog10dat)
call plotpix(sfcpflog10dat,nx,ny,nstep,irecc,amin,amax,87)
   
irecc = irecc+1
call close_files()

if(masterproc) print*, 'appending OUT_MOVIES/*.raw files. nstep=', nstep

call t_stopf ('movies')
return
  
end subroutine mvmovies

!---------------------------------------------------------------------

subroutine close_files()

! Close the movie files after last timestep of integration

  close(unit=80)
  close(unit=81)
  close(unit=82)
  close(unit=83)
  close(unit=84)
  close(unit=85)
  close(unit=86)
  close(unit=87)
  close(unit=88)

  if(masterproc) then
   open(79,file='./OUT_MOVIES/'//trim(case)//'_'//trim(caseid)//'.info.raw',form='formatted')
   write(79,*) nsubdomains_x,nsubdomains_y
   write(79,*) nx_gl/nsubdomains_x,ny_gl/nsubdomains_y
   write(79,*) irecc
   close(79)
  end if

  return

end subroutine close_files


!---------------------------------------------------------------------

subroutine pixplt_raw(data,nx,ny,amin,amax,iun,irec)

  implicit none

  real amin, amax
  integer nx, ny, i, j, idata, iun, ncolor, irec, count
  parameter(ncolor = 254)
  character*1 iclrs(nx,ny)
  real data(nx,ny), rdata
  
  count=0
  do j=1,ny
    do i=1,nx
      rdata=data(i,j)
      idata=int((rdata-amin)/(amax-amin)*float(ncolor))+1
      idata=max(0,min(idata,ncolor+1))
      iclrs(i,j)=char(idata)
      if (idata .ge. 254) then
        count=count+1
      endif
!      write(*,*) idata, data(i,j), iclrs(i,j)
    enddo
  enddo

  write(iun,rec=irec) iclrs(:,:)

  return

end subroutine pixplt_raw



!----------------------------------------------------------------------

subroutine plotpix(data,nx,ny,n,irec,amin,amax,unitnum)

  implicit none

  integer nx, ny, n, unitnum, irec
  real data(nx,ny), amin, amax

    call pixplt_raw(data,nx,ny,amin,amax,unitnum,irec)

  return

end subroutine plotpix

!---------------------------------------------------------------------

subroutine convertlog10(data,nx,nz,minval,retdat)

  implicit none

  integer i, j, nx, nz
  real minval, data(nx,nz), retdat(nx,nz)

  do j=1,nz
    do i=1,nx
      if (data(i,j) .eq. 0 .or. alog10(data(i,j)) .lt. minval) then
        retdat(i,j)=minval
      else
        retdat(i,j)=alog10(data(i,j))
      endif
    enddo
  enddo

  return

end subroutine convertlog10

!----------------------------------------------------------------------

subroutine getvtrr(nx,ny,nz,qr,tabs,rhopro,flux)

  implicit none
  
  integer i, j, nx, ny, nz
  real qrho, act2, vconr, rho, vtrr, rhofac
  real qr(nx,ny), tabs(nx,ny), p(nx,ny), flux(nx,ny), rhopro(nz)

  real, parameter :: pie=3.141593
  real, parameter :: rnzr=8.e6
  real, parameter :: alin=841.99667
  real, parameter :: gam480=17.8379
  real, parameter :: rhor=1.e3

  act2=pie*rnzr*rhor
  vconr=alin*gam480/(6.*act2**(0.2))

  do j=1,ny
    do i=1,nx
      rho=rhopro(1)
      qrho=rho*qr(i,j)
      rhofac=sqrt(1.226/rho)
      vtrr=amin1(vconr*qrho**(0.2)*rhofac,10.0)
      flux(i,j)=rho*vtrr*qr(i,j)
    enddo
  enddo

  return

end subroutine getvtrr


!----------------------------------------------------------------------

subroutine pathvar(nx,ny,nz,z,var,path)

  implicit none

  integer i, j, k, nx, ny, nz
  real var(nx,ny,nz), path(nx,ny), z(nz+1), rho

  rho=1000.
  path(:,:)=0.
  do k=2,nz
    do j=1,ny
      do i=1,nx
        path(i,j)=path(i,j)+(rho*var(i,j,k))*(z(k)-z(k-1))
      enddo
    enddo
  enddo

  return

end subroutine pathvar


!----------------------------------------------------------------------

subroutine cldtoptmp()

use vars
use params
implicit none

integer chkind, chkind2, tcount, i, j, k
real coef1, downint

cldttmp_xy(:,:)=0.
cldttmpl_xy(:,:)=0.
tcount=0
do j=1,ny
  do i=1,nx
    
    downint=0.
    chkind=1.
    chkind2=1.
    do k=nzm,1,-1  ! Integrate downward
      coef1=rho(k)*dz*adz(k)*dtfactor
      downint=downint+(qcl(i,j,k)+qci(i,j,k))*coef1
      if (downint .ge. 0.1 .and. chkind .eq. 1) then 
        cldttmp_xy(i,j)=tabs(i,j,k)
        chkind=2.
      endif      
    end do

    do k=1,nzm-1
      if (qcl(i,j,k)+qci(i,j,k) .gt. 0 .and. qcl(i,j,k+1)+qci(i,j,k+1) .eq. 0 .and. chkind2 .eq. 1) then
        cldttmpl_xy(i,j)=tabs(i,j,k)
!        write(*,*) z
        chkind2=2.
      endif
    enddo

    if (chkind2 .eq. 1) then 
      cldttmpl_xy(i,j)=299.88
    endif

    if (chkind .eq. 1) then 
      cldttmp_xy(i,j)=299.88
      tcount=tcount+1
    endif
   
  end do
end do

return 

end subroutine cldtoptmp

end module movies
