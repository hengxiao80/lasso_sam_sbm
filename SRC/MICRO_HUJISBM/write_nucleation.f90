     
subroutine write_nucleation
	
use vars
use rad, only: qrad
use microphysics
implicit none
character *80 filename,long_name
character *8 name
character *10 timechar
character *4 rankchar
character *6 filetype
character *10 units
character *10 c_z(nzm),c_p(nzm),c_dx, c_dy, c_time
integer i,j,k,nfields,nfields1
real tmp(nx,ny,nzm)

nfields=2 ! number of 3D fields to save
nfields1=0

if(masterproc) then

  write(rankchar,'(i4)') nsubdomains
  write(timechar,'(i10)') nstep
  do k=1,11-lenstr(timechar)-1
    timechar(k:k)='0'
  end do

  if(RUN3D) then
      filetype = '.com3D'
  else
      filetype = '.com2D'
  endif
    filename='./OUT_3D/'//trim(case)//'_'//trim(caseid)//'_nucl_'// &
        rankchar(5-lenstr(rankchar):4)//filetype
      if(nrestart.eq.0 .and.notopened3D) then
         open(46,file=filename,status='unknown',form='unformatted')	
      else
         open(46,file=filename,status='unknown', &
                              form='unformatted', position='append')
      end if
      notopened3D=.false.

  if(save3Dbin) then

    write(46) nx,ny,nzm,nsubdomains,nsubdomains_x,nsubdomains_y,nfields
    do k=1,nzm
      write(46) z(k) 
    end do
    do k=1,nzm
      write(46) pres(k)
    end do
    write(46) dx
    write(46) dy
    write(46) nstep*dt/(3600.*24.)+day0

  else

    write(long_name,'(8i4)') nx,ny,nzm,nsubdomains, &
                                   nsubdomains_x,nsubdomains_y,nfields
    do k=1,nzm
       write(c_z(k),'(f10.3)') z(k)
    end do
    do k=1,nzm
       write(c_p(k),'(f10.3)') pres(k)
    end do
    write(c_dx,'(f10.0)') dx
    write(c_dy,'(f10.0)') dy
    write(c_time,'(f10.5)') nstep*dt/(3600.*24.)+day0
	
    write(46) long_name(1:32)
    write(46) c_time,c_dx,c_dy, (c_z(k),k=1,nzm),(c_p(k),k=1,nzm)

  end if ! save3Dbin

end if ! masterproc

! output ice nucleation rate
  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=rnfreez(i,j,k)*1.e6
    end do
   end do
  end do
  name='RNFR'
  long_name='Bigg freezing rate'
  units='#/m3/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)


  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=rnicr(i,j,k)*1.e6
    end do
   end do
  end do
  name='RNIC'
  long_name='Ice nucleation rate'
  units='#/m3/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  call task_barrier()

  if(nfields.ne.nfields1) then
    if(masterproc) print*,'write_nucleation error: nfields'
    call task_abort()
  end if
  if(masterproc) then
    close (46)
    if(RUN3D.or.save3Dsep) then
       if(dogzip3D) call systemf('gzip -f '//filename)
       print*, 'Writting nucleation data. file:'//filename
    else
       print*, 'Appending nucleation data. file:'//filename
    end if
  endif

 
end
