     
subroutine write_fields3D_micro
	
use vars
use microphysics
 
implicit none
character *80 filename,long_name
character *8 name
character *10 timechar
character *4 rankchar
character *6 filetype
character *10 units
character *10 c_z(nzm),c_p(nzm),c_rho(nzm),c_dx, c_dy, c_time
character *10 c_rn(ncn),c_rnr(ncd)
integer i,j,k,nfields,nfields1
real tmp(nx,ny,nzm)

character *3 binumber
integer m


nfields= 16+1+3*33 ! number of 3D fields to save

! if(.not.docloud) nfields=nfields-1
! if(.not.doprecip) nfields=nfields-1
nfields1=0

if(masterproc) then

  write(rankchar,'(i4)') nsubdomains
  write(timechar,'(i10)') nstep
  do k=1,11-lenstr(timechar)-1
    timechar(k:k)='0'
  end do

  if(RUN3D) then
    if(save3Dbin) then
      filetype = '.bin3D'
    else
      filetype = '.com3D'
    end if
    filename='./OUT_3D/'//trim(case)//'_'//trim(caseid)//'_micro_'// &
        rankchar(5-lenstr(rankchar):4)//'_'//timechar(1:10)//filetype
    open(46,file=filename,status='unknown',form='unformatted')

  else
    if(save3Dbin) then
     if(save3Dsep) then
       filetype = '.bin3D'
     else
       filetype = '.bin2D'
     end if
    else
     if(save3Dsep) then
       filetype = '.com3D'
     else
       filetype = '.com2D'
     end if
    end if
    if(save3Dsep) then
      filename='./OUT_3D/'//trim(case)//'_'//trim(caseid)//'_micro_'// &
        rankchar(5-lenstr(rankchar):4)//'_'//timechar(1:10)//filetype
      open(46,file=filename,status='unknown',form='unformatted')	
    else
      filename='./OUT_3D/'//trim(case)//'_'//trim(caseid)//'_micro_'// &
        rankchar(5-lenstr(rankchar):4)//filetype
      if(nrestart.eq.0.and.notopened3D) then
         open(46,file=filename,status='unknown',form='unformatted')	
      else
         open(46,file=filename,status='unknown', &
                              form='unformatted', position='append')
      end if
      notopened3D=.false.
    end if  

  end if

  if(save3Dbin) then

    write(46) nx,ny,nzm,ncn,ncd,   &
            nsubdomains,nsubdomains_x,nsubdomains_y,nfields
    do k=1,nzm
      write(46) z(k) 
    end do
    do k=1,nzm
      write(46) pres(k)
    end do
    do k=1,nzm
      write(46) rho(k)
    end do
    write(46) dx
    write(46) dy
    write(46) nstep*dt/(3600.*24.)+day0
    do k=1,ncn
      write(46) rccn(k)
    end do
    do k=1,ncn
      write(46) DROPRADII(k)
    end do

  else

    write(long_name,'(10i4)') nx,ny,nzm,ncn,ncd, &
             nsubdomains,nsubdomains_x,nsubdomains_y,nfields
    do k=1,nzm
       write(c_z(k),'(f10.3)') z(k)
    end do
    do k=1,nzm
       write(c_p(k),'(f10.3)') pres(k)
    end do
    do k=1,nzm
       write(c_rho(k),'(f10.3)') rho(k)
    end do
    write(c_dx,'(f10.0)') dx
    write(c_dy,'(f10.0)') dy
    write(c_time,'(f10.4)') nstep*dt/(3600.*24.)+day0
    do k=1,ncn
      write(c_rn(k),'(e10.4)') rccn(k)
    end do
    do k=1,ncn
      write(c_rnr(k),'(e10.4)') DROPRADII(k)
    end do
	
    write(46) long_name(1:40)
    write(46) c_time,c_dx,c_dy,                                    &
         (c_z(k),k=1,nzm),(c_p(k),k=1,nzm),(c_rho(k),k=1,nzm)
    write(46)                                                      &
         (c_rn(k),k=1,ncn),(c_rnr(k),k=1,ncn)
  end if ! save3Dbin

end if ! masterproc

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qt(i,j,k)*1000.*rho(k)
    end do
   end do
  end do
  name='QT'
  long_name='Total Water Vaper Content'
  units='g/m^3'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=ssatw(i,j,k)
    end do
   end do
  end do
  name='SSW'
  long_name='Saturation ratio (water)'
  units=' '
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=ssati(i,j,k)
    end do
   end do
  end do
  name='SSI'
  long_name='Saturation ratio (ice)'
  units=' '
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qc(i,j,k)*1000.*rho(k)
    end do
   end do
  end do
  name='QL'
  long_name='Cloud Water Content'
  units='g/m^3'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qr(i,j,k)*1000.*rho(k)
    end do
   end do
  end do
  name='QR'
  long_name='Rain Water Content'
  units='g/m^3'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qi(i,j,k)*1000.*rho(k)
    end do
   end do
  end do
  name='QI'
  long_name='Ice crystal Content'
  units='g/m^3'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qs(i,j,k)*1000.*rho(k)
    end do
   end do
  end do
  name='QS'
  long_name='Snow Content'
  units='g/m^3'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qg(i,j,k)*1000.*rho(k)
    end do
   end do
  end do
  name='QG'
  long_name='Graupel Content'
  units='g/m^3'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qh(i,j,k)*1000.*rho(k)
    end do
   end do
  end do
  name='QH'
  long_name='Hail Content'
  units='g/m^3'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qna(i,j,k)
    end do
   end do
  end do
  name='NCN'
  long_name='CCN Number Concentration'
  units='#/cm3'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qnin(i,j,k)*1.e3
    end do
   end do
  end do
  name='NIN'
  long_name='IN Number Concentration'
  units='#/L'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qnc(i,j,k)
    end do
   end do
  end do
  name='NCD'
  long_name='Droplet Number Concentration'
  units='#/cm3'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qnr(i,j,k)
    end do
   end do
  end do
  name='NCR'
  long_name='Drop Number Concentration'
  units='#/cm3'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qni(i,j,k)
    end do
   end do
  end do
  name='NCI'
  long_name='Ice Crystal Number Concentration'
  units='#/cm3'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qns(i,j,k)*1000.0
    end do
   end do
  end do
  name='NCS'
  long_name='Snow Number Concentration'
  units='#/L'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qng(i,j,k)*1000.0
    end do
   end do
  end do
  name='NCG'
  long_name='Graupel Number Concentration'
  units='#/L'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qnh(i,j,k)*1000.0
    end do
   end do
  end do
  name='NCH'
  long_name='Hail Number Concentration'
  units='#/L'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

! Write CCN spectrum

  do m=1,ncn

    write(binumber,'(i3)') m
    do k=1,3-lenstr(binumber)
      binumber(k:k)='0'
    end do

    nfields1=nfields1+1
    do k=1,nzm
     do j=1,ny
      do i=1,nx
        tmp(i,j,k)=fncn(i,j,k,m)/(xccn(m)*1.e-3)*rho(k)*col*1.e-6
      end do
     end do
    end do
    name='FNCN'//binumber(1:3)
    long_name='CCN Number Concentration bin '//binumber(1:3)
    units='#/cm3/bin'
    call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
  end do    ! m

! Write dropelt spectrum

  do m=1,ncd

    write(binumber,'(i3)') m
    do k=1,3-lenstr(binumber)
      binumber(k:k)='0'
    end do

    nfields1=nfields1+1
    do k=1,nzm
     do j=1,ny
      do i=1,nx
!        tmp(i,j,k)=fmcd(i,j,k,m)*1000.                ! g/kg/bin
        tmp(i,j,k)=ffcd(i,j,k,m)/(xl(m)*1.e-3)*rho(k)*col*1.e-6   ! #/cm3/bin
      end do
     end do
    end do
    name='FFCD'//binumber(1:3)
    long_name='Droplet Number Concentration per bin '//binumber(1:3)
!    long_name='Liquid Water Mixing Ratio bin '//binumber(1:3)
!    units='g/kg/bin'
    units='#/cm3/bin'
    call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
  end do    ! m

! Write ice partical spectrum

  do m=1,ncd

    write(binumber,'(i3)') m
    do k=1,3-lenstr(binumber)
      binumber(k:k)='0'
    end do

    nfields1=nfields1+1
    do k=1,nzm
     do j=1,ny
      do i=1,nx
!        tmp(i,j,k)=fmcd(i,j,k,m)*1000.                ! g/kg/bin
        tmp(i,j,k)=(ffic(i,j,k,m)/(Xi(m,1)*1.e-3)*rho(k)*col +    &
   &                ffip(i,j,k,m)/(Xi(m,2)*1.e-3)*rho(k)*col +    &
   &                ffiD(i,j,k,m)/(Xi(m,3)*1.e-3)*rho(k)*col +    &
   &                ffsn(i,j,k,m)/(XS(m)*1.e-3)*rho(k)*col +      &
   &                ffgl(i,j,k,m)/(XG(m)*1.e-3)*rho(k)*col +      &
   &                ffhl(i,j,k,m)/(Xh(m)*1.e-3)*rho(k)*col)*1.0e-3        ! #/L/bin
      end do
     end do
    end do
    name='FFIN'//binumber(1:3)
    long_name='Ice partical Number Concentration per bin '//binumber(1:3)
!    long_name='Liquid Water Mixing Ratio bin '//binumber(1:3)
!    units='g/kg/bin'
    units='#/L/bin'
    call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
  end do    ! m

  call task_barrier()
  print*, nfields, nfields1

  if(nfields.ne.nfields1) then
    if(masterproc) print*,'write_fields3D_micro error: nfields'
    call task_abort()
  end if
  if(masterproc) then
    close (46)
    if(RUN3D.or.save3Dsep) then
       if(dogzip3D) call systemf('gzip -f '//filename)
       print*, 'Writting 3D data. file:'//filename
    else
       print*, 'Appending 3D data. file:'//filename
    end if
  endif

 
end subroutine write_fields3D_micro
