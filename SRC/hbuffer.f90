module hbuffer

use grid
use params

implicit none

	integer hbuf_length, HBUF_MAX_LENGTH
	parameter(HBUF_MAX_LENGTH = 1000)
	
	real(8) hbuf(HBUF_MAX_LENGTH*nzm+30)	
	character *8 namelist(HBUF_MAX_LENGTH)	
	character *80 deflist(HBUF_MAX_LENGTH)	
	character *10 unitlist(HBUF_MAX_LENGTH)	
	integer status(HBUF_MAX_LENGTH)
	integer average_type(HBUF_MAX_LENGTH)

CONTAINS

!------------------------------------------------------------

subroutine hbuf_average(n)

!       Average the profiles in buffer

use vars        

implicit none
integer l, n 
real(8) coef

coef=1./dble(n)
do l = 1,hbuf_length*nzm
 hbuf(l) = hbuf(l) *coef
end do  

end subroutine hbuf_average

!------------------------------------------------------------


subroutine hbuf_avg_put(name, f, dimx1,dimx2,dimy1,dimy2,dimz,factor)

!       Write to buffer an averaged 3D field (not to file yet)

implicit none

integer dimx1, dimx2, dimy1, dimy2, dimz
real f(dimx1:dimx2, dimy1:dimy2, dimz), fz(nzm), factor
character *(*) name
integer l, begin, k
logical flag

flag=.false.
do  l = 1,hbuf_length
  if(.not.(lgt(name,namelist(l)).or.llt(name,namelist(l)))) then
     flag=.true.
     if(status(l).gt.0) then
         status(l) = 999
         call averageXY(f,dimx1,dimx2,dimy1,dimy2,dimz,fz)
         begin = (l-1)*nzm
         do k = 1,nzm
           hbuf(begin+k) = hbuf(begin+k) + fz(k) * factor
         end do
     endif
  endif
end do
if(.not.flag.and.masterproc) print*,name

end subroutine hbuf_avg_put

!------------------------------------------------------------

subroutine hbuf_flush

!       Flush the buffer


use vars
implicit none
integer l,k,n

do l=1,hbuf_length*nzm
  hbuf(l) = 0.
end do
do n=1,ncondavg
   do k=1,nzm
      condavg_factor(k,n) = 0.
   end do
end do

s_acld=0.
s_acldl=0.
s_acldm=0.
s_acldh=0.
s_acldcold=0.
s_acldisccp=0.
s_acldlisccp=0.
s_acldmisccp=0.
s_acldhisccp=0.
s_acldmodis=0.
s_acldlmodis=0.
s_acldmmodis=0.
s_acldhmodis=0.
s_acldliqmodis=0.
s_acldicemodis=0.
s_acldmisr=0.
s_relmodis=0.
s_reimodis=0.
s_lwpmodis=0.
s_iwpmodis=0.
s_tbisccp=0.
s_tbclrisccp=0.
s_cldtauisccp=0.
s_cldalbisccp=0.
s_cldtaumodis=0.
s_cldtaulmodis=0.
s_cldtauimodis=0.
s_ptopisccp=0.
s_ptopmodis=0.
s_ztopmisr=0.
s_ar=0.
s_arthr=0.
w_max=0.
u_max=0.
s_flnt=0.
s_flntoa=0.
s_flntoac=0.
s_flns=0.
s_flnsc=0.
s_flds=0.
s_fsnt=0.
s_fsntoa=0.
s_fsntoac=0.
s_fsns=0.
s_fsnsc=0.
s_fsds=0.
s_solin=0.
lhobs=0.
shobs=0.
s_sst = 0.
z_inv=0.
z2_inv=0.
z_ct=0.
z_cb=z(nzm)
z_ctmn=0.
z_cbmn=0.
z2_cb=0. 
z2_ct=0. 
cwpmean=0.
cwp2=0.
prec2=0.
precmax=0.
precmean=0.
precmean=0.
ncmn=0.
nrmn=0.
ncloudy=0.
nrainy=0.

zinv=0.
zctl=0.
zcbl=0.
lwpa=0.
iwpa=0.
lwp_var=0.
iwp_var=0.

end subroutine hbuf_flush

!------------------------------------------------------------
subroutine hbuf_init

!       Read list of vertical profile names to store in file

use grid, only: case, masterproc
use tracers, only: tracers_hbuf_init

implicit none
character *8 nm
character *80 def
character *10 un
integer stat,count,type,i,trcount
integer, external :: lenstr
character*3 filestatus
integer n,m
logical duplicate_entries

if (masterproc) write(*,*) ' lst filepath', trim(rundatadir)  ! MO 4/13/16

open(66,file=trim(rundatadir)//'/lst',&
                          status='old',form='formatted')

! first determine number of entries:

hbuf_length = 0
111    continue
read(66,err=111,end=222,fmt=*) stat,type,nm
if(stat.gt.0) hbuf_length = hbuf_length+1
goto 111
222    continue

if(hbuf_length.gt.HBUF_MAX_LENGTH) then
   print *,'Error: hbuf_length > HBUF_MAX_LENGTH.'
   call task_abort()
endif

! fill the buffers:

rewind(66)
count = 0
333    continue
read(66,err=333,end=444,fmt=*) stat,type,nm,def,un
if(stat.gt.0) then
   count = count + 1
   namelist(count) = nm
   deflist(count) = def
   unitlist(count) = un
   status(count) = stat
   average_type(count) = type
endif
goto 333
444    continue
trcount=0
call hbuf_conditionals_init(count,trcount)
hbuf_length = hbuf_length+trcount
trcount=0
if(dotracers) call tracers_hbuf_init(namelist,deflist,unitlist,status,average_type,count,trcount)
hbuf_length = hbuf_length+trcount
trcount=0
if(docloud.or.dosmoke) call hbuf_micro_init(namelist,deflist,unitlist,status,average_type,count,trcount)
hbuf_length = hbuf_length+trcount
trcount=0                          ! MO 4/13/16
if(dosgs) call hbuf_sgs_init(namelist,deflist,unitlist,status,average_type,count,trcount)
hbuf_length = hbuf_length+trcount

! check if there are dublicate entries in the stat list:
duplicate_entries = .false.
do n = 1,count-1
 do m = n+1,count
  if (trim(namelist(n)).eq.trim(namelist(m))) then
   duplicate_entries=.true.
   if(masterproc) then
    print*,'Error: Multiple definition of '//namelist(n)// ' variable in stat list'
   end if
  end if
 end do
end do

! Halt simulation if duplicate entries appear in namelist.
if(duplicate_entries) call task_abort()

if(masterproc) then
         print *,'Number of statistics profiles:', hbuf_length
         print *,'Statistics profiles to save:'
         write(*,'(8(a,3x))')(namelist(i),i=1,hbuf_length)

! make sure that the stat file doesn't exist if a new run to prevent
! accidental overwrite

  filestatus='old'
  if(nrestart.eq.0.or.nrestart.eq.2) then
    filestatus='new'
  end if

  open (55,file='./OUT_STAT/'// &
                  case(1:lenstr(case))//'_'// &
                  caseid(1:lenstr(caseid))//'.stat', &
                  status=filestatus,form='unformatted')
  close(55)

end if

close (66)

call hbuf_flush()

end subroutine hbuf_init


!-----------------------------------------------------------------

subroutine hbuf_put(name, f, factor)

!       Write to buffer (not to file yet)

use grid, only: masterproc
implicit none

real f(nzm), factor
character *(*) name
integer l, begin, k
logical flag

flag=.false.
do  l = 1,hbuf_length
   if(.not.(lgt(name,namelist(l)).or.llt(name,namelist(l)))) then
       flag=.true.
       if(status(l).gt.0) then
          status(l) = 999
          begin = (l-1)*nzm
          do k = 1,nzm
            hbuf(begin+k) = hbuf(begin+k) + f(k)*factor
          end do
       endif
   endif
end do
if(.not.flag.and.masterproc) print*,'name ', name,' is missing in "lst" file.'

end subroutine hbuf_put

!----------------------------------------------------------------

subroutine hbuf_write(n)

!       Write buffer to file

use vars
implicit none
integer l, k, n, ntape, length
real(8) coef,aver,factor
real tmp(nzm),tmp1(nzm),dummy
real(8) hbuf1(HBUF_MAX_LENGTH*nzm+100)
real(4) hbuf_real4(HBUF_MAX_LENGTH*nzm+100)
integer nbuf
real cloud_f(nzm,ncondavg), tmp_f(nzm,ncondavg) !bloss
integer ncond !bloss
integer nsteplast

real(4) dummy_r4                !MO 4/19/16
real(8) dummy_r8                !MO 4/19/16
character(8) dummy_ch8          !MO 4/19/16
character(10) dummy_ch10        !MO 4/19/16
character(80) dummy_ch80        !MO 4/19/16
integer dummy_int               !MO 4/19/16

data ntape/56/

aver=1./dble(n)
factor=1./dble(nx*ny)

if(dompi) then
  ! average condavg_factor across domains.  This will sum the
  ! weighting of all of the conditional statistics across the
  ! processors and allow us to normalize the statistics on each
  ! processor accordingly.  In the end, these normalized statistics
  ! on each processor will sum to give the value of the conditional
  ! statistic.

  ! sum across processors
  call task_sum_real(condavg_factor,cloud_f,nzm*ncondavg)

else
   cloud_f(1:nzm,1:ncondavg) = condavg_factor(1:nzm,1:ncondavg)
end if

! create normalization/weighting factor for conditional statistics
!   Here, cloud_f holds the sum of condavg_factor across all processors.
condavg_factor(:,:) = 0.
do ncond = 1,ncondavg
   do k=1,nzm
      if(ABS(cloud_f(k,ncond)).gt.EPSILON(1.)) then
         condavg_factor(k,ncond)=float(nsubdomains)*float(n)/cloud_f(k,ncond)
      end if
   end do
end do

! compute each processor's component of the total conditional average.
length = 0
do l = 1,hbuf_length
  if(status(l).eq. 999) then
    length = length+1
    do ncond = 1,ncondavg
       if(average_type(l).eq.ncond) then
          do k=1,nzm
             hbuf((l-1)*nzm+k) = hbuf((l-1)*nzm+k)*condavg_factor(k,ncond)
          end do
          do k=1,nzm
             ! insert a missing value if there are no samples of this
             !   conditional statistic at this level.
             if(ABS(cloud_f(k,ncond)).lt.EPSILON(1.)) hbuf((l-1)*nzm+k) = -9999.
          end do
       endif
    end do
  endif
end do

!  Get statistics buffer from different processes, add them together
!  and average

if(dompi) then

   tmp1(1)=w_max
   tmp1(2)=u_max
   tmp1(3)=precmax
   tmp1(4)=z_ct
   call task_max_real(tmp1,tmp,4)
   w_max=tmp(1)
   u_max=tmp(2)
   precmax=tmp(3)
   z_ct=tmp(4)
   tmp1(1)=z_cb
   call task_min_real(tmp1,tmp,1)
   z_cb=tmp(1)
   k=hbuf_length*nzm
   hbuf(k+1)=s_acld
   hbuf(k+2)=s_acldl
   hbuf(k+3)=s_acldm
   hbuf(k+4)=s_acldh
   hbuf(k+5)=s_acldcold
   hbuf(k+6)=s_ar
   hbuf(k+7)=s_flns
   hbuf(k+8)=s_flnt
   hbuf(k+9)=s_flntoa
   hbuf(k+10)=s_flnsc
   hbuf(k+11)=s_flntoac
   hbuf(k+12)=s_flds
   hbuf(k+13)=s_fsns
   hbuf(k+14)=s_fsnt
   hbuf(k+15)=s_fsntoa
   hbuf(k+16)=s_fsnsc
   hbuf(k+17)=s_fsntoac
   hbuf(k+18)=s_fsds
   hbuf(k+19)=s_solin
   hbuf(k+20)=s_acldisccp
   hbuf(k+21)=s_acldlisccp
   hbuf(k+22)=s_acldmisccp
   hbuf(k+23)=s_acldhisccp
   hbuf(k+24)=s_sst
   hbuf(k+25)=s_acldmodis
   hbuf(k+26)=s_acldlmodis
   hbuf(k+27)=s_acldmmodis
   hbuf(k+28)=s_acldhmodis
   hbuf(k+29)=s_acldmisr
   hbuf(k+30)=s_relmodis
   hbuf(k+31)=s_reimodis
   hbuf(k+32)=s_lwpmodis
   hbuf(k+33)=s_iwpmodis
   hbuf(k+34)=s_tbisccp
   hbuf(k+35)=s_tbclrisccp
   hbuf(k+36)=s_acldliqmodis
   hbuf(k+37)=s_acldicemodis
   hbuf(k+38)=s_cldtaumodis
   hbuf(k+39)=s_cldtaulmodis
   hbuf(k+40)=s_cldtauimodis
   hbuf(k+41)=s_cldtauisccp
   hbuf(k+42)=s_cldalbisccp
   hbuf(k+43)=s_ptopisccp
   hbuf(k+44)=s_ptopmodis
   hbuf(k+45)=s_ztopmisr
   hbuf(k+46)=z_inv
   hbuf(k+47)=z2_inv
   hbuf(k+48)=ncloudy
   if(ncloudy.gt.0) then
      coef = 1./dble(ncloudy) 
      ncloudy = 1.
   else
      coef = 0.
   end if
   hbuf(k+49)=z_cbmn*coef
   hbuf(k+50)=z2_cb*coef
   hbuf(k+51)=z_ctmn*coef
   hbuf(k+52)=z2_ct*coef
   hbuf(k+53)=cwpmean
   hbuf(k+54)=cwp2
   hbuf(k+55)=precmean
   hbuf(k+56)=prec2
   hbuf(k+57)=ncmn*coef
   hbuf(k+58)=nrainy
   if(nrainy.gt.0) then
      coef = 1./dble(nrainy) 
      nrainy = 1.
   else
      coef = 0.
   end if
   hbuf(k+59)=nrmn*coef
   hbuf(k+60)=s_arthr

   nbuf = k+60

!Add for ISDAC (MO: mikhail@pnnl.gov)
   hbuf(k+61)=zinv
   hbuf(k+62)=zctl
   hbuf(k+63)=zcbl
   hbuf(k+64)=lwpa
   hbuf(k+65)=iwpa
   hbuf(k+66)=lwp_var
   hbuf(k+67)=iwp_var
   nbuf = k+67
!END Add for ISDAC

   call task_sum_real8(hbuf,hbuf1,nbuf)
   coef = 1./dble(nsubdomains)
   hbuf(1:nbuf) = hbuf1(1:nbuf)*coef

   s_acld=hbuf(k+1)
   s_acldl=hbuf(k+2)
   s_acldm=hbuf(k+3)
   s_acldh=hbuf(k+4)
   s_acldcold=hbuf(k+5)
   s_ar=hbuf(k+6)
   s_flns=hbuf(k+7)
   s_flnt=hbuf(k+8)
   s_flntoa=hbuf(k+9)
   s_flnsc=hbuf(k+10)
   s_flntoac=hbuf(k+11)
   s_flds=hbuf(k+12)
   s_fsns=hbuf(k+13)
   s_fsnt=hbuf(k+14)
   s_fsntoa=hbuf(k+15)
   s_fsnsc=hbuf(k+16)
   s_fsntoac=hbuf(k+17)
   s_fsds=hbuf(k+18)
   s_solin=hbuf(k+19)
   s_acldisccp=hbuf(k+20)
   s_acldlisccp=hbuf(k+21)
   s_acldmisccp=hbuf(k+22)
   s_acldhisccp=hbuf(k+23)
   s_sst=hbuf(k+24)
   s_acldmodis=hbuf(k+25)
   s_acldlmodis=hbuf(k+26)
   s_acldmmodis=hbuf(k+27)
   s_acldhmodis=hbuf(k+28)
   s_acldmisr=hbuf(k+29)
   s_relmodis=hbuf(k+30)
   s_reimodis=hbuf(k+31)
   s_lwpmodis=hbuf(k+32)
   s_iwpmodis=hbuf(k+33)
   s_tbisccp=hbuf(k+34)
   s_tbclrisccp=hbuf(k+35)
   s_acldliqmodis=hbuf(k+36)
   s_acldicemodis=hbuf(k+37)
   s_cldtaumodis=hbuf(k+38)
   s_cldtaulmodis=hbuf(k+39)
   s_cldtauimodis=hbuf(k+40)
   s_cldtauisccp=hbuf(k+41)
   s_cldalbisccp=hbuf(k+42)
   s_ptopisccp=hbuf(k+43)
   s_ptopmodis=hbuf(k+44)
   s_ztopmisr=hbuf(k+45)
   z_inv=hbuf(k+46)
   z2_inv=hbuf(k+47)
   if(hbuf(k+48).gt.0) then
     coef=1./hbuf(k+48)
   else
     coef=1.
   end if
   z_cbmn=hbuf(k+49)*coef
   z2_cb=hbuf(k+50)*coef
   z_ctmn=hbuf(k+51)*coef
   z2_ct=hbuf(k+52)*coef
   cwpmean=hbuf(k+53)
   cwp2=hbuf(k+54)
   precmean=hbuf(k+55)
   prec2=hbuf(k+56)
   ncmn=hbuf(k+57)*coef
   if(hbuf(k+58).gt.0) then
     coef=1./hbuf(k+58)
   else
     coef=1.
   end if
   nrmn=hbuf(k+59)*coef
   s_arthr=hbuf(k+60)

!Add for ISDAC (MO: mikhail@pnnl.gov)
   zinv=hbuf(k+61)
   zctl=hbuf(k+62)
   zcbl=hbuf(k+63)
   lwpa=hbuf(k+64)
   iwpa=hbuf(k+65)
   lwp_var=hbuf(k+66)
   iwp_var=hbuf(k+67)
!END Add for ISDAC


   z2_inv=z2_inv*factor*aver-(z_inv*factor*aver)**2
   z2_cb=z2_cb*factor-(z_cbmn*factor)**2
   z2_ct=z2_ct*factor-(z_ctmn*factor)**2
   cwp2=cwp2*factor*aver-(cwpmean*factor*aver)**2
   prec2=prec2*factor*aver-(precmean*factor*aver)**2

endif
if(masterproc) then

  open (ntape,file='./OUT_STAT/'// &
                  case(1:lenstr(case))//'_'// &
                  caseid(1:lenstr(caseid))//'.stat', &
                  status='unknown',form='unformatted')
!mo  if(nstep.ne.nstat) then
  if(nstep.ne.nstat .and. nstep.ne.0) then  ! Added nstep.ne.0 for writing initial
!  if(nstep.ne.0) then           !MO 4/13/16 ! Added nstep.ne.0 for writing initial
                                            ! profiles for ISDAC (MO:mikhail@pnnl.gov)
    do while(.true.)
!MO      read(ntape, end=222) dummy
      read(ntape, end=222) dummy_ch80
!      print*, '  caseID=',dummy_ch80    !MO 4/13/16
!      print*, '  nstep=',nstep,' hbuf_length =',hbuf_length    !MO 4/13/16

      read(ntape)  dummy,dummy,nsteplast
      print*,' nsteplast= ',nsteplast                          !MO 4/13/16

      if(nsteplast.ge.nstep) then
        backspace(ntape)
        backspace(ntape)  ! these two lines added because of 
!MO        read(ntape)       ! a bug in gfrotran compiler
      read(ntape, end=222) dummy_ch80
!      print*, '  caseID1=',dummy_ch80    !MO 4/13/16
!         print*,'stat file at nstep ',nsteplast
        goto 222   ! yeh, I know, it's bad ....
      end if    ! (nsteplast.ge.nstep)
!MO      read(ntape) 
      read(ntape)  dummy_int
!      print*,' dummy_int (length)= ',dummy_int  !MO 4/13/16

      do l = 1,hbuf_length
         if(status(l).eq. 999) then
!MO           read(ntape) 
!MO           read(ntape) 
!MO           read(ntape) 
!MO           read(ntape)
           read(ntape) dummy_ch8
!           print*,' l=',l,'  dummy_ch8=',dummy_ch8       !MO 4/13/16
           read(ntape) dummy_ch80
!           print*,'      dummy_ch80=',dummy_ch80     !MO 4/13/16
           read(ntape) dummy_ch10
!           print*,'      dummy_ch10=',dummy_ch10     !MO 4/13/16 
           read(ntape) dummy_r4
!           print*,'      dummy_r4=',dummy_r4         !MO 4/13/16

        end if
      end do
    end do
222 continue
     backspace(ntape)
     backspace(ntape)
     read(ntape)
  endif

  print *,'Writting history file ',caseid
  write(ntape)  caseid, version
  write(ntape)  day-nstat*dt/2./86400.,dt,nstep,nx,ny,nz,nzm, &
                dx,dy,dz,adz,z,pres,real(s_sst*aver*factor)+t00,pres0, &
                real(s_acld*aver*factor),real(s_ar*aver*factor), &
                real(s_acldcold*aver*factor),w_max,u_max,-1.,&
                real(s_flns*aver*factor),real(s_flnt*aver*factor),real(s_flntoa*aver*factor),&
                real(s_flnsc*aver*factor),real(s_flntoac*aver*factor),real(s_flds*aver*factor), &
                real(s_fsns*aver*factor),real(s_fsnt*aver*factor),real(s_fsntoa*aver*factor), &
                real(s_fsnsc*aver*factor),real(s_fsntoac*aver*factor), &
                real(s_fsds*aver*factor),real(s_solin*aver*factor), &
                real(sstobs),real(lhobs*aver),real(shobs*aver), &
                real(s_acldl*aver*factor),real(s_acldm*aver*factor),real(s_acldh*aver*factor), &
                s_acldisccp,s_acldlisccp,s_acldmisccp,s_acldhisccp, &
                s_acldmodis,s_acldlmodis,s_acldmmodis,s_acldhmodis,s_acldmisr, &
                s_relmodis, s_reimodis, s_lwpmodis, s_iwpmodis, s_tbisccp, s_tbclrisccp, &
                s_acldliqmodis, s_acldicemodis, s_cldtauisccp, s_cldalbisccp,  s_ptopisccp, &
                s_cldtaumodis, s_cldtaulmodis, s_cldtauimodis, s_ptopmodis, s_ztopmisr, &
                real(z_inv*aver*factor), z2_inv, &
                real(z_ctmn*factor), z2_ct, z_ct, &
                real(z_cbmn*factor), z2_cb, z_cb,&
                real(cwpmean*aver*factor), cwp2, &
                real(precmean*aver*factor), prec2, precmax, &
!mo                ncmn, nrmn, real(s_arthr*aver*factor) 
                ncmn, nrmn, real(s_arthr*aver*factor), &          ! mo
                zinv, zctl, zcbl, lwpa, iwpa, lwp_var, iwp_var    ! mo
 
  write(ntape) length
  
!  print*,' hbuf_length =',hbuf_length              !MO 4/13/16
  hbuf_real4(1:hbuf_length*nzm) = hbuf(1:hbuf_length*nzm)
  do l = 1,hbuf_length
     if(status(l).eq. 999) then
         write(ntape) namelist(l)
         write(ntape) deflist(l)
         write(ntape) unitlist(l)
         write(ntape) (hbuf_real4((l-1)*nzm+k),k=1,nzm)
     end if
  end do
  close (ntape)

end if

end subroutine hbuf_write



end module hbuffer
