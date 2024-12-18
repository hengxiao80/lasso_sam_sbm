module microphysics

! This module adapts the SBM from wrf, originally from A. Khain. It should be located in the
! MICRO_HUJISBM directory.
! The following flags must be set
! docloud 	= .true.,
! doprecip 	= .true.,  

! 'grid' is SAM module which contains the required grid information

use grid, only: nx,ny,nzm,nz, masterproc,RUN3D, &  ! grid dimensions; nzm=nz-1 - # of levels for all scalars
              & dimx1_s,dimx2_s,dimy1_s,dimy2_s ! actual scalar-array dimensions

use params, only: cp, ggr, rgas,rv,lsub 
! v6/9/4
use params, only: doprecip, docloud

use module_hujisbm

implicit none

! Allocate the required memory for all the prognostic microphysics arrays:

real micro_field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, nmicro_fields)

! We assume that our prognostic variables are alligned as follows:


! 1                      - water mixing ratio
! 2       -- (ncn+1)     - CCN concentrations for NCN bins
! (ncn+2) -- (ncd+ncd+1) - liquid water mixing ratio in NCD bins 

! For many reasons, for example, to compute water budget, we may need to
! know which variables among the prognostic ones represent water mixing ratio, regardless
! of water species. We use a simple array of flags, with 1 marking the water mass
! variable:

!integer, parameter :: flag_wmass(nmicro_fields) = (/0                           &
!        &       ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0                          &
!        &       ,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/) 
integer, parameter :: flag_wmass(nmicro_fields) = (/1,                          &
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,         &
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,         &
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,         &
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,         &
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,         &
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,         &
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,         &
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,         &
! add for IN
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)

! To implement large-scale forcing, surface fluxes, etc, SAM needs to know
! which variable has a water vapor information. In our example, it is variable #3:

integer, parameter :: index_water_vapor = 1 ! index for variable that contains water vapor

! Now, we need to specify which variables describe precipitation. This is needed because
! SAM has two logical flags to deal with microphysics proceses - docloud and doprecip.
! docloud set to 'true' means that condensation/sublimation processes are allowed to
! form clouds. However, the possibility of rain, snow, heil, etc., is controled by
! a second flag: doprecip. If doprecip=.false. than no precipitation is allowed, hence 
! no advection, diffusion, and fallout of corresponding variables should be done; 
! therefore, SAM needs an array of flags that mark the prognostic variables which
! only make sense when doprecip=.true. :

!integer, parameter :: flag_precip(nmicro_fields) = (/0,                         &
!     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,         &
!     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,         &
!     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,         &
!     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,         &
!     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,         &
!     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,         &
!     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,         &
!     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,         &
!! add for IN
!     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)

!26 oct 2009: added sedimentation for all liquid and ice particles (MO)
integer, parameter :: flag_precip(nmicro_fields) = (/0,                         &
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,         &
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,         &
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,         &
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,         &
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,         &
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,         &
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,         &
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,         &
! add for IN
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)


! Sometimes, cloud ice (or even cloud water) is allowed to be a subject of
! gravitational sedimentation, usually quite slow compared to the precipitation
! drops. SAM calls a special routine, ice_fall() that computes sedimentation of cloud ice.
! However, it is a rudiment from SAM's original single-moment microphysics.
! Instead, you may want to handle sedimentation of cloud water/ice yourself similarly
! to precipitation variables. In this case, set the index for falling cloud ice to -1, which
! means that no default ice mixing ratio sedimentation is done. 

integer, parameter :: index_cloud_ice = -1   ! index for cloud ice (sedimentation)

! The following arrays are needed to set the turbulent surface and top fluxes 
! for the microphysics prognostic variables:

real fluxbmk (nx, ny, 1:nmicro_fields) ! surface fluxes
real fluxtmk (nx, ny, 1:nmicro_fields) ! top boundary fluxes 

!!! these arrays are needed for output statistics from advection and diffusion routines:

real mkwle(nz,1:nmicro_fields)  ! resolved vertical flux
real mkwsb(nz,1:nmicro_fields)  ! SGS vertical flux
real mkadv(nz,1:nmicro_fields)  ! tendency due to vertical advection
real mkdiff(nz,1:nmicro_fields) ! tendency due to vertical diffusion
! add for version 6.7.3
real mklsadv(nz,1:nmicro_fields) ! tendency due to large-scale vertical advection
!------------------------------------------------------------------

! It would be quite inconvenient to work with the micro_field array itself. Besides,
! your original microphysics routines use some specific names for the prognostic variables
! that you don't wanna change. Therefore, you need a way of using portions of the 
! make aliases for prognostic variables. You can use the aliases as you would ordinary arrays.

!--------------------------------------------------------------------
! prognostic variables:

real qt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm)       ! total water mixing ratio [kg/kg]
real fncn(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm,ncn) ! CCN number distribution  
                                                   !              function [#/kg/bin]
real ffcd(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm,ncd) !
real ffic(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm,ncd) !
real ffip(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm,ncd) !
real ffid(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm,ncd) !
real ffsn(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm,ncd) !
real ffgl(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm,ncd) !
real ffhl(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm,ncd) !
! add for in
real ffin(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm,ncd) !

! map prognostic variables onto micro_field array:

equivalence (qt(dimx1_s,dimy1_s,1),micro_field(dimx1_s,dimy1_s,1,1))
equivalence (fncn(dimx1_s,dimy1_s,1,1),micro_field(dimx1_s,dimy1_s,1,2))
equivalence (ffcd(dimx1_s,dimy1_s,1,1),micro_field(dimx1_s,dimy1_s,1,ncn+2))
equivalence (ffic(dimx1_s,dimy1_s,1,1),micro_field(dimx1_s,dimy1_s,1,ncn*2+2))
equivalence (ffip(dimx1_s,dimy1_s,1,1),micro_field(dimx1_s,dimy1_s,1,ncn*3+2))
equivalence (ffid(dimx1_s,dimy1_s,1,1),micro_field(dimx1_s,dimy1_s,1,ncn*4+2))
equivalence (ffsn(dimx1_s,dimy1_s,1,1),micro_field(dimx1_s,dimy1_s,1,ncn*5+2))
equivalence (ffgl(dimx1_s,dimy1_s,1,1),micro_field(dimx1_s,dimy1_s,1,ncn*6+2))
equivalence (ffhl(dimx1_s,dimy1_s,1,1),micro_field(dimx1_s,dimy1_s,1,ncn*7+2))
!add for in
equivalence (ffin(dimx1_s,dimy1_s,1,1),micro_field(dimx1_s,dimy1_s,1,ncn*8+2))

!--------------------------------------------------------------------
! diagnostic variables (3D arrays):

!real qv(nx, ny, nzm)   ! water vapor mixing ratio
real qc(nx, ny, nzm)   ! liquid water mixing ratio    [kg/kg]
real qr(nx, ny, nzm)   ! rain/drizzle mixing ratio    [kg/kg]
real qi(nx, ny, nzm)   ! total ice mixing ratio    [kg/kg]
real qic(nx, ny, nzm)  ! ice1 mixing ratio    [kg/kg]
real qip(nx, ny, nzm)  ! ice2
real qid(nx, ny, nzm)  ! ice3
real qs(nx, ny, nzm)   ! snow mixing ratio    [kg/kg]
real qg(nx, ny, nzm)   ! graupel mixing ratio    [kg/kg]
real qh(nx, ny, nzm)   ! hail mixing ratio    [kg/kg]
real qna(nx, ny, nzm)  ! CCN number concentration [/cm^3]
! add for in
real qnin(nx, ny, nzm)  ! IN number concentration [/cm^3]
!
real qnc(nx, ny, nzm)  ! droplet number concentration [/cm^3]
real qnr(nx, ny, nzm)  ! rain drop number concentration [/cm^3]
real qni(nx, ny, nzm)  ! total ice number concentration [/L]
real qnic(nx, ny, nzm)  ! ice1 number concentration [/L]
real qnip(nx, ny, nzm)  ! ice2 number concentration [/L]
real qnid(nx, ny, nzm)  ! ice3 number concentration [/L]
real qns(nx, ny, nzm)  ! snow number concentration  [/L]
real qng(nx, ny, nzm)  ! graupel number concentration [/L]
real qnh(nx, ny, nzm)  ! hail number concentration [/L]
real rainnc(nx, ny)    ! accumulated total grid scale precipitation [/L] 
 
real ssatw(nx, ny, nzm) ! supersaturation over liquid water [%]
real ssati(nx, ny, nzm) ! supersaturation over ice water [%]

!v6.9.4
! effective radius for instrument simulators
real reffc(nx, ny, nzm)
real reffi(nx, ny, nzm)

! added for SHEA output variable vfice_mw
real vfice_mw(nx, ny, nzm) ! mass weighted ice fall velocity
real vfic(nx, ny, nzm)  ! mass*fall velocity
real vfip(nx, ny, nzm) 
real vfid(nx, ny, nzm) 
real vfs(nx, ny, nzm) 
real vfg(nx, ny, nzm) 
real diffui_tend(nx, ny, nzm)
real difful_tend(nx, ny, nzm)
real qlbf_sed(dimx1_s:dimx2_s,dimy1_s:dimy2_s, nzm)
real qlaf_sed(dimx1_s:dimx2_s,dimy1_s:dimy2_s, nzm)
real qibf_sed(dimx1_s:dimx2_s,dimy1_s:dimy2_s, nzm)
real qiaf_sed(dimx1_s:dimx2_s,dimy1_s:dimy2_s, nzm)
real sedl_tend(dimx1_s:dimx2_s,dimy1_s:dimy2_s, nzm)
real sedi_tend(dimx1_s:dimx2_s,dimy1_s:dimy2_s, nzm)
real frzl_tend(nx, ny, nzm)
real qvbf(nx, ny, nzm)
real qvaf(nx, ny, nzm)
real qvtend_adv(nx, ny, nzm)
real qltend_adv(nx, ny, nzm)
real qitend_adv(nx, ny, nzm)
real qvst(nx, ny, nzm)
real qvend(nx, ny, nzm)
real qvtend(nx, ny, nzm)
real qlst(nx, ny, nzm)
real qlend(nx, ny, nzm)
real qltend(nx, ny, nzm)
real qist(nx, ny, nzm)
real qiend(nx, ny, nzm)
real qitend(nx, ny, nzm)
! for qni100
real qni100(nx, ny, nzm)  ! total ice number concentration [/L]
real qnic100(nx, ny, nzm)  ! ice1 number concentration [/L]
real qnip100(nx, ny, nzm)  ! ice2 number concentration [/L]
real qnid100(nx, ny, nzm)  ! ice3 number concentration [/L]
real qns100(nx, ny, nzm)  ! snow number concentration  [/L]
real qng100(nx, ny, nzm)  ! graupel number concentration [/L]
real qnh100(nx, ny, nzm)  ! hail number concentration [/L]

!--------------------------------------------------------------------
! other variables:
real t_old(nx, ny, nzm) 
real qv_old(nx, ny, nzm) 
! J. Fan 08/2008
!real, parameter    :: ql_prec = 1.e-9   ! min ql for precipitation
real, parameter    :: ql_prec = 1.e-15   ! min ql for precipitation

double precision, parameter     :: ss_max = 0.003d+0    ! maximum supersaturation ratio for droplet activation
! double precision, parameter     :: ss_max = 0.03d+0    ! maximum supersaturation ratio for droplet activation

! add for nucleation rate (J. Fan)
      REAL rnfreez(nx, ny, nzm), rnicr(nx, ny, nzm)
      REAL  rndropr(nx, ny, nzm) ! drop evaporation freezing rate (=0 if evapor IN mechanism is used)
      REAL rnfr_hom(nx, ny, nzm), rnfr_imm(nx, ny, nzm)

! adpated from wrf-sbm interface
!-----------------------------------------------------------------------
!
      CONTAINS

!-----------------------------------------------------------------------
! add to be compatible with version 6.7.3
subroutine micro_setparm()
  use vars
  implicit none
end subroutine micro_setparm
                                                                                                        
subroutine micro_print()
  implicit none
  integer :: k
end subroutine micro_print

!-----------------------------------------------------------------------
      SUBROUTINE micro_proc()

      use vars, only: dudt,dvdt,dwdt,tabs,t,pres,rho,nrestart,qv,gamaz
      use grid, only: nc,dx,dy,dz,dt,nstep,icycle,z, &
                     & dimx1_w,dimx2_w, dimy1_w, dimy2_w 
      use params
      IMPLICIT NONE
!-----------------------------------------------------------------------
      INTEGER, PARAMETER :: its=1, ite=nx, jts=1,     &
                            jte=ny, kts=1, kte=nzm    
      INTEGER, PARAMETER :: hujisbm_unit1=22
!
!-----------------------------------------------------------------------
!     LOCAL VARS
!-----------------------------------------------------------------------

!     SOME VARS WILL BE USED FOR DATA ASSIMILATION (DON'T NEED THEM NOW). 
!     THEY ARE TREATED AS LOCAL VARS, BUT WILL BECOME STATE VARS IN THE 
!     FUTURE. SO, WE DECLARED THEM AS MEMORY SIZES FOR THE FUTURE USE

      REAL,  DIMENSION(kts:kte):: rhocgs, pcgs, zcgs

      INTEGER :: I,J,K,KFLIP
       REAL &
     &        SUP2_OLD, DSUPICEXZ,TFREEZ_OLD,DTFREEZXZ, &
     &        DTIME,DTCOND
! SBM VARIABLES
      REAL,DIMENSION (nkr) :: FF1IN,FF3IN,FF4IN,FF5IN,&
     &              FF1R,FF3R,FF4R,FF5R,FCCN
! Add for IN
      REAL :: fin(nkr)

      REAL,DIMENSION (nkr,icemax) :: FF2IN,FF2R
!!! for ccn regeneration
      real :: fccn0(nkr)
      real :: ndrop, subtot  ! for diagnostic CCN
!!!

! add for nucleation rate
      real rnfr, rnic,rndrop,frzfract
      real fr_hom, fr_imm
!
!mo add for SS limit 
      double precision del1in_limit
      DOUBLE PRECISION DEL1NR,DEL2NR,DEL12R,DEL12RD,ES1N,ES2N,EW1N,EW1PN
      DOUBLE PRECISION DELSUP1,DELSUP2,DELDIV1,DELDIV2
      DOUBLE PRECISION TT,QQ,TTA,QQA,PP,DPSA,DELTATEMP,DELTAQ
      DOUBLE PRECISION DIV1,DIV2,DIV3,DIV4,DEL1IN,DEL2IN,DEL1AD,DEL2AD
      REAL DEL_BB,DEL_BBN,DEL_BBR
      REAL FACTZ,CONCCCN_XZ,CONCDROP
       REAL SUPICE(KTE),AR1,AR2, &
     & DERIVT_X,DERIVT_Y,DERIVT_Z,DERIVS_X,DERIVS_Y,DERIVS_Z, &
     & ES2NPLSX,ES2NPLSY,EW1NPLSX,EW1NPLSY,UX,VX, &
! WX -add it by J Fan 08/2008
     & WX,  &
     & DEL2INPLSX,DEL2INPLSY,DZZ(KTE)
       INTEGER KRR
   
       REAL DTFREEZ_XYZ(ITE,JTE,KTE),DSUPICE_XYZ(ITE,JTE,KTE)

       REAL DXHUCM,DYHUCM
       REAL FMAX1,FMAX2,FMAX3,FMAX4,FMAX5
       INTEGER ISYM1,ISYM2,ISYM3,ISYM4,ISYM5
       INTEGER DIFFU
       REAL DTHALF
       REAL DELTAW
       integer imax,kmax
       real gmax
       real tmax,qmax,divmax,rainmax
       real qnmax,inmax,knmax
       real hydro
       real difmax,tdif,tt_old
       real teten,es
       real maxw
!!! For CCN regeration
       real tot_reg
       real totibf_diffu, totiaf_diffu,totibf_sed, totiaf_sed 
       real totlbf_frz, totlaf_frz, totlbf_diffu, totlaf_diffu, totlbf_sed, totlaf_sed
       real fccn_before(nkr), ff1r_before(nkr), fccnin(nkr)
       real totnum, totccn, totdrop,totnum_beforem,totccn_beforem, totdrop_beforem
       real ccnbefore, ccnafter

! For topfrz
       real tempin

!    For CCN verical distribution
       real factor

       integer itimestep, kr, ikl, ice, nkro, nkre 

!
       real inbeforenuc,inafternuc, icebeforenuc, iceafternuc,iceafterdiff

! v6.9.4
! for effcs calculation
       real top, bottom

!MO flag for coagulation (docoag=.false.  for ISDAC intercomparison)
       logical docoag
       docoag = .false.
  if (dt*(nstep-1).ge.tprecip) then !MO 4/19/16: For VOCALS intercomparison
       docoag = .true.              ! coagulation is ON after 1 h.
  endif                             ! tprecip (in s) is set in params.f90

       itimestep=nstep
!

      difmax = 0
!       print*,'itimestep = ',itimestep
!        if (itimestep.gt.150)return
        if (itimestep.eq.1.and.icycle.eq.1.and.masterproc)then
         if (iceprocs.eq.1)print*,'ICE PROCESES ACTIVE'
         if (iceprocs.eq.0)print*,'LIQUID PROCESES ONLY'
        end if
       tmax = 0
! nucleation timeteps
       IF (NCOND/2.eq.0)then
        dthalf = dt
!        print*,'dthalf not half = ',dthalf
       else
         dthalf = dt/2
!        print*,'dthalf really half = ',dthalf
       END IF
       DTCOND=DT/REAL(NCOND)
!       print*,'dt,dtcond = ',dt,dtcond

! print the height of the max velocity
!     if(mod(itimestep,30).eq.0) then
!      maxw = MAXVAL(dwdt(:,:,:,nc))
!      DO j = 1, ny
!      DO k = 1,nzm
!      DO i = 1,nx
!        if (dwdt(i,j,k,nc) .eq. maxw) print*,'W and height', maxw, z(k)
!      enddo
!      enddo
!      enddo
!     endif

!test
     qmax=0
     imax=0
     kmax=0
     DO j = jts,jte
     DO k = kts,kte
     DO i = its,ite
     if (qc(i,j,k).gt.qmax)imax=i
     if (qc(i,j,k).gt.qmax)kmax=k
     if (qc(i,j,k).gt.qmax)qmax=qc(i,j,k)
     end do
     end do
     end do
!test 
       if (itimestep.eq.1.and.icycle.eq.1.and.masterproc)then
             do kr = 1,nkr
              print*,'xl = ',xl(kr),vr1(kr),RLEC(kr),RO1BL(kr)
              print*,'xi = ',xi(kr,1),vr2(kr,1),RIEC(KR,1),RO2BL(KR,1)
              print*,'xi = ',xi(kr,2),vr2(kr,2),RIEC(KR,2),RO2BL(KR,2)
              print*,'xi = ',xi(kr,3),vr2(kr,3),RIEC(KR,3),RO2BL(KR,3)
              print*,'xs = ',xs(kr),vr3(kr),RSEC(kr),RO3BL(kr)
              print*,'xg = ',xg(kr),vr4(kr),RGEC(kr),RO4BL(kr)
              print*,'xh = ',xh(kr),vr5(kr),RHEC(kr),RO5BL(kr)
             end do
          print*,'tcrit = ',tcrit   
          print*,'ttcoal = ',ttcoal 

        end if


!       call micro_init()
      
!        OPEN (UNIT=1,FILE="eta_micro_lookup.dat",FORM="UNFORMATTED")
!         OPEN(UNIT=hujisbm_unit1,FILE="hujisbm_DATA",                  &
!    &        FORM="UNFORMATTED",STATUS="OLD",ERR=9061)
!
!         READ(hujisbm_unit1) VENTR1
!         READ(hujisbm_unit1) VENTR2
!         READ(hujisbm_unit1) ACCRR
!         READ(hujisbm_unit1) MASSR
!         READ(hujisbm_unit1) VRAIN
!         READ(hujisbm_unit1) RRATE
!         READ(hujisbm_unit1) VENTI1
!         READ(hujisbm_unit1) VENTI2
!         READ(hujisbm_unit1) ACCRI
!         READ(hujisbm_unit1) MASSI
!         READ(hujisbm_unit1) VSNOWI
!         READ(hujisbm_unit1) VEL_RF
!        read(hujisbm_unit1) my_growth    ! Applicable only for DTPH=180 s
!         CLOSE (hujisbm_unit1)
!       ENDIF
!
      DEL_BB=BB2_MY-BB1_MY
      DEL_BBN=BB2_MYN-BB1_MYN
      DEL_BBR=BB1_MYN/DEL_BBN
!
      DXHUCM=100.*DX
      DYHUCM=100.*DY
!      print*,'dxhucm = ',dxhucm
!      print*,'dyhucm = ',dyhucm

!
!-----------------------------------------------------------------------
!**********************************************************************
!-----------------------------------------------------------------------
!

!     MY_GROWTH(MY_T1:MY_T2)=MP_RESTART_STATE(MY_T1:MY_T2)
!
!     C1XPVS0=MP_RESTART_STATE(MY_T2+1)
!     C2XPVS0=MP_RESTART_STATE(MY_T2+2)
!     C1XPVS =MP_RESTART_STATE(MY_T2+3)
!     C2XPVS =MP_RESTART_STATE(MY_T2+4)
!     CIACW  =MP_RESTART_STATE(MY_T2+5)
!     CIACR  =MP_RESTART_STATE(MY_T2+6)
!     CRACW  =MP_RESTART_STATE(MY_T2+7)
!     CRAUT  =MP_RESTART_STATE(MY_T2+8) !
!     TBPVS(1:NX) =TBPVS_STATE(1:NX)
!     TBPVS0(1:NX)=TBPVS0_STATE(1:NX)
!
!      print*, its, ite, jts, jte, kts,kte

      DO k = kts,kte
      DO j = jts,jte
      DO i = its,ite

!=== 26oct2009 (MO)
! begin update qc, qr, qi, qs, qg, qh

        QC(I,j,k)=0
        QR(I,j,k)=0
        QI(I,j,k)=0
        QIC(I,j,k)=0
        QIP(I,j,k)=0
        QID(I,j,k)=0
        QS(I,j,k)=0
        QG(I,j,k)=0
        QH(I,j,k)=0

        DO  KR=1,NKR
          IF (KR.LT.KRDROP)THEN
            QC(I,j,k)=QC(I,j,k) &
     &        +ffcd(I,j,k,KR)
          ELSE
            QR(I,j,k)=QR(I,j,k) &
     &        +ffcd(I,j,k,KR)
          END IF
        END DO
        IF (QC(I,j,k).LT.1.E-8)QC(I,J,k)=0.0
        IF (QR(I,j,k).LT.1.E-8)QR(I,j,k)=0.0
        QC(I,j,k)= QC(I,j,k)* col
        QR(I,j,k)= qr(i,j,k)* col
   
        IF (ICEPROCS.EQ.1)THEN
         DO  KR=1,33
          QIC(I,j,k)=QIC(I,j,k) &
     &      +ffic(I,j,k,KR)
          QIP(I,j,k)=QIP(I,j,k) &
     &      +ffip(I,j,k,KR)
          QID(I,j,k)=QID(I,j,k) &
     &      +ffid(I,j,k,KR)
          QS(I,j,k)=QS(I,j,k) &
     &      +ffsn(I,j,k,KR)
          QG(I,j,k)=QG(I,j,k) &
     &      +ffgl(I,j,k,KR)
          QH(I,j,k)=QH(I,j,k) &
     &      +ffhl(I,j,k,KR)
         END DO
         QI(I,j,k) = QID(I,j,k)+QIP(I,j,k)+QIC(I,j,k)

         QI(I,j,k)= QI(I,j,k) * col 
         QIC(I,j,k)= QIC(I,j,k) * col
         QIP(I,j,k)= QIP(I,j,k) * col
         QID(I,j,k)= QID(I,j,k) * col
         QS(I,j,k)= QS(I,j,k) * col
         QG(I,j,k)= QG(I,j,k) * col
         QH(I,j,k)= QH(I,j,k) * col
        ENDIF

! end update qc, qr, qi, qs, qg, qh 
!=========================

! get rid of negative qv
        qv(i,j,k) = qt(i,j,k) - (qc(i,j,k)+qr(i,j,k)+qi(i,j,k)+qs(i,j,k)+    &
                   qg(i,j,k)+qh(i,j,k))
        if (qv(i,j,k) .le. 1.e-15) qv(i,j,k)=1.e-15

        IF (ICEPROCS.EQ.1)THEN
           tabs(i,j,k) = t(i,j,k) - gamaz(k) + fac_cond * (qc(i,j,k)+ qr(i,j,k)) &
   &                  + fac_sub * (qi(i,j,k)+ qs(i,j,k)+qg(i,j,k)+qh(i,j,k))
         ELSE
           tabs(i,j,k) = t(i,j,k) - gamaz(k) + fac_cond * (qc(i,j,k)+ qr(i,j,k))
        ENDIF
        if (tabs(i,j,k) .eq.0.) print*,'tempstart', i,k,t(i,j,k), (qc(i,j,k)+ qr(i,j,k))
      ENDDO
      ENDDO
      ENDDO
      
      DO k = kts,kte
          pcgs(K)=pres(K)*1000.      ! mb to dyne/cm^2
          rhocgs(K)=rho(K)*0.001      ! kgm-3 to gcm-3
          zcgs(K)=z(k)*100.           !
      ENDDO

!      print*,'after loop', pcgs(1),rhocgs(1),tabs(2,2,1)

! calculate the total number in the domain
!      totnum_beforem = 0.
!      totccn_beforem = 0.0
!      totdrop_beforem = 0.0
!       DO i = its,ite 
!       DO j = jts,jte
!       DO k = kts,kte
!         totnum_beforem = totnum_beforem + qna(i,j,k)+ qnc(i,j,k)+qnr(i,j,k)
!         totccn_beforem = totccn_beforem+ qna(i,j,k)
!         totdrop_beforem = totdrop_beforem + qnc(i,j,k)+qnr(i,j,k)
!       END DO
!       END DO
!       END DO
!       print*, 'domain totnum_beforem', totnum_beforem, totccn_beforem, totdrop_beforem

!
      do j = jts,jte
      do k = kts,kte
      do i = its,ite

! Initalize the ice nucleation rates
      rnfreez(i,j,k) = 0.0
      rnicr(i,j,k) = 0.0
      rnfr_hom(i,j,k) = 1.e-36     ! set a value large than 0 to
      rnfr_imm(i,j,k) = 1.e-36     ! get the log values for output
      rndropr(i,j,k) = 0.0 
! Initalize VARIABLES FOR TENDENCY
      difful_tend(i,j,k) = 0.0
      diffui_tend(i,j,k) = 0.0
      frzl_tend(i,j,k) = 0.0
!
! added for sheba output
     diffui_tend(i,j,k)=0.0
! test
   inbeforenuc=0.0
   inafternuc=0.0
   icebeforenuc=0.0
   iceafternuc=0.0
   iceafterdiff=0.0
!

! LIQUID
        DO KR=1,NKR
          FF1R(KR)=ffcd(I,j,k,KR)*rhocgs(k)/xl(kr)/xl(kr)/3.0
          IF (FF1R(KR).LT.0) FF1R(KR)=0.
!test
!         if (i.eq.20.and.k.eq.kmax)then
!                          print*,'ff1r(kr) above = ',kr,ff1r(kr)
!         end if
!test
        END DO   
! CCN

         DO KR=1,NKR
!mo          FCCN0(KR)=FCCNR_mp(KR)    ! constant concentration with height
          FCCN0(KR)=FCCNR_mp(KR)*rhocgs(k)/rhocgs(1)       ! constant mixing ratio with height
         END DO

        if (itimestep.eq.1.and.icycle.eq.1)then
! print*, diagCCN, iceform
         DO KR=1,NKR
!mo        FCCN0(KR)=FCCNR_mp(KR)    ! constant concentration with height
!          FCCN0(KR)=FCCNR_mp(KR)*rhocgs(k)/rhocgs(1)       ! constant mixing ratio with height
!!!block          if (z(k) < 800.) then
            fncn(I,J,k,KR)=FCCN0(KR)/rhocgs(k)*xccn(kr)
! Add for IN
            ffin(I,J,k,KR)=fracin*FCCN0(KR)/rhocgs(k)*xccn(kr) ! constant mixing ratio with height
!!!block          else
!!!block            fncn(I,J,k,KR)=0.4 * FCCN0(KR)/rhocgs(k)*xccn(kr)
!!!block            ffin(I,j,k,kr)=0.                            ! above boundary layer
!!!block
!!!block          end if

! For Sheba comparison
!           ffin(I,J,k,KR)=1.7e-3/col/33.0/rhocgs(k)*xccn(kr)
!         chem_new(I,j,k,KR)=FCCNR2(KRR)
         end do
        end if
! Add for diagnostic CCN
        if (diagCCN) then
           ndrop=0.0
           DO KR=1,NKR
           ndrop=ndrop+COL*ffcd(I,j,k,KR)/XL(KR)*rhocgs(k)
           ENDDO
           FCCN(:) = FCCN0(:)
           if (ndrop >= FCCN0(NKR)*COL) then
             subtot=0.0
             do kr = nkr, 1, -1
              subtot=subtot+FCCN0(kr)*COL
              FCCN(kr)=0.0
              if (subtot >= ndrop) then
               FCCN(kr)= (subtot-ndrop)/COL
               exit
              endif
             enddo
           endif                            
        else      ! if (diagCCN) not true, i.e., prognostic       
            
         DO KR=1,NKR
          FCCN(KR)=fncn(I,J,k,KR)*rhocgs(k)/xccn(kr)
          if (fccn(kr).lt.0) fccn(kr)=0.
         END DO
       endif      ! if (diagCCN)

! add for IN
         DO KR=1,NKR
          FIN(KR)=ffin(I,J,k,KR)*rhocgs(k)/xccn(kr)
          if (fin(kr).lt.0) fin(kr)=0.
         END DO


! check the balance of the particles
!         print*, 'check ccn', sum(fccnr_mp(:)), sum(fccn(:))
!         totnum_before = 0.0
!         do kr  = 1, nkr
!            totnum_before = totnum_before + fccn(kr)*col+ &
!      &                        ff1r(kr)*3.0*col*xl(kr)
!            fccn_before(kr) = fccn(kr)
!            ff1r_before(kr) = ff1r(kr)
!         enddo
                                                                                              
        IF (ICEPROCS.EQ.1)THEN
! COLUMNS!       
         DO KR=1,NKR
          FF2R(kr,1)=ffic(I,J,k,KR)*rhocgs(k)/xi(kr,1)/xi(kr,1)/3.0
          if (ff2r(kr,1).lt.0) ff2r(kr,1)=0.
         END DO
! PLATES!
         
         DO KR=1,NKR
          FF2R(KR,2)=ffip(I,J,k,KR)*rhocgs(k)/xi(kr,2)/xi(kr,2)/3.0
          if (ff2r(kr,2).lt.0) ff2r(kr,2)=0.
         END DO
! DENDRITES!
         DO KR=1,NKR
          FF2R(KR,3)=ffid(I,J,k,KR)*rhocgs(k)/xi(kr,3)/xi(kr,3)/3.0
          if (ff2r(kr,3).lt.0) ff2r(kr,3)=0.
         END DO
! SNOW
         DO KR=1,NKR
          FF3R(KR)=ffsn(I,J,k,KR)*rhocgs(k)/xs(kr)/xs(kr)/3.0
          if (ff3r(kr).lt.0) ff3r(kr)=0.
         END DO
! Graupel
         DO KR=1,NKR
          FF4R(KR)=ffgl(I,J,k,KR)*rhocgs(k)/xg(kr)/xg(kr)/3.0
          if (ff4r(kr).lt.0) ff4r(kr)=0.
         END DO
! Hail
         DO KR=1,NKR
          FF5R(KR)=ffhl(I,J,k,KR)*rhocgs(k)/xh(kr)/xh(kr)/3.0
          if (ff5r(kr).lt.0) ff5r(kr)=0.
         END DO

! check units for FF1R, FF2R....

!    if (i==20.and.j==1.and.k==kmax) then
!    do kr = 1,nkr
!        print*,'ff1r(kr) after melt ',kr,ff1r(kr)
!    end do
!    endif
!
       IF(K.EQ.KTE)THEN
        DZZ(K)=(zcgs(K)-zcgs(k-1))
       ELSE IF(K.EQ.1)THEN
        DZZ(K)=(zcgs(k+1)-zcgs(K))
       ELSE
        DZZ(K)=(zcgs(k+1)-zcgs(k-1))
       END IF
       ES2N=AA2_MY*EXP(-BB2_MY/t_old(I,J,K))
       EW1N=qv_old(I,J,K)*pcgs(K)/(0.622+0.378*qv_old(I,J,K))
       SUPICE(K)=EW1N/ES2N-1.
       IF(SUPICE(K).GT.0.5) SUPICE(K)=.5

!     if (i.eq.95.and.k.eq.17.and.j.eq.2)print*,' here at a'
!New changes about the calculation of SS - Fan 08/2008
!             IF (I.LT.ITE.AND.J.LT.JTE)THEN
!              UX=25.*(dudt(I,j,k,nc)+dudt(I+1,j,k,nc)+dudt(I,J+1,k,nc)+dudt(I+1,J+1,k,nc))
!              VX=25.*(dvdt(I,j,k,nc)+dvdt(I+1,j,k,nc)+dvdt(I,J+1,k,nc)+dvdt(I+1,J+1,k,nc))
              UX=50.*(dudt(I,j,k,nc)+dudt(I+1,j,k,nc))
              if (RUN3D) then
              VX=50.*(dvdt(I,J,k,nc)+dvdt(I,J+1,k,nc))
              else
              VX=0.
              endif
              WX=50.*(dwdt(I,J,k,nc)+dwdt(I,J,k+1,nc))
!             ELSE
!              UX=U(I,j,k)
!              VX=V(I,j,k)
!              UX=dudt(I,j,k,nc)*100       ! ?Originally, no conversion here, but i think it is obviously wrong.
!                                     The unit of this UX is not consistent with the one above.
!              VX=dvdt(I,j,k,nc)*100
!             END IF
! skip the caculation of  DTFREEZ_XYZ since it is not used
!       GOTO 1111
                                                                                             
       IF(t_old(I,J,k).GE.238.15.AND.t_old(I,J,k).LT.274.15) THEN
             IF(K.EQ.1) DERIVT_Z=(t_old(I,j,k+1)-t_old(I,j,k))/DZZ(K)
             IF(K.EQ.KTE) DERIVT_Z=(t_old(I,j,k)-t_old(I,j,k-1))/DZZ(K)
             IF(K.GT.1.AND.K.LT.KTE) DERIVT_Z= &
     &        (t_old(I,j,K+1)-t_old(I,j,K-1))/DZZ(K)
             IF (I.EQ.1)THEN
              DERIVT_X=(t_old(I+1,J,K)-t_old(I,J,K))/(DXHUCM)
             ELSE IF (I.EQ.ITE)THEN
              DERIVT_X=(t_old(I,J,K)-t_old(I-1,J,K))/(DXHUCM)
             ELSE
              DERIVT_X=(t_old(I+1,J,K)-t_old(I-1,J,K))/(2.*DXHUCM)
             END IF
             IF (J.EQ.1.and.J.EQ.JTE)THEN
              DERIVT_Y=0.0
             else if (J.EQ.1.and.JTE.gt.1) then
              DERIVT_Y=(t_old(I,J+1,k)-t_old(I,J,k))/(DYHUCM)
             ELSE IF (J.EQ.JTE.and.JTE.gt.1)THEN
              DERIVT_Y=(t_old(I,J,k)-t_old(I,J-1,k))/(DYHUCM)
             ELSE
              DERIVT_Y=(t_old(I,J+1,k)-t_old(I,J-1,k))/(2.*DYHUCM)
             END IF
!             DTFREEZ_XYZ(I,j,k)=DT*(VX*DERIVT_Y+ &
!     &            UX*DERIVT_X+100.*dwdt(I,j,k,nc)*DERIVT_Z)
! use WX to replace original code
             DTFREEZ_XYZ(I,j,k)=DT*(VX*DERIVT_Y+ &
     &            UX*DERIVT_X+WX*DERIVT_Z)
          ELSE
             DTFREEZ_XYZ(I,j,k)=0.
          ENDIF
          IF(SUPICE(K).GE.0.02.AND.t_old(I,j,k).LT.268.15) THEN
            IF (I.LT.ITE)THEN
             ES2NPLSX=AA2_MY*EXP(-BB2_MY/t_old(I+1,j,k))
             EW1NPLSX=qv_old(I+1,j,k)*pcgs(k)/ &
     &               (0.622+0.378*qv_old(I+1,j,k))
            ELSE
             ES2NPLSX=AA2_MY*EXP(-BB2_MY/t_old(I,j,k))
             EW1NPLSX=qv_old(I,j,k)*pcgs(k)/ &
     &               (0.622+0.378*qv_old(I,j,k))
            END IF
            IF (ES2NPLSX.EQ.0)THEN
             DEL2INPLSX=0.5
            ELSE
             DEL2INPLSX=EW1NPLSX/ES2NPLSX-1.
            END IF
            IF(DEL2INPLSX.GT.0.5) DEL2INPLSX=.5
            IF (I.GT.1)THEN
             ES2N=AA2_MY*EXP(-BB2_MY/t_old(I-1,j,k))
             EW1N=qv_old(I-1,j,k)*pcgs(k)/(0.622+0.378*qv_old(I-1,j,k))
            ELSE
             ES2N=AA2_MY*EXP(-BB2_MY/t_old(I,j,k))
             EW1N=qv_old(I,j,k)*pcgs(k)/(0.622+0.378*qv_old(I,j,k))
            END IF
            DEL2IN=EW1N/ES2N-1.
            IF(DEL2IN.GT.0.5) DEL2IN=.5
            IF (I.GT.1.AND.I.LT.ITE)THEN
             DERIVS_X=(DEL2INPLSX-DEL2IN)/(2.*DXHUCM)
            ELSE
             DERIVS_X=(DEL2INPLSX-DEL2IN)/(DXHUCM)
            END IF
            IF (J.LT.JTE)THEN
             ES2NPLSY=AA2_MY*EXP(-BB2_MY/t_old(I,J+1,k))
             EW1NPLSY=qv_old(I,J+1,k)*pcgs(k)/(0.622+0.378*qv_old(I,J+1,k))
            ELSE
             ES2NPLSY=AA2_MY*EXP(-BB2_MY/t_old(I,j,k))
             EW1NPLSY=qv_old(I,j,k)*pcgs(k)/(0.622+0.378*qv_old(I,j,k))
            END IF
            DEL2INPLSY=EW1NPLSY/ES2NPLSY-1.
            IF(DEL2INPLSY.GT.0.5) DEL2INPLSY=.5
            IF (J.GT.1)THEN
             ES2N=AA2_MY*EXP(-BB2_MY/t_old(I,J-1,k))
             EW1N=qv_old(I,J-1,k)*pcgs(k)/(0.622+0.378*qv_old(I,J-1,k))
            ELSE
             ES2N=AA2_MY*EXP(-BB2_MY/t_old(I,J,k))
             EW1N=qv_old(I,j,k)*pcgs(k)/(0.622+0.378*qv_old(I,J,k))
            END IF
             DEL2IN=EW1N/ES2N-1.
            IF(DEL2IN.GT.0.5) DEL2IN=.5
            IF (J.GT.1.AND.J.LT.JTE)THEN
             DERIVS_Y=(DEL2INPLSY-DEL2IN)/(2.*DYHUCM)
            ELSE
             DERIVS_Y=(DEL2INPLSY-DEL2IN)/(DYHUCM)
            END IF
!
            IF (K.EQ.1)DERIVS_Z=(SUPICE(K+1)-SUPICE(K))/DZZ(K)
            IF (K.EQ.KTE)DERIVS_Z=(SUPICE(K)-SUPICE(K-1))/DZZ(K)
            IF(K.GT.1.and.K.LT.KTE) DERIVS_Z=(SUPICE(K+1)-SUPICE(K-1))/DZZ(K)

!            IF (I.LT.ITE.AND.J.LT.JTE)THEN
!             UX=25.*(dudt(I,j,k,nc)+dudt(I+1,j,k,nc)+dudt(I,J+1,k,nc)+dudt(I+1,J+1,k,nc))
!             VX=25.*(dvdt(I,J,k,nc)+dvdt(I+1,J,k,nc)+dvdt(I,J+1,k,nc)+dvdt(I+1,J+1,k,nc))
!            ELSE
!!             UX=U(I,j,k)
!!             VX=V(I,j,k)
!             UX=dudt(I,J,k,nc)*100     ! see above note
!             VX=dvdt(I,J,k,nc)*100
!            END IF  
!            DSUPICE_XYZ(I,J,k)=(UX*DERIVS_X+VX*DERIVS_Y+ &
!      &                        100.*dwdt(I,j,k,nc)*DERIVS_Z)*DTCOND
            DSUPICE_XYZ(I,J,k)=(UX*DERIVS_X+VX*DERIVS_Y+ &
      &                        WX*DERIVS_Z)*DT
          ELSE
            DSUPICE_XYZ(I,j,k)=0.0
          END IF
!           DSUPICE_XYZ(I,J,k)=0.0
!           DTFREEZ_XYZ(I,j,k)=0.
!       if(i==20.and.j==1.and.k==20) print*,'DSUPICE(20,1,20) = ',DSUPICE_XYZ(I,j,k),DSUPICE_XYZ(I,j,k)
!        print*,' DTFREEZ_XYZ(I,J,k) = ', DTFREEZ_XYZ(202,1,1),DTFREEZ_XYZ(201,1,1)
     

! check units for FF1R, FF2R....
    totlbf_frz=sum(ff1r(:)/rhocgs(k)*xl(:)*xl(:)*3)

       if (fixice .ne.1) then
          CALL FREEZ &
     &     (FF1R,XL,FF2R,XI,FF3R,XS,FF4R,XG,FF5R,XH, &
     &      tabs(I,J,k),DT,rhocgs(k), &
     &      COL,AFREEZMY,BFREEZMY,BFREEZMAX, &
     &      KRFREEZ,ICEMAX,NKR,rnfr)
         rnfreez(i,j,k) = rnfr
   totlaf_frz=sum(ff1r(:)/rhocgs(k)*xl(:)*xl(:)*3)
   frzl_tend(i,j,k) = (totlaf_frz-totlbf_frz)/dt   !kg/kg/s
       end if                                                                                           
         CALL MELT &
     &    (FF1R,XL,FF2R,XI,FF3R,XS,FF4R,XG,FF5R,XH, &
     &     tabs(I,j,k),DT,rhocgs(k),COL,ICEMAX,NKR)
        ENDIF  ! if iceprocess == 1

        IF (t_old(I,j,k).GT.233)THEN
! TEMPMIN_MICRO changes for different cloud cases
!        IF (t_old(I,j,k).GT.TEMPMIN_MICRO) THEN      
         TT=t_old(I,j,k)
         QQ=QV_OLD(I,j,k)
         PP=pcgs(k)
!         TTA=T_NEW(I,j,k)
         TTA=tabs(i,j,k)
!         QQA=QV(I,j,k)
         QQA=qv(i,j,k)
         ES1N=AA1_MY*DEXP(-BB1_MY/TT)
         ES2N=AA2_MY*DEXP(-BB2_MY/TT)
         EW1N=QQ*PP/(0.622+0.378*QQ)
         DIV1=EW1N/ES1N

         DEL1IN=EW1N/ES1N-1.
         DIV2=EW1N/ES2N
         DEL2IN=EW1N/ES2N-1.
         ES1N=AA1_MY*DEXP(-BB1_MY/TTA)
         ES2N=AA2_MY*DEXP(-BB2_MY/TTA)
         EW1N=QQA*PP/(0.622+0.378*QQA)
         DIV3=EW1N/ES1N
         DEL1AD=EW1N/ES1N-1.
         DIV4=EW1N/ES2N
         DEL2AD=EW1N/ES2N-1.
         SUP2_OLD=DEL2IN
         DELSUP1=(DEL1AD-DEL1IN)/NCOND
         DELSUP2=(DEL2AD-DEL2IN)/NCOND
         DELDIV1=(DIV3-DIV1)/NCOND
         DELDIV2=(DIV4-DIV2)/NCOND
         if (DIV4 .le. 1.e-15) print*,'div2 a', DIV4, DIV2, NCOND, DELDIV2
         DELTATEMP=0
         DELTAQ=0
         tt_old = TT
! added for sheba output (tendency of diffusion growth)
        totibf_diffu =0.0
        totiaf_diffu =0.0
        totlbf_diffu =0.0
        totlaf_diffu =0.0


         DO IKL=1,NCOND
          DEL1IN=DEL1IN+DELSUP1
          DEL2IN=DEL2IN+DELSUP2
          DIV1=DIV1+DELDIV1
          DIV2=DIV2+DELDIV2
         if (DIV2 == 0.) print*,'div2 b', DIV2, DELDIV2,DIV4, NCOND
!959       format (' ',i3,1x,f7.1,1x,f6.1,1x,f6.4,1x,f6.2,1x,f6.3)
          DIFFU=1
          IF (DIV1.GT.DIV2.AND.TT.LE.265)THEN
           print*,'div1 > div2',div1,div2
           print*,'STOP'
           print*,'RESET'
           DIV2=0.99999*DIV1
           DEL2IN=0.99999*DEL2IN
!          STOP
           DIFFU=0
          END IF
          DEL1NR=A1_MYN*(100.*DIV1)
          DEL2NR=A2_MYN*(100.*DIV2)
          if (DEL2NR.EQ.0) print*, 'DEL2NR',tt,qq, tta,qqa
          IF (DEL2NR.EQ.0)PRINT*,'DEL2NR = 0'
          IF (DEL2NR.EQ.0)STOP
          DEL12R=DEL1NR/DEL2NR
          DEL12RD=DEL12R**DEL_BBR
          EW1PN=AA1_MY*100.*DIV1*DEL12RD/100.
          IF (DEL12R.EQ.0)PRINT*,'DEL12R = 0'
          IF (DEL12R.EQ.0)STOP
          TT=-DEL_BB/DLOG(DEL12R)
          QQ=0.622*EW1PN/(PP-0.378*EW1PN)
          DO KR=1,NKR
            FF1IN(KR)=FF1R(KR)
            DO ICE=1,ICEMAX
             FF2IN(KR,ICE)=FF2R(KR,ICE)
            ENDDO
          ENDDO

          IF (BULKNUC.eq.1)THEN
            IF (DEL1IN.GT.0)THEN
              IF (zcgs(k).LE.500.E2)THEN
                FACTZ=0.
              ELSE
                FACTZ=1
!               FACTZ=EXP(-(zcgs(k)-2.E5)/Z0IN)
              END IF
             CONCCCN_XZ=FACTZ*ACCN*(100.*DEL1IN)**BCCN

             CONCDROP=0.D0

             DO KR=1,NKR
               CONCDROP=CONCDROP+FF1IN(KR)*XL(KR)
             ENDDO

             CONCDROP=CONCDROP*3.D0*COL
!            print*,'factz,accn,bccn,del1in,conccnz_xz,concdrop= '
!            print*,factz,accn,bccn,del1in,concccn_xz,concdrop
             IF(CONCCCN_XZ.GT.CONCDROP) &
     &       FF1IN(1)=FF1IN(1)+(CONCCCN_XZ-CONCDROP)/(3.D0*COL*XL(1))
            END IF
          ELSE
!        inbeforenuc = 1.e3*sum(fin)*col
!        icebeforenuc = 1.e3*(sum(FF2IN(:,1)*3*col*xi(:,1))+sum(FF2IN(:,2)*3*col*xi(:,2))+sum(FF2IN(:,3)*3*col*xi(:,3)))
!        if (i > 10 .and. i < 20 ) print*, 'before nucl',i,k,inbeforenuc, &
!     &     icebeforenuc,sum(FCCN),sum(FF1IN),TT, DEL1IN, DEL2IN
            IF(DEL1IN.GT.0.OR.DEL2IN.GT.0)THEN
!mo Limit SS for the droplet nucleation
             del1in_limit = min(DEL1IN,ss_max)
!             CONCDROP=0.D0
!             DO KR=1,NKR
!               CONCDROP=CONCDROP+FF1IN(KR)*XL(KR)
!             ENDDO
!             CONCDROP=CONCDROP*3.D0*COL
!             if (CONCDROP > ) del1in_limit = 0.
!mo
!
!MO: account for condensation of vapor during droplet nucleation and 
!            for deposition of vapor on ice during ice nucleation
        totlbf_diffu = totlbf_diffu+sum(ff1in(:)/rhocgs(k)*xl(:)*xl(:)*3)

        totibf_diffu = totibf_diffu+ sum(ff2in(:,1)/rhocgs(k)*xi(:,1)*xi(:,1)*3)+ &
   &                 sum(ff2in(:,2)/rhocgs(k)*xi(:,2)*xi(:,2)*3)+ &
   &                 sum(ff2in(:,3)/rhocgs(k)*xi(:,3)*xi(:,3)*3)+ &
   &                 sum(ff3in(:)/rhocgs(k)*xs(:)*xs(:)*3)+       &
   &                 sum(ff4in(:)/rhocgs(k)*xg(:)*xg(:)*3) 

             CALL JERNUCL01(FF1IN,FF2IN,FCCN &
     &       ,XL,XI,TT,QQ &
     &       ,rhocgs(k),pcgs(k) &
!     &       ,DEL1IN,DEL2IN &
     &       ,DEL1IN_limit,DEL2IN &
     &       ,COL,AA1_MY, BB1_MY, AA2_MY,BB2_MY &
     &       ,C1_MEY,C2_MEY,SUP2_OLD,DSUPICE_XYZ(I,j,k) &
     &       ,RCCN,DROPRADII,NKR,ICEMAX,DT,ICEPROCS,ICEFLAG,fin,rnic &
     &       ,ff3r,xs,ff4r,xg)
! output ice nucleation rate
             rnicr(i,j,k) = rnic

        totlaf_diffu = totlaf_diffu+sum(ff1in(:)/rhocgs(k)*xl(:)*xl(:)*3)
        totiaf_diffu = totiaf_diffu+sum(ff2in(:,1)/rhocgs(k)*xi(:,1)*xi(:,1)*3)+ & 
   &                 sum(ff2in(:,2)/rhocgs(k)*xi(:,2)*xi(:,2)*3)+ &
   &                 sum(ff2in(:,3)/rhocgs(k)*xi(:,3)*xi(:,3)*3)+ &
   &                 sum(ff3in(:)/rhocgs(k)*xs(:)*xs(:)*3)+       &
   &                 sum(ff4in(:)/rhocgs(k)*xg(:)*xg(:)*3)



!        inafternuc =1.e3* sum(fin)*col
!        iceafternuc = 1.e3*(sum(FF2IN(:,1)*3*col*xi(:,1))+sum(FF2IN(:,2)*3*col*xi(:,2))+sum(FF2IN(:,3)*3*col*xi(:,3)))
!        if (i > 10 .and. i < 20 ) print*, 'diff nucl',i,k, inafternuc, &
!     &     iceafternuc,(inafternuc-inbeforenuc), (iceafternuc-icebeforenuc),sum(FCCN),sum(FF1IN)
       
            END IF
          END IF
          DO KR=1,NKR
            FF1R(KR)=FF1IN(KR)
            DO ICE=1,ICEMAX
             FF2R(KR,ICE)=FF2IN(KR,ICE)
            ENDDO
          ENDDO
          FMAX1=0.
          FMAX2=0.
          FMAX3=0.
          FMAX4=0.
          FMAX5=0.
          DO KR=1,NKR
            FF1IN(KR)=FF1R(KR)
            FMAX1=AMAX1(FF1R(KR),FMAX1)
            FF3IN(KR)=FF3R(KR)
            FMAX3=AMAX1(FF3R(KR),FMAX3)
            FF4IN(KR)=FF4R(KR)
            FMAX4=AMAX1(FF4R(KR),FMAX4)
            FF5IN(KR)=FF5R(KR)
            FMAX5=AMAX1(FF5R(KR),FMAX5)
            DO ICE=1,ICEMAX
             FF2IN(KR,ICE)=FF2R(KR,ICE)
             FMAX2=AMAX1(FF2R(KR,ICE),FMAX2)
            END DO
          END DO
          ISYM1=0
          ISYM2=0
          ISYM3=0
          ISYM4=0
          ISYM5=0
!mo          IF(FMAX1.GT.0)ISYM1=1
          IF(FMAX1.GT.1.e-6)ISYM1=1
          IF (ICEPROCS.EQ.1)THEN
! Reduced threshold from 1.e-4 to 1.e-6 for Arctic stratus cases (MO; 03-Mar-2010)
           IF(FMAX2.GT.1.E-6)ISYM2=1
           IF(FMAX3.GT.1.E-6)ISYM3=1
           IF(FMAX4.GT.1.E-6)ISYM4=1
           IF(FMAX5.GT.1.E-6)ISYM5=1
          END IF

! check the balance of the particles
!         totnum_afternul = 0.0
!         do kr  = 1, nkr
!            totnum_afternul = totnum_afternul + fccn(kr)*col+ &
!      &                        ff1r(kr)*3.0*col*xl(kr)
!         enddo

! added for SHEBA output
        totlbf_diffu = totlbf_diffu+sum(ff1r(:)/rhocgs(k)*xl(:)*xl(:)*3)

        totibf_diffu = totibf_diffu+ sum(ff2r(:,1)/rhocgs(k)*xi(:,1)*xi(:,1)*3)+ &
   &                 sum(ff2r(:,2)/rhocgs(k)*xi(:,2)*xi(:,2)*3)+ &
   &                 sum(ff2r(:,3)/rhocgs(k)*xi(:,3)*xi(:,3)*3)+ &
   &                 sum(ff3r(:)/rhocgs(k)*xs(:)*xs(:)*3)+       &
   &                 sum(ff4r(:)/rhocgs(k)*xg(:)*xg(:)*3) 

!        if (i > 10 .and. i < 20 ) print*,'before diffu',i,k,sum(fccn), sum(fin), sum(ff1in), sum(ff1in)
          ccnreg=0.
          inreg=0.
        
          IF (DIFFU.NE.0)THEN
          IF(ISYM1.EQ.1.AND.((TT-273.15).GT.-0.187.OR. &
     &     (ISYM2.EQ.0.AND. &
     &     ISYM3.EQ.0.AND.ISYM4.EQ.0.AND.ISYM5.EQ.0)))THEN
           CALL ONECOND1(TT,QQ,PP,rhocgs(k) &
     &      ,VR1,pcgs(k) &
     &      ,DEL1IN,DEL2IN,DIV1,DIV2 &
     &      ,FF1R,FF1IN,XL,RLEC,RO1BL &
     &      ,AA1_MY,BB1_MY,AA2_MY,BB2_MY &
     &      ,C1_MEY,C2_MEY &
     &      ,COL,DTCOND,ICEMAX,NKR)
          ELSE IF(ISYM1.EQ.0.AND.(TT-273.15).LE.-0.187.AND. &
     &     (ISYM2.EQ.1.OR.ISYM3.EQ.1.OR.ISYM4.EQ.1.OR.ISYM5.EQ.1))THEN
           CALL ONECOND2(TT,QQ,PP,rhocgs(k) &
     &      ,VR2,VR3,VR4,VR5,pcgs(k) &
     &      ,DEL1IN,DEL2IN,DIV1,DIV2 &
     &      ,FF2R,FF2IN,XI,RIEC,RO2BL &
     &      ,FF3R,FF3IN,XS,RSEC,RO3BL &
     &      ,FF4R,FF4IN,XG,RGEC,RO4BL &
     &      ,FF5R,FF5IN,XH,RHEC,RO5BL &
     &      ,AA1_MY,BB1_MY,AA2_MY,BB2_MY &
     &      ,C1_MEY,C2_MEY &
     &      ,COL,DTCOND,ICEMAX,NKR &
     &      ,ISYM2,ISYM3,ISYM4,ISYM5)
          ELSE IF(ISYM1.EQ.1.AND.(TT-273.15).LE.-0.187.AND. &
     &     (ISYM2.EQ.1.OR.ISYM3.EQ.1.OR.ISYM4.EQ.1 &
     &     .OR.ISYM5.EQ.1))THEN
           CALL ONECOND3(TT,QQ,PP,rhocgs(k) &
     &      ,VR1,VR2,VR3,VR4,VR5,pcgs(k) &
     &      ,DEL1IN,DEL2IN,DIV1,DIV2 &
     &      ,FF1R,FF1IN,XL,RLEC,RO1BL &
     &      ,FF2R,FF2IN,XI,RIEC,RO2BL &
     &      ,FF3R,FF3IN,XS,RSEC,RO3BL &
     &      ,FF4R,FF4IN,XG,RGEC,RO4BL &
     &      ,FF5R,FF5IN,XH,RHEC,RO5BL &
     &      ,AA1_MY,BB1_MY,AA2_MY,BB2_MY &
     &      ,C1_MEY,C2_MEY &
     &      ,COL,DTCOND,ICEMAX,NKR &
     &      ,ISYM1,ISYM2,ISYM3,ISYM4,ISYM5)
          END IF

          END IF
! added for SHEBA output
        totlaf_diffu = totlaf_diffu+sum(ff1r(:)/rhocgs(k)*xl(:)*xl(:)*3)
        totiaf_diffu = totiaf_diffu+sum(ff2r(:,1)/rhocgs(k)*xi(:,1)*xi(:,1)*3)+ & 
   &                 sum(ff2r(:,2)/rhocgs(k)*xi(:,2)*xi(:,2)*3)+ &
   &                 sum(ff2r(:,3)/rhocgs(k)*xi(:,3)*xi(:,3)*3)+ &
   &                 sum(ff3r(:)/rhocgs(k)*xs(:)*xs(:)*3)+       &
   &                 sum(ff4r(:)/rhocgs(k)*xg(:)*xg(:)*3)

          DIFFU=1

!        iceafterdiff = 1.e3*sum(FF2r(:,1)*3*col*xi(:,1))+sum(FF2r(:,2)*3*col*xi(:,2))+sum(FF2r(:,3)*3*col*xi(:,3))
!        if (i > 10 .and. i < 20 .and. inreg > 0.) print*, 'after diff',i,k,iceafterdiff,(iceafterdiff-iceafternuc), inreg*1.e3, ccnreg, sum(FF1r)

!!! For ccn regeneration from evaporation (J. Fan Oct 2007)
!         if (inreg*1.e3 > 5.e-5) print*, 'inreg', inreg*1.e3, ccnreg 
      if (fixice .ne.1 .and. iceflag == 1) then 
       if (iceform == 1) then   !!!! For drop freezing from evaporation (J. Fan Oct 2007)
         if (ccnreg > 0.0 .and. TT < (273.15-5.0)) then
           frzfract = 0.8e-5
           call evapfrz(ccnreg,frzfract,NKR,ICEMAX, TT, DT, xi,ff2r,rndrop)
             rndropr(i,j,k) = rndrop
         endif
         if (inreg > 0.) then
           fin(nkr)=fin(nkr)+inreg/col
         endif
       else if (iceform == 2) then     !!!! part of drop evaporating residuals back to inreg
!         if (itimestep > 3600 ) then
          inreg= inreg + ccnreg*0.5e-5
!         endif
!         if (inreg > 0.) then
           fin(nkr)=fin(nkr)+inreg/col
!         endif
       endif
     endif   !if (.not.fixice .and. iceflag == 1)

! Put ccnreg back to aerosol (CCN regeneration) if diagCCN = false
        if (diagCCN .eq. .false.) then
          tot_reg  = ccnreg
          if (tot_reg > 0.0) then
          fccnin(:) = fccn(:)  
          call ccn_reg(fccn0,fccnin,fccn,nkr,tot_reg,ff1r,xl,  &
     &     ff2r,xi,ff3r,xs,ff4r,xg)
!MO begin diagnose CCN regeneration
!         if (tot_reg > 10. .AND. masterproc) then
!           write (*,*) '  Tot_reg=', tot_reg
!             write (*,*) '  k   rcn   fccn0    fccn_in     fccn_out'
!           do KR=1,NKR 
!             write (*,*)  kr,rccn(kr),fccn0(kr),fccnin(kr),fccn(kr)
!           end do   !  KR=1,NKR
!         end if
!MO end diagnose CCN regeneration
          endif
         endif

         END DO    ! end NCOND

        diffui_tend(i,j,k) = (totiaf_diffu-totibf_diffu)/DT*col ! kg/kg/s
        difful_tend(i,j,k) = (totlaf_diffu-totlbf_diffu)/DT*col ! kg/kg/s

! collision is called every 3*dt
!mo          IF(mod(itimestep,3).eq.0) THEN
          IF(docoag .and. mod(itimestep,3).eq.0) THEN
             CALL COAL_BOTT_NEW(FF1R,XL,FF2R,XI,FF3R,XS, &
     &       FF4R,XG,FF5R,Xh,TT,QQ,PP,rhocgs(k),dthalf,TCRIT,TTCOAL)
          END IF

!f         th_phy(i,j,k) = tt/pi_phy(i,j,k)
         t_old(i,j,k)= tt
         qv_old(i,j,k)=qq

         tabs(i,j,k) = tt
!         if (i==20.and.j==1.and.k==1) print*, 'temp after', tt, qq, pp
!        print*,'tt_dif = ', i,j,k,tt,tt_old,tt-tt_old
!        if (abs(tt-tt_old).gt.difmax)difmax=tt-tt_old
!        print*,'tt,th_phys = ',i,k,j,tt,th_phy(i,k,j),qq
         qv(i,j,k)=qq
! check
         ssatw(i,j,k)=DEL1IN*1.e2  ! %
         ssati(i,j,k)=DEL2IN*1.e2
        END IF
!      if(i >10 .and. i< 20) print*, 'after micro', i, k, sum(fin)
! CCN
        DO KR=1,NKR
          fncn(I,j,k,KR)=FCCN(KR)/rhocgs(k)*xccn(kr)
        END DO
! Add for IN 
        DO KR=1,NKR
          ffin(I,j,k,KR)=fin(KR)/rhocgs(k)*xccn(kr)
        END DO
! LIQIUD
        DO KR=1,NKR
          ffcd(I,j,k,KR)=FF1R(KR)/rhocgs(k)*xl(kr)*xl(kr)*3.0
        END DO 
  
        IF (ICEPROCS.EQ.1)THEN
! COLUMNS!
         DO KR=1,NKR
          ffic(I,j,k,KR)=FF2R(KR,1)/rhocgs(k)*xi(kr,1)*xi(kr,1)*3
         END DO
! PLATES!
         DO KR=1,NKR
          ffip(I,j,k,KR)=FF2R(KR,2)/rhocgs(k)*xi(kr,2)*xi(kr,2)*3
         END DO
! DENDRITES!
         DO KR=1,NKR 
          ffid(I,j,k,KR)=FF2R(KR,3)/rhocgs(k)*xi(kr,3)*xi(kr,3)*3
         END DO
! SNOW
         DO KR=1,NKR
          ffsn(I,j,k,KR)=FF3R(KR)/rhocgs(k)*xs(kr)*xs(kr)*3
         END DO
! Graupel
         DO KR=1,NKR
          ffgl(I,j,k,KR)=FF4R(KR)/rhocgs(k)*xg(kr)*xg(kr)*3
         END DO
! Hail
         DO KR=1,NKR
          ffhl(I,j,k,KR)=FF5R(KR)/rhocgs(k)*xh(kr)*xh(kr)*3
         END DO
        END IF
!    if (rnfreez(i,j,k) .lt. 0.0) print*,'neg rate', i, k, itimestep, rnfreez

      END DO
      END DO
      END DO

!    print*, 'ic rate', MAXVAL(rnfreez(:,:,:)), MAXVAL(rnicr(:,:,:))
!     if(mod(itimestep,150).eq.0) then
!       call write_nucleation
!     endif 

      gmax=0
      qmax=0
      imax=0
      kmax=0
      qnmax=0
      inmax=0
      knmax=0
      DO j = jts,jte
      DO k = kts,kte
      DO i = its,ite
      QC(I,j,k)=0
      QR(I,j,k)=0
      QI(I,j,k)=0
      QIC(I,j,k)=0
      QIP(I,j,k)=0
      QID(I,j,k)=0
      QS(I,j,k)=0
      QG(I,j,k)=0
      QH(I,j,k)=0
      QNC(I,j,k)=0
      QNR(I,j,k)=0
      QNI(I,j,k)=0
      QNIC(I,j,k)=0
      QNIP(I,j,k)=0
      QNID(I,j,k)=0
      QNS(I,j,k)=0
      QNG(I,j,k)=0
      QNH(I,j,k)=0
      QNA(I,j,k)=0
! initialize effcs - Fandec11 for v6.9.4
      reffc(I,j,k)=0.
      reffi(I,j,k)=0.
! add for IN
      QNIN(I,j,k)=0
! added for SHEBA -mass weighted ice particel fall velocity
      vfice_mw(I,j,k)=0
      vfic(I,j,k)=0
      vfip(I,j,k)=0
      vfid(I,j,k)=0
      vfs(I,j,k)=0
      vfg(I,j,k)=0
! added for qni100
      QNIC100(I,j,k)=0
      QNIP100(I,j,k)=0
      QNID100(I,j,k)=0
      QNS100(I,j,k)=0
      QNG100(I,j,k)=0
      QNH100(I,j,k)=0
      qni100(I,j,k)=0

      top = 0.0
      bottom=0.0
      DO  KR=1,NKR
        IF (KR.LT.KRDROP)THEN
          QC(I,j,k)=QC(I,j,k) &
     &      +ffcd(I,j,k,KR)*COL
          QNC(I,j,k)=QNC(I,j,k) &
     &      +COL*ffcd(I,j,k,KR)/XL(KR)*rhocgs(k)
! calculate effcs - Fandec11 for v6.9.4
        top = top+ffcd(I,j,k,KR)*rhocgs(k)/XL(KR)*DROPRADII(KR)**3
        bottom = bottom+ffcd(I,j,k,KR)*rhocgs(k)/XL(KR)*DROPRADII(KR)**2

        ELSE
          QR(I,j,k)=QR(I,j,k) &
     &      +COL*ffcd(I,j,k,KR)
          QNR(I,j,k)=QNR(I,j,k) &
     &      +COL*ffcd(I,j,k,KR)/XL(KR)*rhocgs(k)
        END IF
      END DO
! calculate effcs - Fandec11 for v6.9.4
!        if (bottom > 0.) then
        if (QC(I,j,k) > 1.e-6 .and. bottom > 0.) then
        reffc(I,j,k) = top/bottom
        endif
! transform from cm to um
        reffc(i,j,k)  = reffc(i,j,k)*1.e4
!        EFFCS(I,K,J)     = MIN(EFFCS(I,K,J),59.9)
!        EFFCS(I,K,J)     = MAX(EFFCS(I,K,J),1.51)
! end - Fandec11

     IF (QC(I,j,k).LT.1.E-8)QC(I,J,k)=0.0
     IF (QR(I,j,k).LT.1.E-8)QR(I,j,k)=0.0
      if (qc(i,j,k).gt.qmax)imax=i
      if (qc(i,j,k).gt.qmax)kmax=k
      if (qc(i,j,k).gt.qmax)qmax=qc(i,j,k)
      if (qnc(i,j,k).gt.qnmax)inmax=i
      if (qnc(i,j,k).gt.qnmax)knmax=k
      if (qnc(i,j,k).gt.qnmax)qnmax=qnc(i,j,k)
   
      IF (ICEPROCS.EQ.1)THEN
       DO  KR=1,33
        QIC(I,j,k)=QIC(I,j,k) &
     &   +COL*ffic(I,j,k,KR)
        QNIC(I,j,k)=QNIC(I,j,k) &
     &   +COL*ffic(I,j,k,KR)/Xi(KR,1)*rhocgs(k)
       END DO
       DO  KR=1,33
        QIP(I,j,k)=QIP(I,j,k) &
     &   +COL*ffip(I,j,k,KR)
        QNIP(I,j,k)=QNIP(I,j,k) &
     &   +COL*ffip(I,j,k,KR)/Xi(KR,2)*rhocgs(k)
       END DO
       DO  KR=1,33
        QID(I,j,k)=QID(I,j,k) &
     &   +COL*ffid(I,j,k,KR)
        QNID(I,j,k)=QNID(I,j,k) &
     &   +COL*ffid(I,j,k,KR)/Xi(KR,3)*rhocgs(k)
       END DO
       QI(I,j,k) = QID(I,j,k)+QIP(I,j,k)+QIC(I,j,k)
       QNI(I,j,k) = QNID(I,j,k)+QNIP(I,j,k)+QNIC(I,j,k)
! added for sheba vfice_mw
       DO  KR=1,33
        vfic(I,j,k)=vfic(I,j,k) &
     &   +COL*ffic(I,j,k,KR)*vr2(kr,1)*0.01  ! cm/s to m/s
       END DO
       DO  KR=1,33
        vfip(I,j,k)=vfip(I,j,k) &
     &   +COL*ffip(I,j,k,KR)*vr2(kr,2)*0.01  ! cm/s to m/s
       END DO
       DO  KR=1,33
        vfid(I,j,k)=vfid(I,j,k) &
     &   +COL*ffid(I,j,k,KR)*vr2(kr,3)*0.01  ! cm/s to m/s
       END DO
     
       DO  KR=1,33 
        QS(I,j,k)=QS(I,j,k) &
     &   +COL*ffsn(I,j,k,KR)
        QNS(I,j,k)=QNS(I,j,k) &
     &   +COL*ffsn(I,j,k,KR)/xs(kr)*rhocgs(k)
       END DO
       DO  KR=1,33
        QG(I,j,k)=QG(I,j,k) &
     &   +COL*ffgl(I,j,k,KR)
        QNG(I,j,k)=QNG(I,j,k) &
     &   +COL*ffgl(I,j,k,KR)/xg(kr)*rhocgs(k)
       END DO
      if (qg(i,j,k).gt.gmax)imax=i
      if (qg(i,j,k).gt.gmax)kmax=k
      if (qg(i,j,k).gt.gmax)gmax=qg(i,j,k)
       DO  KR=1,33
        QH(I,j,k)=QH(I,j,k) &
     &   +COL*ffhl(I,j,k,KR)
        QNH(I,j,k)=QNH(I,j,k) &
     &   +COL*ffhl(I,j,k,KR)/xh(kr)*rhocgs(k)
       END DO
! qni100
       DO  KR=1,33
        if (RADXXO(KR,2) > 0.5e-2) then
         QNIC100(I,j,k)=QNIC100(I,j,k) &
     &   +COL*ffic(I,j,k,KR)/Xi(KR,1)*rhocgs(k)
        endif
        if (RADXXO(KR,3) > 0.5e-2) then
         QNIP100(I,j,k)=QNIP100(I,j,k) &
     &   +COL*ffip(I,j,k,KR)/Xi(KR,2)*rhocgs(k)
        endif
        if (RADXXO(KR,4) > 0.5e-2) then
         QNID100(I,j,k)=QNID100(I,j,k) &
     &   +COL*ffid(I,j,k,KR)/Xi(KR,3)*rhocgs(k)
        endif
        if (RADXXO(KR,5) > 0.5e-2) then
        QNS100(I,j,k)=QNS100(I,j,k) &
     &   +COL*ffsn(I,j,k,KR)/xs(kr)*rhocgs(k)
        endif
        if (RADXXO(KR,6) > 0.5e-2) then
        QNG100(I,j,k)=QNG100(I,j,k) &
     &   +COL*ffgl(I,j,k,KR)/xg(kr)*rhocgs(k)
        endif
        if (RADXXO(KR,7) > 0.5e-2) then
        QNH100(I,j,k)=QNH100(I,j,k) &
     &   +COL*ffhl(I,j,k,KR)/xh(kr)*rhocgs(k)
       endif
      
       END DO
      qni100(I,j,k)= QNIC100(I,j,k)+QNIP100(I,j,k)+QNID100(I,j,k)+ &
     &               QNS100(I,j,k)+QNG100(I,j,k)+ QNH100(I,j,k)     ! cm-3

! added for sheba vfice_mw
       DO  KR=1,33
        vfs(I,j,k)=vfs(I,j,k) &
     &   +COL*ffsn(I,j,k,KR)*vr3(kr)*0.01  ! cm/s to m/s
       END DO
       DO  KR=1,33
        vfg(I,j,k)=vfg(I,j,k) &
     &   +COL*ffgl(I,j,k,KR)*vr4(kr)*0.01  ! cm/s to m/s
       END DO

     if ((qi(I,j,k)+qs(I,j,k)+qg(I,j,k)) > 0) then
     vfice_mw(I,j,k) = (vfic(I,j,k)+vfip(I,j,k)+vfid(I,j,k)+ &
     &     vfs(I,j,k)+vfg(I,j,k))/(qi(I,j,k)+qs(I,j,k)+qg(I,j,k))  ! mean over 5 species
     endif

      ENDIF
   
       DO  KR=1,33
        QNA(I,j,k)=QNA(I,j,k) &
     &   +COL*fncn(I,j,k,KR)/xccn(kr)*rhocgs(k) 
       END DO

! Add for IN
       DO  KR=1,33
        QNIN(I,j,k)=QNIN(I,j,k) &
     &   +COL*ffin(I,j,k,KR)/xccn(kr)*rhocgs(k)
       END DO

! add artificil sources for CCN
!       hydro =QC(I,j,k)+QR(I,j,k)+QI(I,j,k)+QS(I,j,k)+QG(I,j,k)+QH(I,j,k)
!!       hydro =QC(I,j,k)+QR(I,j,k)+QS(I,j,k)
!       DO KR=1,NKR
!       if (hydro .LE. 1.e-8 ) then
!         fncn(I,J,k,KR)=FCCNR_mp(KR)/rhocgs(k)*xccn(kr)
!       else
!         fncn(I,J,k,KR)=0.5*FCCNR_mp(KR)/rhocgs(k)*xccn(kr)
!         FCCN(KR)= FCCNR_mp(kr)
!       end if
!       END DO

      END DO
      END DO
      END DO
! print*, 'qni100', maxval(qni100)
! calculate the total number in the domain
!      totnum = 0.
!      totccn = 0.0
!      totdrop = 0.0
!       DO i = its,ite
!       DO j = jts,jte
!       DO k = kts,kte
!         totnum = totnum + qna(i,j,k)+ qnc(i,j,k)+qnr(i,j,k)
!         totccn = totccn+ qna(i,j,k)
!         totdrop = totdrop + qnc(i,j,k)+qnr(i,j,k)
!       END DO
!       END DO
!       END DO
!       print*, 'domain totnum', totnum, totccn, totdrop

! Since SAM calculates precipitation somewhere else, I just comment out the rain calculation
!      DO j = jts,jte
!      DO i = its,ite
!       DO KR=1,NKR
!        DELTAW=VR1(KR)
!        RAINNC(I,J)=RAINNC(I,J) &
!     &  +10*(rhocgs(1)/RO1BL(KR))*COL*DT*DELTAW* &
!     &           ffcd(I,j,1,KR)
!       END DO
                                                                                                    
!       DO KR=1,NKR
!        DELTAW=VR2(KR,1)
!        RAINNC(I,J)=RAINNC(I,J) &
!     &  +10*(rhocgs(1)/RO1BL(KR))*COL*DT*DELTAW* &
!     &           ffic(I,j,1,KR)
!       END DO
                                                                                                    
!       DO KR=1,NKR
                                                                                                    
!        DELTAW=VR2(KR,2)
!        RAINNC(I,J)=RAINNC(I,J) &
!     &  +10*(rhocgs(1)/RO1BL(KR))*COL*DT*DELTAW* &
!     &           ffip(I,J,1,KR)
!       END DO
                                                                                                    
!       DO KR=1,NKR
                                                                                                    
!        DELTAW=VR2(KR,3)
!        RAINNC(I,J)=RAINNC(I,J) &
!     &  +10*(rhocgs(1)/RO1BL(KR))*COL*DT*DELTAW* &
!     &           ffid(I,J,1,KR)
!       END DO
                                                                                                    
!       DO KR=1, NKR
                                                                                                    
!        DELTAW=VR3(KR)
!        RAINNC(I,J)=RAINNC(I,J) &
!     &  +10*(rhocgs(1)/RO1BL(KR))*COL*DT*DELTAW* &
!     &           ffsn(I,J,1,KR)
!       END DO
!       DO KR=1,NKR
                                                                                                    
!        DELTAW=VR4(KR)
!        RAINNC(I,J)=RAINNC(I,J) &
!     &  +10*(rhocgs(1)/RO1BL(KR))*COL*DT*DELTAW* &
!     &           ffgl(I,J,1,KR)
!       END DO
                                                                                                    
!       DO KR=1,NKR
                                                                                                    
!        DELTAW=VR5(KR)
!        RAINNC(I,J)=RAINNC(I,J) &
!     &  +10*(rhocgs(1)/RO1BL(KR))*COL*DT*DELTAW* &
!     &           ffhl(I,J,1,KR)
!       END DO
!      END DO
!      END DO

      if (docloud)  call micro_diagnose()   ! leave this line here
     
      RETURN
  END SUBROUTINE micro_proc

!---------------------------------------------------------------------------------
      SUBROUTINE micro_init()

      use grid, only: nx, ny, nzm, dt,case   
      use vars
      use params   ! change for 6.9.4 
  
      IMPLICIT NONE
      INTEGER IKERN_0,IKERN_Z,L0_REAL,L0_INTEGER,INEWMEY,INEST
      INTEGER I,J,K,KR
      REAL rho_phys
      INTEGER,parameter :: hujisbm_unit1 = 22
      LOGICAL, PARAMETER :: PRINT_diag=.FALSE.
      LOGICAL :: opened 
      CHARACTER*80 errmess
!      REAL PI     ! block for v6.9.4
      double precision ax
!      data pi/3.141592654/  ! block for v6.9.4

! dtime - timestep of integration (calculated in main program) :
! ax - coefficient used for masses calculation 
! ima(i,j) - k-category number, c(i,j) - courant number 

        REAL C1(NKR,NKR)
!add T_old
!        REAL t_old(nx,ny,nzm), qv_old(nx,ny,nzm)
! DON'T NEED ALL THESE VARIABLES: STILL NEED EDITING
       INTEGER ICE,KGRAN,IPRINT01
       REAL TWSIN,TWCIN,TWNUC,XF5,XF4,XF3,CONCHIN,CONCGIN,CONCSIN, &
     & CONCCLIN,TWHIN,RADH,RADS,RADG,RADL,CONCLIN,A1_MY,A2,A2_MY,XLK, &
     & A1N,A3_MY,A3,A1_MYN,R0CCN,X0DROP,DEG01,CONTCCNIN,CONCCCNIN, &
     & A,B,X0CCN,S_KR,RCCNKR,R0,X0,TWCALLIN,A1,RCCNKR_CM,SUMIIN,TWGIN, &
     & XF1N,XF1,WC1N,RF1N,WNUC,RNUC,WC5,RF5, &
     & WC4,RF4,WC3,RF3,WC1,RF1,SMAX
       REAL TWIIN(ICEMAX)

       real  rccn1(nkr)
!       REAL RO_SOLUTE      
!       PARAMETER (RO_SOLUTE=2.16)


       if (masterproc) then
        PRINT*, 'INITIALIZING HUCM'  
	print *, ' ****** HUCM *******'
       end if

! INPUT :
        dlnr=dlog(2.d0)/(3.d0*scal)
!
       if (masterproc) print*, trim(case)

!--- Read in various lookup tables
!
          OPEN(UNIT=hujisbm_unit1,FILE="./sbm_input/capacity.asc",  &
!          OPEN(UNIT=hujisbm_unit1,FILE="./ARM970621_sbm/sbm_input/capacity.asc",  &
     &        FORM="FORMATTED",STATUS="OLD")

  900	FORMAT(6E13.5)
	READ(hujisbm_unit1,900) RLEC,RIEC,RSEC,RGEC,RHEC
	CLOSE(hujisbm_unit1)
        if (masterproc) print*,'here at 2'
! MASSES :
!
! read ice nuclei size distribution in m^-3
          OPEN(UNIT=hujisbm_unit1,FILE="./sbm_input/ice_nuclei.asc", &
!          OPEN(UNIT=hujisbm_unit1,FILE="./ARM970621_sbm/sbm_input/ice_nuclei.asc", &      
     &        FORM="FORMATTED",STATUS="OLD")
          READ(hujisbm_unit1,*) ICEN
          CLOSE(hujisbm_unit1)
          if (masterproc) print *, ' ***** ice_nuclei file: succesful *******'

!          OPEN(UNIT=hujisbm_unit1,FILE="./'//trim(case)//'/sbm_input/masses.asc", &
           OPEN(UNIT=hujisbm_unit1,FILE="./sbm_input/masses.asc", &
     &        FORM="FORMATTED",STATUS="OLD")
	READ(hujisbm_unit1,900) XL,XI,XS,XG,XH          
	CLOSE(hujisbm_unit1)
	if (masterproc) print *, ' ***** file2: succesfull *******'
! TERMINAL VELOSITY :
!
!          OPEN(UNIT=hujisbm_unit1,FILE="./'//trim(case)//'/sbm_input/termvels.asc",  &
           OPEN(UNIT=hujisbm_unit1,FILE="./sbm_input/termvels.asc", &
     &        FORM="FORMATTED",STATUS="OLD")
	READ(hujisbm_unit1,900) VR1,VR2,VR3,VR4,VR5     
	CLOSE(hujisbm_unit1)
	if (masterproc) print *, ' ***** file3: succesfull *******'
! CONSTANTS :
!          OPEN(UNIT=hujisbm_unit1,FILE="./'//trim(case)//'/sbm_input/constants.asc", &
          OPEN(UNIT=hujisbm_unit1,FILE="./sbm_input/constants.asc", &
     &        FORM="FORMATTED",STATUS="OLD")
	READ(hujisbm_unit1,900) SLIC,TLIC,COEFIN,C2,C3,C4
	CLOSE(hujisbm_unit1)
	if (masterproc) print *, ' ***** file4: succesfull *******'
! CONSTANTS :
! KERNELS DEPENDING ON PRESSURE :
!
!          OPEN(UNIT=hujisbm_unit1,FILE="./'//trim(case)//'/sbm_input/kernels_z.asc",  &
          OPEN(UNIT=hujisbm_unit1,FILE="./sbm_input/kernels_z.asc",  &
     &        FORM="FORMATTED",STATUS="OLD")
        READ(hujisbm_unit1,900)  &
     &  YWLL_1000MB,YWLL_750MB,YWLL_500MB
	CLOSE(hujisbm_unit1)
!
!          OPEN(UNIT=hujisbm_unit1,FILE="./'//trim(case)//'/sbm_input/kernels.asc_s_0_03_0_9",  &
          OPEN(UNIT=hujisbm_unit1,FILE="./sbm_input/kernels.asc_s_0_03_0_9",  &
     &        FORM="FORMATTED",STATUS="OLD")
! KERNELS NOT DEPENDING ON PRESSURE :
	READ(hujisbm_unit1,900) &
     &  YWLL,YWLI,YWLS,YWLG,YWLH, &
     &  YWIL,YWII,YWIS,YWIG,YWIH, &
     &  YWSL,YWSI,YWSS,YWSG,YWSH, &
     &  YWGL,YWGI,YWGS,YWGG,YWGH, &
     &  YWHL,YWHI,YWHS,YWHG,YWHH
       close (hujisbm_unit1)
! BULKDENSITY :
!          OPEN(UNIT=hujisbm_unit1,FILE="./'//trim(case)//'/sbm_input/bulkdens.asc_s_0_03_0_9", & 
          OPEN(UNIT=hujisbm_unit1,FILE="./sbm_input/bulkdens.asc_s_0_03_0_9", &
     &        FORM="FORMATTED",STATUS="OLD")
	READ(hujisbm_unit1,900) RO1BL,RO2BL,RO3BL,RO4BL,RO5BL
	CLOSE(hujisbm_unit1)
	if (masterproc) print *, ' ***** file6: succesfull *******'
! BULKRADIUS
!
!          OPEN(UNIT=hujisbm_unit1,FILE="./'//trim(case)//'/sbm_input/bulkradii.asc_s_0_03_0_9", & 
          OPEN(UNIT=hujisbm_unit1,FILE="./sbm_input/bulkradii.asc_s_0_03_0_9", &
     &        FORM="FORMATTED",STATUS="OLD")
	READ(hujisbm_unit1,*) RADXXO
	CLOSE(hujisbm_unit1)
	if (masterproc) print *, ' ***** file7: succesfull *******'
	if (masterproc) PRINT *, '******* Hebrew Univ Cloud model-HUCM *******'

! calculation of the mass(in mg) for categories boundaries :
        ax=2.d0**(1.0/scal)
        xl_mg(1)=0.3351d-7
	do i=2,nkr
           xl_mg(i)=ax*xl_mg(i-1)
!         if (i.eq.22)print*,'printing xl_mg = ',xl_mg(22)
        enddo
	do i=1,nkr
           xs_mg(i)=xs(i)*1.e3
           xg_mg(i)=xg(i)*1.e3
           xh_mg(i)=xh(i)*1.e3
           xi1_mg(i)=xi(i,1)*1.e3
           xi2_mg(i)=xi(i,2)*1.e3
           xi3_mg(i)=xi(i,3)*1.e3
        enddo
! calculation of c(i,j) and ima(i,j) :
! ima(i,j) - k-category number, c(i,j) - courant number 
        if (masterproc) print*, 'calling courant_bott'
        call courant_bott
        if (masterproc) print*, 'called courant_bott'
 

	DEG01=1./3.

!------------------------------------------------------------------

!        print*,'XL(ICCN) = ',ICCN,XL
	X0DROP=XL(ICCN)
!        print*,'X0DROP = ',X0DROP
	X0CCN =X0DROP/(2.**(NKR-1))
	R0CCN =(3.*X0CCN/4./3.141593/ROCCN0)**DEG01
!------------------------------------------------------------------
! THIS TEXT FROM TWOINITM.F_203
!------------------------------------------------------------------
! TEMPERATURA IN SURFACE LAYER EQUAL 15 Celsius(288.15 K)  
        A=3.3E-05/tsfc_aero
        B=ions*4.3/mwaero
        B=B*(4./3.)*3.14*RO_SOLUTE
        A1=2.*(A/3.)**1.5/SQRT(B)
        A2=A1*100.
!------------------------------------------------------------------
	CONCCCNIN=0.
	CONTCCNIN=0.
	DO KR=1,NKR
           DROPRADII(KR)=(3.*XL(KR)/4./3.141593/1.)**DEG01
        ENDDO
	DO KR=1,NKR
!           print*,'ROCCN0 = ',ROCCN0
!           print*, 'X0CCN = ',X0CCN 
!           print*, 'DEG01 = ',DEG01
	   ROCCN(KR)=ROCCN0
	   X0=X0CCN*2.**(KR-1)
	   R0=(3.*X0/4./3.141593/ROCCN(KR))**DEG01
	   XCCN(KR)=X0
	   RCCN(KR)=R0
           RCCNKR_CM=R0
! CCN SPECTRUM 

! UNCOMMENT END IF  FOR CONTINENTAL/COMMENT FOR MARITIME
!          IF(R0.LE.(0.4E-4)) THEN

           S_KR=A2/RCCNKR_CM**1.5
!           print*,'accn, bccn,S_KR = ',accn,bccn,S_KR
           FCCNR0(KR)=1.5*ACCN*BCCN*S_KR**BCCN
!           print*,'fccnr0(kr) = ',fccnr0(kr)
!          read(6,*)

! for maritime clouds in case of lack of small aerosol(ccn)
! particles we introduce maximum of supersaturation(SMAX) that
! if S > SMAX, % then FCCNR(KR)=0.

! UNCOMENT THESE TWO LINES FOR MARITIME/COMMENT FOR CONTINENTAL
           SMAX=1.1
!          IF(S_KR.GE.SMAX) FCCNR(KR)=0.

	     CONTCCNIN=CONTCCNIN+COL*FCCNR0(KR)*R0*R0*R0
             CONCCCNIN=CONCCCNIN+COL*FCCNR0(KR)
          
! UNCOMMENT FOR CONTINENTAL/COMMENT FOR MARITIME
!         ELSE
!           FCCNR(KR)=0.
!         END IF
	ENDDO
        if (masterproc) PRINT*,    'RCCN(KR)*1.E4,FCCNR0(KR),KR=1,NKR'
        if (masterproc) PRINT 101, (RCCN(KR)*1.E4,FCCNR0(KR),KR=1,NKR)
	if (masterproc) PRINT *, '********* MAR CCN CONCENTRATION & MASS *******'
	if (masterproc) PRINT 200, CONCCCNIN,CONTCCNIN
! CONTINENTAL
	DO KR=1,NKR
	   ROCCN(KR)=ROCCN02
	   X0=X0CCN*2.**(KR-1)
	   R0=(3.*X0/4./3.141593/ROCCN(KR))**DEG01
	   XCCN(KR)=X0
	   RCCN(KR)=R0
           RCCNKR_CM=R0
! CCN SPECTRUM 

! UNCOMMENT END IF  FOR CONTINENTAL/COMMENT FOR MARITIME
           IF(R0.LE.(0.4E-4)) THEN

           S_KR=A2/RCCNKR_CM**1.5
           FCCNR2(KR)=1.5*ACCN2*BCCN2*S_KR**BCCN2

! for maritime clouds in case of lack of small aerosol(ccn)
! particles we introduce maximum of supersaturation(SMAX) that
! if S > SMAX, % then FCCNR(KR)=0.

! UNCOMMENT THESE TWO LINES FOR MARITIME/COMMENT FOR CONTINENTAL
!          SMAX=1.1
!          IF(S_KR.GE.SMAX) FCCNR(KR)=0.

	     CONTCCNIN=CONTCCNIN+COL*FCCNR2(KR)*R0*R0*R0
             CONCCCNIN=CONCCCNIN+COL*FCCNR2(KR)
          
! UNCOMMENT FOR CONTINENTAL/COMMENT FOR MARITIME
          ELSE
            FCCNR2(KR)=0.
          END IF
        END DO
        if (masterproc) PRINT*,    'RCCN(KR)*1.E4,FCCNR2(KR),KR=1,NKR'
        if (masterproc) PRINT 101, (RCCN(KR)*1.E4,FCCNR2(KR),KR=1,NKR)
	if (masterproc) PRINT 200, CONCCCNIN,CONTCCNIN
	DO KR=1,NKR
	   ROCCN(KR)=ROCCN03
	   X0=X0CCN*2.**(KR-1)
	   R0=(3.*X0/4./3.141593/ROCCN(KR))**DEG01
	   XCCN(KR)=X0
	   RCCN(KR)=R0
           RCCNKR_CM=R0
! CCN SPECTRUM 

! UNCOMMENT END IF  FOR CONTINENTAL/COMMENT FOR MARITIME
           IF(R0.LE.(0.4E-4)) THEN

           S_KR=A2/RCCNKR_CM**1.5
           FCCNR3(KR)=1.5*ACCN3*BCCN3*S_KR**BCCN3

! for maritime clouds in case of lack of small aerosol(ccn)
! particles we introduce maximum of supersaturation(SMAX) that
! if S > SMAX, % then FCCNR(KR)=0.

! UNCOMENT THESE TWO LINES FOR MARITIME/COMMENT FOR CONTINENTAL
!          SMAX=1.1
!          IF(S_KR.GE.SMAX) FCCNR(KR)=0.

	     CONTCCNIN=CONTCCNIN+COL*FCCNR3(KR)*R0*R0*R0
             CONCCCNIN=CONCCCNIN+COL*FCCNR3(KR)
          
! UNCOMMENT FOR CONTINENTAL/COMMENT FOR MARITIME
          ELSE
            FCCNR3(KR)=0.
          END IF
        END DO
        if (masterproc) PRINT*,    'RCCN(KR)*1.E4,FCCNR3(KR),KR=1,NKR'
        if (masterproc) PRINT 101, (RCCN(KR)*1.E4,FCCNR3(KR),KR=1,NKR)
	if (masterproc) PRINT 200, CONCCCNIN,CONTCCNIN

        fccnr_mp(:) = 0.
!        fcinr_mp(:) = 0.


! For ISDAC F31 Liu & Earle (used for ISDAC intercomparison)
!       DO KR=1,NKR
!        fccnr_mp(kr) = 207./(sqrt(2*3.1416)*log(1.5))    &
!     &                *exp(-(log(rccn(kr)*1.E4/0.1))**2/2.0/(log(1.5))**2)
!        fccnr_mp(kr) = fccnr_mp(kr)+ 8.5/(sqrt(2*3.1416)*log(2.45))      &
!     &                *exp(-(log(rccn(kr)*1.E4/0.35))**2/2.0/(log(2.45))**2)
!       enddo

!MO 4/19/16: For Case 4 (VOCALS) for the Int. Cloud Modeling Workshop 2016
!       DO KR=1,NKR
!        fccnr_mp(kr) = 300./(sqrt(2*3.1416)*log(1.4))    &
!    &                *exp(-(log(rccn(kr)*1.E4/0.04))**2/2.0/(log(1.4))**2)
!!        fccnr_mp(kr) = 50./(sqrt(2*3.1416)*log(1.5))    &
!!     &                *exp(-(log(rccn(kr)*1.E4/0.04))**2/2.0/(log(1.5))**2)
!       enddo


! For DYCOMS-II RF01 Low (Andrejczuk, personal communication and Andrejczuk et al., 2008) 
!       DO KR=1,NKR
!        fccnr_mp(kr) = 125./(sqrt(2*3.1416)*log(1.2))    &
!     &                *exp(-(log(rccn(kr)*1.E4/0.011))**2/2.0/(log(1.2))**2)
!        fccnr_mp(kr) = fccnr_mp(kr)+ 65.0/(sqrt(2*3.1416)*log(1.7))      &
!     &                *exp(-(log(rccn(kr)*1.E4/0.06))**2/2.0/(log(1.7))**2)
!       enddo

! For DYCOMS-II RF01 High (Andrejczuk, personal communication and Andrejczuk et al., 2008) 
!       DO KR=1,NKR
!        fccnr_mp(kr) = 125./(sqrt(2*3.1416)*log(1.2))    &
!     &                *exp(-(log(rccn(kr)*1.E4/0.011))**2/2.0/(log(1.2))**2)
!        fccnr_mp(kr) = fccnr_mp(kr)+ 1170.0/(sqrt(2*3.1416)*log(1.7))      &
!     &                *exp(-(log(rccn(kr)*1.E4/0.06))**2/2.0/(log(1.7))**2)
!       enddo

! For CHAPS
!        DO KR=1,NKR
!         fccnr_mp(kr) = 250./(sqrt(2*3.1416)*log(1.4))    &
!      &                *exp(-(log(rccn(kr)*2.E4/0.25))**2/2.0/(log(1.4))**2)
!         fccnr_mp(kr) =fccnr_mp(kr)+ 2400.0/(sqrt(2*3.1416)*log(1.7))      &
!      &                *exp(-(log(rccn(kr)*2.E4/0.075))**2/2.0/(log(1.7))**2)
! !        fcinr_mp(kr) =550.0/(sqrt(2*3.1416)*log(1.35))      &
! !     &                *exp(-(log(rccn(kr)*2.E4/0.023))**2/2.0/(log(1.35))**2)
!        enddo

! Adapted for LASSO-ENA test case runs
! For DYCOMS-II RF01 Low (Andrejczuk, personal communication and Andrejczuk et al., 2008) 
      DO KR=1,NKR
       fccnr_mp(kr) = 125./(sqrt(2*3.1416)*log(1.2))    &
    &                *exp(-(log(rccn(kr)*1.E4/0.011))**2/2.0/(log(1.2))**2)
       fccnr_mp(kr) = fccnr_mp(kr)+ 65.0/(sqrt(2*3.1416)*log(1.7))      &
    &                *exp(-(log(rccn(kr)*1.E4/0.06))**2/2.0/(log(1.7))**2)
      enddo

        if (masterproc) PRINT*,    'RCCN(KR)*1.E4,FCCNR_mp(KR),KR=1,NKR'
        if (masterproc) PRINT 101, (RCCN(KR)*1.E4,FCCNR_mp(KR),KR=1,NKR)

!        STOP
         CALL BREAKINIT
!        CALL TWOINITMXVAR

        if (masterproc) then
	PRINT *, '**** MIN CCN RADIUS,MASS & DENSITY ***'
	PRINT 200, R0CCN,X0CCN,ROCCN0
	PRINT *, '*********  CONT CCN CONCENTRATION & MASS *******'
	PRINT 200, CONCCCNIN,CONTCCNIN
	PRINT *, '*********  DROP RADII *******'
	PRINT 200, DROPRADII
	PRINT *, '*********  CCN RADII *******'
	PRINT 200, RCCN
	PRINT *, '********* CCN MASSES *******'
	PRINT 200, XCCN
	PRINT *, '********* INITIAL CCN DISTRIBUTION *******'
	PRINT 200, FCCNR0
	PRINT *, '********* INITIAL CCN2 DISTRIBUTION *******'
	PRINT 200, FCCNR2
	PRINT *, '********* INITIAL CCN3 DISTRIBUTION *******'
	PRINT 200, FCCNR3


  100	FORMAT(10I4)
  101   FORMAT(3X,F7.5,E13.5)
  102	FORMAT(4E12.4)
  105	FORMAT(A48)
  106	FORMAT(A80)
  123	FORMAT(3E12.4,3I4)
  200	FORMAT(6E13.5)
  201   FORMAT(6D13.5)
  300	FORMAT(8E14.6) 
  301   FORMAT(3X,F8.3,3X,E13.5)
  302   FORMAT(5E13.5)
        print*, 'dt = ',dt
!       if (IFREST)THEN
!       dtime=dt*0.5
!       else
!       END IF
        print*, 'dtime = ',dt/ncond
        end if ! masterproc

!check
!        call kernals(dt)
!  collision was called every 3*dt
        call kernals(dt*3.0)
! add

  if(nrestart.eq.0) then     ! initialize arrays for the initial run

     micro_field = 0.
     qc = 0.
     qr = 0.
     qi = 0.
     qg = 0.
     qs = 0.
     qh = 0.
     do k=1,nzm
       do j=1,ny
         do i=1,nx
           qt(i,j,k) = q0(k)             ! Initialize total water mixing ratio
           tabs(i,j,k) = t(i,j,k) - gamaz(k) + fac_cond * (qc(i,j,k)+ qr(i,j,k)) &
    &                  + fac_sub * (qi(i,j,k)+ qs(i,j,k)+qg(i,j,k)+qh(i,j,k))        
           if (i==20.and.j==1.and.k==1) print*,'ptemp,temp', t(i,j,k), tabs(i,j,k)
         end do    
       end do     
     end do       

     t_old = tabs
     qv_old = qv

       fluxbmk = 0.
       fluxtmk = 0.
                                                                                                             
        if(docloud) then
         call micro_diagnose()
        endif
  else                          !  nrestart /= 0  then this is a continuation run
       do k=1,nzm
        do j=1,ny
         do i=1,nx
          DO  KR=1,NKR
          IF (KR.LT.KRDROP)THEN
          QC(I,j,k)=QC(I,j,k) &
     &      +COL*ffcd(I,j,k,KR)
        ELSE
          QR(I,j,k)=QR(I,j,k) &
     &      +COL*ffcd(I,j,k,KR)
         END IF
         END DO
         IF (ICEPROCS.EQ.1)THEN
         DO  KR=1,33
         QIC(I,j,k)=QIC(I,j,k) &
     &   +COL*ffic(I,j,k,KR)
         END DO
         DO  KR=1,33
         QIP(I,j,k)=QIP(I,j,k) &
     &   +COL*ffip(I,j,k,KR)
         END DO
         DO  KR=1,33
         QID(I,j,k)=QID(I,j,k) &
     &   +COL*ffid(I,j,k,KR)
         END DO
         QI(I,j,k) = QID(I,j,k)+QIP(I,j,k)+QIC(I,j,k)
                                                                                                         
         DO  KR=1,33
         QS(I,j,k)=QS(I,j,k) &
     &   +COL*ffsn(I,j,k,KR)
         END DO
         DO  KR=1,33
         QG(I,j,k)=QG(I,j,k) &
     &   +COL*ffgl(I,j,k,KR)
         END DO
         DO  KR=1,33
         QH(I,j,k)=QH(I,j,k) &
     &   +COL*ffhl(I,j,k,KR)
         END DO
         END IF

         end do
        end do
       end do                  

     do k=1,nzm
       do j=1,ny
         do i=1,nx
           tabs(i,j,k) = t(i,j,k) - gamaz(k) + fac_cond * (qc(i,j,k)+ qr(i,j,k)) &
    &                  + fac_sub * (qi(i,j,k)+ qs(i,j,k)+qg(i,j,k)+qh(i,j,k))
         end do
       end do
     end do


     t_old = tabs
     qv_old= qv
  end if


        return
      END SUBROUTINE micro_init

!--------------------------------------------------------------------------------------
! adapted end


! Below, the required subroutines and functions are given:
!!! fill-in surface and top boundary fluxes:
!
! Obviously, for liquid/ice water variables those fluxes are zero. They are not zero
! only for water vapor variable and, possibly, for CCN and IN if you have those.

subroutine micro_flux()

  use vars, only: fluxbq, fluxtq

  fluxbmk(:,:,index_water_vapor) = fluxbq(:,:)
  fluxtmk(:,:,index_water_vapor) = fluxtq(:,:)

end subroutine micro_flux

subroutine micro_diagnose()
                                                                                                             
   use vars
   integer i,j,k
                                                                                                             
   do k=1,nzm
    do j=1,ny
     do i=1,nx
       qv(i,j,k) = qt(i,j,k) - (qc(i,j,k)+qr(i,j,k)+qi(i,j,k)+qs(i,j,k)+    &
                   qg(i,j,k)+qh(i,j,k))
!       qcl(i,j,k) = qc(i,j,k)
!       qci(i,j,k) = qi(i,j,k)
!       qpl(i,j,k) = qr(i,j,k)
!       qpi(i,j,k) = qs(i,j,k) + qg(i,j,k) + qh(i,j,k)

!mo For DYCOMS and ISDAC get all condenced liquid water into the "cloud water" 
!mo   for radiation and other purposes
!
       qcl(i,j,k) = qc(i,j,k) + qr(i,j,k)
       qci(i,j,k) = qi(i,j,k)
       qpl(i,j,k) = 0
       qpi(i,j,k) = qs(i,j,k) + qg(i,j,k) + qh(i,j,k)
                                                                                               
     end do
    end do
   end do                                                                                                             
  return
end subroutine micro_diagnose

!!! functions to compute terminal velocity for precipitating variables:
!
! you need supply functions to compute terminal velocity for all of your 
! precipitating prognostic variables. Note that all functions should
! compute vertical velocity given two microphysics parameters var1, var2, 
! and temperature, and air density (single values, not arrays). Var1 and var2 
! are some microphysics variables like water content and concentration.
! Don't change the number of arguments or their meaning!

!real function term_vel_ql(qr,nr,tabs,rho)
!
!  real, intent(in) ::  qr, nr, tabs, rho
!  term_vel_ql = 0.
!end function term_vel_ql

real function term_vel_ffcd(i,j,k,nbin)  ! check1:pres !change for version 6.7.5
                                                                                             
  real, intent(in) ::  i,j,k
  integer nbin
!  nbin = int(bin)                         ! bin index
!  if (qbin > ql_prec) then
                                                                                             
    term_vel_ffcd = VR1(nbin)*0.01
!   else
!    term_vel_ffcd = 0.
!  end if
end function term_vel_ffcd
                                                                                             
real function term_vel_ffic(i,j,k, nbin)
                                                                                             
  real, intent(in) ::  i,j,k
  integer nbin
!  nbin = int(bin)                         ! bin index
!  if (qbin > ql_prec) then
                                                                                             
!    term_vel_ffic = VR2(nbin,1)*0.01*SQRT(1.E5/(rho*287.1*tabs))
    term_vel_ffic = VR2(nbin,1)*0.01  ! the vertical correction has been done in precip_fall
!   else
!    term_vel_ffic = 0.
!  end if
end function term_vel_ffic
                                                                                             
real function term_vel_ffip(i,j,k,nbin)
                                                                                             
  real, intent(in) ::  i,j,k
  integer nbin
!  nbin = int(bin)                         ! bin index
!  if (qbin > ql_prec) then
    term_vel_ffip = VR2(nbin,2)*0.01
!   else
!    term_vel_ffip = 0.
!  end if
end function term_vel_ffip
                                                                                             
real function term_vel_ffid(i,j,k,nbin)
                                                                                             
  real, intent(in) ::  i,j,k
  integer nbin
!  nbin = int(bin)                         ! bin index
!  if (qbin > ql_prec) then
                                                                                             
    term_vel_ffid = VR2(nbin,3)*0.01
!   else
!    term_vel_ffid = 0.
!  end if
end function term_vel_ffid
                                                                                             
real function term_vel_ffsn(i,j,k,nbin)
                                                                                             
  real, intent(in) ::  i,j,k
  integer nbin
!  nbin = int(bin)                         ! bin index
!  if (qbin > ql_prec) then
                                                                                             
    term_vel_ffsn = VR3(nbin)*0.01
!   else
!    term_vel_ffsn = 0.
!  end if
end function term_vel_ffsn
                                                                                             
real function term_vel_ffgl(i,j,k,nbin)
                                                                                             
  real, intent(in) ::  i,j,k
  integer nbin
!  nbin = int(bin)                         ! bin index
!  if (qbin > ql_prec) then
    term_vel_ffgl = VR4(nbin)*0.01
!   else
!    term_vel_ffgl = 0.
!  end if
end function term_vel_ffgl
                                                                                             
real function term_vel_ffhl(i,j,k,nbin)
                                                                                             
  real, intent(in) ::  i,j,k
  integer nbin
!  nbin = int(bin)                         ! bin index
!  if (qbin > ql_prec) then
                                                                                             
    term_vel_ffhl = VR5(nbin)*0.01
!   else
!    term_vel_ffhl = 0.
!  end if
end function term_vel_ffhl
                                                                                             
!real function term_vel_Nr(qr,nr,tabs,rho)
!  real, intent(in) ::  qr, dummy, tabs, rho
!end function term_vel_Nr
                                                                                             
!real function term_vel_qs(qs,ns,tabs,rho)
!  real, intent(in) ::  qr, dummy, tabs, rho
!end function term_vel_qs
                                                                                             
! etc.
                                                                                             
                                                                                             
subroutine micro_precip_fall()
                                                                                             
use vars
! use vars, only : s_ar
                                                                                             
! before calling precip_fall() for each of falling prognostic variables,
! you need to set hydro_type and omega(:,:,:) variables.
! hydro_type can have four values:
! 0 - variable is liquid water mixing ratio
! 1 - hydrometeor is ice mixing ratio ! 2 - hydrometeor is mixture-of-liquid-and-ice mixing ratio. (As in original SAM microphysics).
! 3 - variable is not mixing ratio, but, for example, rain drop concentration
! OMEGA(:,:,:) is used only for hydro_type=2, and is the fraction of liquid phase (0-1).
! for our hypothetical case, there is no mixed hydrometeor, so omega is not actually used.
                                                                                             
  integer hydro_type
!  real omega(nx,ny,nzm)
  real omega(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
  real dummy(1)
  integer i,j,k,m,kr
  integer, parameter :: its=1, ite=nx, jts=1,     &
                         jte=ny, kts=1, kte=nzm    

                                                                                             
! Initialize arrays that accumulate surface precipitation flux
                                                                                             
 if(mod(nstep-1,nstatis).eq.0.and.icycle.eq.1) then
   do j=1,ny
    do i=1,nx
     precsfc(i,j)=0.
    end do
   end do
   do k=1,nzm
    precflux(k) = 0.
   end do
 end if
                                                                                             
 do k = 1,nzm ! Initialize arrays which hold precipitation fluxes for stats.
    qpfall(k)=0.
    tlat(k) = 0.
 end do

!======================================================
!=== 26oct2009 (MO)

      DO k = kts,kte
      DO j = jts,jte
      DO i = its,ite

        QC(I,j,k)=0
        QR(I,j,k)=0
        QI(I,j,k)=0
        QIC(I,j,k)=0
        QIP(I,j,k)=0
        QID(I,j,k)=0
        QS(I,j,k)=0
        QG(I,j,k)=0
        QH(I,j,k)=0

        DO  KR=1,NKR
          IF (KR.LT.KRDROP)THEN
            QC(I,j,k)=QC(I,j,k) &
     &        +ffcd(I,j,k,KR)
          ELSE
            QR(I,j,k)=QR(I,j,k) &
     &        +ffcd(I,j,k,KR)
          END IF
        END DO
        IF (QC(I,j,k).LT.1.E-8)QC(I,J,k)=0.0
        IF (QR(I,j,k).LT.1.E-8)QR(I,j,k)=0.0
        QC(I,j,k)= QC(I,j,k)* col
        QR(I,j,k)= qr(i,j,k)* col
   
        IF (ICEPROCS.EQ.1)THEN
         DO  KR=1,33
          QIC(I,j,k)=QIC(I,j,k) &
     &      +ffic(I,j,k,KR)
          QIP(I,j,k)=QIP(I,j,k) &
     &      +ffip(I,j,k,KR)
          QID(I,j,k)=QID(I,j,k) &
     &      +ffid(I,j,k,KR)
          QS(I,j,k)=QS(I,j,k) &
     &      +ffsn(I,j,k,KR)
          QG(I,j,k)=QG(I,j,k) &
     &      +ffgl(I,j,k,KR)
          QH(I,j,k)=QH(I,j,k) &
     &      +ffhl(I,j,k,KR)
         END DO
         QI(I,j,k) = QID(I,j,k)+QIP(I,j,k)+QIC(I,j,k)

         QI(I,j,k)= QI(I,j,k) * col 
         QIC(I,j,k)= QIC(I,j,k) * col
         QIP(I,j,k)= QIP(I,j,k) * col
         QID(I,j,k)= QID(I,j,k) * col
         QS(I,j,k)= QS(I,j,k) * col
         QG(I,j,k)= QG(I,j,k) * col
         QH(I,j,k)= QH(I,j,k) * col
        ENDIF

        qv(i,j,k) = qt(i,j,k) - (qc(i,j,k)+qr(i,j,k)+qi(i,j,k)+qs(i,j,k)+    &
                   qg(i,j,k)+qh(i,j,k))
!       if (qv(i,j,k) .le. 1.e-15) qv(i,j,k)=1.e-15

      ENDDO
      ENDDO
      ENDDO
!======================================================
                                                                                             
! Compute sedimentation of falling variables:
!
! COL factor is added so the sum of the fluxes computed in PRECIP_FALL 
! corresponds to kg/kg for a proper accounting of teh latent heat
! (26 oct 2009 [MO])
      sedl_tend(:,:,:) = 0.0
      sedi_tend(:,:,:) = 0.0
      DO k = kts,kte
      DO j = jts,jte
      DO i = its,ite
      qlbf_sed(i,j,k)=0.
      qibf_sed(i,j,k)=0.
       DO  KR=1,NKR
        qlbf_sed(i,j,k)=qlbf_sed(i,j,k) &
        &      +ffcd(i,j,k,KR)*COL
        qibf_sed(i,j,k)=qibf_sed(i,j,k) &
        &      +COL*ffic(i,j,k,KR)+     &
        &      +COL*ffip(i,j,k,KR)+     &
        &      +COL*ffid(i,j,k,KR)+     &
        &      +COL*ffsn(i,j,k,KR)+     &
        &      +COL*ffgl(i,j,k,KR)      
       ENDDO     
      ENDDO
      ENDDO
      ENDDO
                                                                                         
 hydro_type=0      ! 0 - all liquid, 1 - all ice, 2 - mixed
 do m=1,ncd
   omega(:,:,:) = ffcd(:,:,:,m) * col
   dummy(1) = real(m)
!   call precip_fall(omega, dummy, term_vel_ffcd, hydro_type, omega, m)
   call precip_fall(omega, term_vel_ffcd, hydro_type, omega, m)  ! for version 6.7.5
   ffcd(:,:,:,m) = omega(:,:,:) / col
 end do
                                                                                             
if (ICEPROCS.EQ.1) then
hydro_type=1
 do m=1,ncd
   omega(:,:,:) = ffic(:,:,:,m) * col
   dummy(1) = real(m)
   call precip_fall(omega, term_vel_ffic, hydro_type, omega,m)
   ffic(:,:,:,m) = omega(:,:,:) / col
 end do
 do m=1,ncd
   omega(:,:,:) = ffip(:,:,:,m) * col
   dummy(1) = real(m)
   call precip_fall(omega, term_vel_ffip, hydro_type, omega,m)
   ffip(:,:,:,m) = omega(:,:,:) / col
 end do
 do m=1,ncd
   omega(:,:,:) = ffid(:,:,:,m) * col
   dummy(1) = real(m)
   call precip_fall(omega, term_vel_ffid, hydro_type, omega,m)
   ffid(:,:,:,m) = omega(:,:,:) / col
 end do
 do m=1,ncd
   omega(:,:,:) = ffsn(:,:,:,m) * col
   dummy(1) = real(m)
   call precip_fall(omega, term_vel_ffsn, hydro_type, omega,m)
   ffsn(:,:,:,m) = omega(:,:,:) / col
 end do
 do m=1,ncd
   omega(:,:,:) = ffgl(:,:,:,m) * col
   dummy(1) = real(m)
   call precip_fall(omega, term_vel_ffgl, hydro_type, omega,m)
   ffgl(:,:,:,m) = omega(:,:,:) / col
 end do
 do m=1,ncd
   omega(:,:,:) = ffhl(:,:,:,m) * col
   dummy(1) = real(m)
   call precip_fall(omega, term_vel_ffhl, hydro_type, omega,m)
   ffhl(:,:,:,m) = omega(:,:,:) / col
 end do
                                                                                             
endif

      DO k = kts,kte
      DO j = jts,jte
      DO i = its,ite
      qlaf_sed(i,j,k)=0.
      qiaf_sed(i,j,k)=0.
       DO  KR=1,NKR
        qlaf_sed(i,j,k)=qlaf_sed(i,j,k) &
        &      +ffcd(i,j,k,KR)*COL
        qiaf_sed(i,j,k)=qiaf_sed(i,j,k) &
        &      +COL*ffic(i,j,k,KR)+     &
        &      +COL*ffip(i,j,k,KR)+     &
        &      +COL*ffid(i,j,k,KR)+     &
        &      +COL*ffsn(i,j,k,KR)+     &
        &      +COL*ffgl(i,j,k,KR)
       ENDDO
      ENDDO
      ENDDO
      ENDDO
   sedl_tend(:,:,:) = (qlaf_sed(:,:,:)-qlbf_sed(:,:,:))/dt
   sedi_tend(:,:,:) = (qiaf_sed(:,:,:)-qibf_sed(:,:,:))/dt

! print*, 'sed tend', maxval(sedl_tend(:,:,:)), maxval(sedi_tend(:,:,:))

! call precip_fall(qr, dummy, term_vel_qr, hydro_type, omega)
! hydro_type=3
! call precip_fall(Nr, dummy, term_vel_Nr, hydro_type, omega)
! hydro_type=1
! call precip_fall(qs, dummy, term_vel_qs, hydro_type, omega)
! hydro_type=3
! call precip_fall(Ns, dummy, term_vel_Ns, hydro_type, omega)
! hydro_type=1
! call precip_fall(qg, dummy, term_vel_qg, hydro_type, omega)
! hydro_type=3
! call precip_fall(Ng, dummy, term_vel_Ng, hydro_type, omega)


!==Update Qt ==================================================
!=== 26oct2009 (MO)

      DO k = kts,kte
      DO j = jts,jte
      DO i = its,ite

        QC(I,j,k)=0
        QR(I,j,k)=0
        QI(I,j,k)=0
        QIC(I,j,k)=0
        QIP(I,j,k)=0
        QID(I,j,k)=0
        QS(I,j,k)=0
        QG(I,j,k)=0
        QH(I,j,k)=0

        DO  KR=1,NKR
          IF (KR.LT.KRDROP)THEN
            QC(I,j,k)=QC(I,j,k) &
     &        +ffcd(I,j,k,KR)
          ELSE
            QR(I,j,k)=QR(I,j,k) &
     &        +ffcd(I,j,k,KR)
          END IF
        END DO
        IF (QC(I,j,k).LT.1.E-8)QC(I,J,k)=0.0
        IF (QR(I,j,k).LT.1.E-8)QR(I,j,k)=0.0
        QC(I,j,k)= QC(I,j,k)* col
        QR(I,j,k)= qr(i,j,k)* col
   
        IF (ICEPROCS.EQ.1)THEN
         DO  KR=1,33
          QIC(I,j,k)=QIC(I,j,k) &
     &      +ffic(I,j,k,KR)
          QIP(I,j,k)=QIP(I,j,k) &
     &      +ffip(I,j,k,KR)
          QID(I,j,k)=QID(I,j,k) &
     &      +ffid(I,j,k,KR)
          QS(I,j,k)=QS(I,j,k) &
     &      +ffsn(I,j,k,KR)
          QG(I,j,k)=QG(I,j,k) &
     &      +ffgl(I,j,k,KR)
          QH(I,j,k)=QH(I,j,k) &
     &      +ffhl(I,j,k,KR)
         END DO
         QI(I,j,k) = QID(I,j,k)+QIP(I,j,k)+QIC(I,j,k)

         QI(I,j,k)= QI(I,j,k) * col 
         QIC(I,j,k)= QIC(I,j,k) * col
         QIP(I,j,k)= QIP(I,j,k) * col
         QID(I,j,k)= QID(I,j,k) * col
         QS(I,j,k)= QS(I,j,k) * col
         QG(I,j,k)= QG(I,j,k) * col
         QH(I,j,k)= QH(I,j,k) * col

        ENDIF

        qt(i,j,k) = qv(i,j,k) + (qc(i,j,k)+qr(i,j,k)+qi(i,j,k)+qs(i,j,k)+    &
                   qg(i,j,k)+qh(i,j,k))
!       if (qv(i,j,k) .le. 1.e-15) qv(i,j,k)=1.e-15

      ENDDO
      ENDDO
      ENDDO
!======================================================

! compute surface precipitation area fraction statistics
 do j=1,ny
   do i=1,nx
!     if(ql(i,j,1).gt.1.e-6) s_ar=s_ar+dtfactor
     if((qr(i,j,1)+qs(i,j,1)+qg(i,j,1)+qh(i,j,1)).gt.1.e-6) s_ar=s_ar+dtfactor
   end do
 end do

end subroutine micro_precip_fall


!----------------------------------------------------------------------
!!! Save fields at the begining of the dynamical time step for 
!!! microphysics calculations. These are needed to compute supersaturation.

!subroutine micro_fields_save()
!use vars, only:  t

!qt_sv = qt
!t_sv  = t

!return
!end subroutine micro_fields_save


!----------------------------------------------------------------------
!!! Initialize the list of microphysics statistics that will be outputted
!!  to *.stat statistics file

subroutine micro_hbuf_init(namelist,deflist,unitlist,status,average_type,count,trcount)


   character(*) namelist(*), deflist(*), unitlist(*)
   integer status(*),average_type(*),count,trcount
   integer ntr

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QTFLUX'
   deflist(count) = 'Nonprecipitating water flux (Total)'
   unitlist(count) = 'W/m2'
   status(count) = 1
   average_type(count) = 0
                                                                                                                
   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QTFLUXS'
   deflist(count) = 'Nonprecipitating-water flux (SGS)'
   unitlist(count) = 'W/m2'
   status(count) = 1
   average_type(count) = 0
                                                                                                                
   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QPFLUX'
   deflist(count) = 'Precipitating-water turbulent flux (Total)'
   unitlist(count) = 'W/m2'
   status(count) = 1
   average_type(count) = 0
                                                                                                                
   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QPFLUXS'
   deflist(count) = 'Precipitating-water turbulent flux (SGS)'
   unitlist(count) = 'W/m2'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'MQC'
   deflist(count) = 'Cloud water (microphysics)'
   unitlist(count) = 'g/m^3'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'MQR'
   deflist(count) = 'Rain water (microphysics)'
   unitlist(count) = 'g/m^3'
   status(count) = 1
   average_type(count) = 0
                                                                                                                
   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'MQI'
   deflist(count) = 'Ice crystal (microphysics)'
   unitlist(count) = 'g/m^3'
   status(count) = 1
   average_type(count) = 0
                                                                                                                
   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'MQS'
   deflist(count) = 'Snow water (microphysics)'
   unitlist(count) = 'g/m^3'
   status(count) = 1
   average_type(count) = 0
                                                                                                                
   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'MQG'
   deflist(count) = 'Graupel water (microphysics)'
   unitlist(count) = 'g/m^3'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'MQH'
   deflist(count) = 'Hail water (microphysics)'
   unitlist(count) = 'g/m^3'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'SSWMICRO'
   deflist(count) = 'Supersaturation with respect to water'
   unitlist(count) = '%'
   status(count) = 1
   average_type(count) = 0
                                                                                                       
   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'SSIMICRO'
   deflist(count) = 'Supersaturation with respect to ice'
   unitlist(count) = '%'
   status(count) = 1
   average_type(count) = 0

!... etc.


   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'NCN'
   deflist(count) = 'CCN concentration'
   unitlist(count) = 'cm-3'
   status(count) = 1    
   average_type(count) = 0

!Add for IN
   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'NIN'
   deflist(count) = 'IN concentration'
   unitlist(count) = 'L-3'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'NCD'
   deflist(count) = 'CD concentration'
   unitlist(count) = 'cm-3'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'NR'
   deflist(count) = 'Rain drop concentration'
   unitlist(count) = 'cm-3'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'NI'
   deflist(count) = 'ice crystal number concentration'
   unitlist(count) = 'cm-3'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'NS'
   deflist(count) = 'Snow number concentration'
   unitlist(count) = 'l-3'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'NG'
   deflist(count) = 'graupel concentration'
   unitlist(count) = 'l-3'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'NH'
   deflist(count) = 'Hail number concentration'
   unitlist(count) = 'l-3'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'VFI_MW'
   deflist(count) = 'mass-weighted ice fall velocity'
   unitlist(count) = 'm/s'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'DITEND'
   deflist(count) = 'Ice mixing ratio tendency due to dep/sub'
   unitlist(count) = 'g/kg/s'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'DLTEND'
   deflist(count) = 'Liquid mixing ratio tendency due to cond/evap'
   unitlist(count) = 'g/kg/s'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'SEDLTEND'
   deflist(count) = 'Liquid mixing ratio tendency due to sed'
   unitlist(count) = 'g/kg/s'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'SEDITEND'
   deflist(count) = 'Ice mixing ratio tendency due to sed'
   unitlist(count) = 'g/kg/s'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'FRZTEND'
   deflist(count) = 'Liquid mixing ratio tendency due to freezing'
   unitlist(count) = 'g/kg/s'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QVADVT'
   deflist(count) = 'Water vapor tendency due to advection'
   unitlist(count) = 'g/kg/s'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QLADVT'
   deflist(count) = 'Liquid water tendency due to advection'
   unitlist(count) = 'g/kg/s'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QIADVT'
   deflist(count) = 'Ice water tendency due to advection'
   unitlist(count) = 'g/kg/s'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QNI100'
   deflist(count) = 'Ice particle conc > 100 um'
   unitlist(count) = 'L-1'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'Ze'
   deflist(count) = 'Radar reflect. (include attenu.) for Qi>0.1 kg/kg'
   unitlist(count) = 'dBZ'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'VI_Ze'
   deflist(count) = 'Ze-weighted ice fall velocity for Qi>0.1 kg/kg'
   unitlist(count) = 'm/s'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QVTENDT'
   deflist(count) = 'Total water vapor tendency'
   unitlist(count) = 'g/kg/s'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QLTENDT'
   deflist(count) = 'Total liquid water tendency'
   unitlist(count) = 'g/kg/s'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QITENDT'
   deflist(count) = 'Total ice water tendency'
   unitlist(count) = 'g/kg/s'
   status(count) = 1
   average_type(count) = 0
!... etc.

end subroutine micro_hbuf_init

!----------------------------------------------------------------------
!!!! Collect microphysics history statistics (vertical profiles)
!! Note that only the fields declared in micro_hbuf_init() are allowed to
! be collected

subroutine micro_statistics()
  
  use vars
  use hbuffer, only: hbuf_put,hbuf_avg_put
  use params

  real tmp(2), factor_xy 
  real qcz(nzm), qiz(nzm), qrz(nzm), qsz(nzm), qgz(nzm), qhz(nzm),omg, sswz(nzm),ssiz(nzm)
  real qcdz(nzm), qcnz(nzm), qnrz(nzm), qniz(nzm), qnsz(nzm), qngz(nzm), qnhz(nzm)
! Add for IN
  real qinz(nzm)
! added for sheba output
  real vfi_mwz(nzm), diffuiratez(nzm), difful_tendz(nzm), sedl_tendz(nzm),sedi_tendz(nzm)
  real frzl_tendz(nzm), qvtend_advz(nzm), qltend_advz(nzm), qitend_advz(nzm)
  real qni100z(nzm), qvtendz(nzm), qltendz(nzm),qitendz(nzm)

  real factor_vol(nzm), ct1
  integer i,j,k,m

  factor_xy = 1./float(nx*ny)

  do k=1,nzm
      tmp(1) = dz/rhow(k)
      tmp(2) = tmp(1) / dtn
      mkwsb(k,1) = mkwsb(k,1) * tmp(1) * rhow(k) * lcond
      mkwle(k,1) = mkwle(k,1)*tmp(2)*rhow(k)*lcond + mkwsb(k,1)
      if(docloud.and.doprecip) then
        mkwsb(k,2) = mkwsb(k,2) * tmp(1) * rhow(k) * lcond
        mkwle(k,2) = mkwle(k,2)*tmp(2)*rhow(k)*lcond + mkwsb(k,2)
      endif
  end do

  do k=1,nzm
      tmp(1) = dz/rhow(k)
      tmp(2) = tmp(1) / dtn
      mkwsb(k,:) = mkwsb(k,:) * tmp(1) * rhow(k) * lcond
      mkwle(k,:) = mkwle(k,:)*tmp(2)*rhow(k)*lcond + mkwsb(k,:)
  end do
!  call hbuf_put('QVFLUX',mkwle(:,3),factor_xy)
!  call hbuf_put('QVFLUXS',mkwsb(:,3),factor_xy)
!  call hbuf_put('QCFLUX',mkwle(:,4),factor_xy)
!  call hbuf_put('QCFLUXS',mkwsb(:,4),factor_xy)
  call hbuf_put('QTFLUX',mkwle(:,1),factor_xy)
  call hbuf_put('QTFLUXS',mkwsb(:,1),factor_xy)
  call hbuf_put('QPFLUX',mkwle(:,2),factor_xy)
  call hbuf_put('QPFLUXS',mkwsb(:,2),factor_xy)
! for the statistics from radar simulator
!  if (masterproc) print*, 'Ze', maxval(Zez), maxval(V_doppiz)
!    call hbuf_put('Ze',Zez,1.)
!    call hbuf_put('VI_Ze',V_doppiz,1.)
! ... etc
do k=1,nzm
  factor_vol(k) = 1.e3*rho(k)  ! from g/g to g/m3
enddo
!-----------------------------------------------------------------------
! Average profiles

  do k=1,nzm
    qcz(k) = 0.
    qiz(k) = 0.
    qrz(k) = 0.
    qsz(k) = 0.
    qgz(k) = 0.
    qhz(k) = 0.
    sswz(k) = 0.
    ssiz(k) = 0.

    qcdz(k) = 0. 
    qcnz(k) = 0.
    qnrz(k) = 0.
    qniz(k) = 0.
    qnsz(k) = 0.
    qngz(k) = 0.
    qnhz(k) = 0.
!Add for IN
    qinz(k) = 0.
! added for sheba output
    vfi_mwz(k) = 0.
    diffuiratez(k) = 0.
    difful_tendz(k) = 0.
    sedl_tendz(k) = 0.
    sedi_tendz(k) = 0.
    frzl_tendz(k) = 0.
    qvtend_advz(k) = 0.
    qltend_advz(k) = 0.
    qitend_advz(k) = 0.
    qvtendz(k)  = 0.
    qltendz(k)  = 0.
    qitendz(k)  = 0.

    qni100z(k) = 0.

    ct1=0.0

    do j=1,ny
    do i=1,nx
      qcz(k)=qcz(k)+qc(i,j,k)*factor_vol(k)
      qrz(k)=qrz(k)+qr(i,j,k)*factor_vol(k)
      qiz(k)=qiz(k)+qi(i,j,k)*factor_vol(k)
      qsz(k)=qsz(k)+qs(i,j,k)*factor_vol(k)
      qgz(k)=qgz(k)+qg(i,j,k)*factor_vol(k)
      qhz(k)=qhz(k)+qh(i,j,k)*factor_vol(k)
      sswz(k)=sswz(k)+ssatw(i,j,k)
      ssiz(k)=ssiz(k)+ssati(i,j,k)
 
      qcnz(k)=qcnz(k)+qna(i,j,k)
      qcdz(k)=qcdz(k)+qnc(i,j,k)
      qnrz(k)=qnrz(k)+qnr(i,j,k)
      qniz(k)=qniz(k)+qni(i,j,k)
      qnsz(k)=qnsz(k)+qns(i,j,k)
      qngz(k)=qngz(k)+qng(i,j,k)
      qnhz(k)=qnhz(k)+qnh(i,j,k)

      qinz(k)=qinz(k)+qnin(i,j,k)
! added for sheba output  
      diffuiratez(k)=diffuiratez(k)+diffui_tend(i,j,k)
      difful_tendz(k)=difful_tendz(k)+difful_tend(i,j,k)
      sedl_tendz(k)=sedl_tendz(k)+sedl_tend(i,j,k)
      sedi_tendz(k)=sedi_tendz(k)+sedi_tend(i,j,k)
      frzl_tendz(k)=frzl_tendz(k)+frzl_tend(i,j,k)
      qvtend_advz(k)=qvtend_advz(k)+qvtend_adv(i,j,k)
      qltend_advz(k)=qltend_advz(k)+qltend_adv(i,j,k)
      qitend_advz(k)=qitend_advz(k)+qitend_adv(i,j,k)
      qvtendz(k)=qvtendz(k)+qvtend(i,j,k)
      qltendz(k)=qltendz(k)+qltend(i,j,k)
      qitendz(k)=qitendz(k)+qitend(i,j,k)

      if ((qi(i,j,k)+qs(i,j,k)+qg(i,j,k)+qh(i,j,k)) >1.e-7) then
      vfi_mwz(k)=vfi_mwz(k)+vfice_mw(i,j,k)
      qni100z(k)=qni100z(k)+qni100(i,j,k)
      ct1=ct1+1.0
      endif

    end do
    end do
     if (ct1 > 0.) then
      vfi_mwz(k)=vfi_mwz(k)/ct1
      qni100z(k)=qni100z(k)/ct1
     endif
  end do

!print*, 'VF', maxval(vfice_mw),maxval(vfi_mwz) 
!print*, 'NI100', maxval(qni100),maxval(qni100z)
!print*, 'DIFFI', maxval(diffui_tend),maxval(diffuiratez)
!print*, 'DIFFL', maxval(difful_tend),maxval(difful_tendz)
!print*, 'SEDL', maxval(sedl_tend),maxval(sedl_tendz)
!print*, 'SEDI', maxval(sedi_tend),maxval(sedi_tendz)
!print*, 'FRZL', maxval(frzl_tend),maxval(frzl_tendz)
!print*, 'QVADV', maxval(qvtend_adv),maxval(qvtend_advz)

  call hbuf_put('MQC',qcz,factor_xy)       ! g/m3
  call hbuf_put('MQR',qrz,factor_xy)
  call hbuf_put('MQI',qiz,factor_xy)
  call hbuf_put('MQS',qsz,factor_xy)
  call hbuf_put('MQG',qgz,factor_xy)
  call hbuf_put('MQH',qhz,factor_xy)
  call hbuf_put('SSWMICRO',sswz,factor_xy)
  call hbuf_put('SSIMICRO',ssiz,factor_xy)

  call hbuf_put('NCN',qcnz,1.0*factor_xy)

  call hbuf_put('NIN',qinz,1.e3*factor_xy)

  call hbuf_put('NCD',qcdz,1.0*factor_xy)
  call hbuf_put('NR',qnrz,1.0*factor_xy)
  call hbuf_put('NI',qniz,1.0*factor_xy)
  call hbuf_put('NS',qnsz,1.e3*factor_xy)
  call hbuf_put('NG',qngz,1.e3*factor_xy)
  call hbuf_put('NH',qnhz,1.e3*factor_xy)
! added for sheba output
  call hbuf_put('VFI_MW',vfi_mwz,1.)
  call hbuf_put('QNI100',qni100z,1.e3)
  call hbuf_put('DITEND',diffuiratez,1.e3*factor_xy)  ! g/kg/s
  call hbuf_put('DLTEND',difful_tendz,1.e3*factor_xy)  ! g/kg/s
  call hbuf_put('SEDLTEND',sedl_tendz,1.e3*factor_xy)  ! g/kg/s
  call hbuf_put('SEDITEND',sedi_tendz,1.e3*factor_xy)  ! g/kg/s
  call hbuf_put('FRZTEND',frzl_tendz,1.e3*factor_xy)  ! g/kg/s
  call hbuf_put('QVADVT',qvtend_advz,1.e3*factor_xy)  ! g/kg/s
  call hbuf_put('QLADVT',qltend_advz,1.e3*factor_xy)  ! g/kg/s
  call hbuf_put('QIADVT',qitend_advz,1.e3*factor_xy)  ! g/kg/s
  call hbuf_put('QVTENDT',qvtendz,1.e3*factor_xy)  ! g/kg/s
  call hbuf_put('QLTENDT',qltendz,1.e3*factor_xy)  ! g/kg/s
  call hbuf_put('QITENDT',qitendz,1.e3*factor_xy)  ! g/kg/s
!  call hbuf_avg_put('QC',qc,1,nx,1,ny,nzm,1.e+3)     ! [g/kg]
!  call hbuf_avg_put('QR',qr,1,nx,1,ny,nzm,1.e+3)     ! [g/kg]
!  call hbuf_avg_put('QI',qi,1,nx,1,ny,nzm,1.e+3)     ! [g/kg]
!  call hbuf_avg_put('QS',qs,1,nx,1,ny,nzm,1.e+3)     ! [g/kg]
!  call hbuf_avg_put('QG',qg,1,nx,1,ny,nzm,1.e+3)     ! [g/kg]
!  call hbuf_avg_put('QH',qh,1,nx,1,ny,nzm,1.e+3)     ! [g/kg]

!  call hbuf_avg_put('NCN',qna,1,nx,1,ny,nzm,1.0)  ! [#/cm^3]
!  call hbuf_avg_put('NCD',qnc,1,nx,1,ny,nzm, 1.0) ! [#/cm^3]
!  call hbuf_avg_put('NR',qnr,1,nx,1,ny,nzm, 1.0) ! [#/cm^3]
!  call hbuf_avg_put('NI',qni,1,nx,1,ny,nzm, 1.0) ! [#/l^3]
!  call hbuf_avg_put('NS',qns,1,nx,1,ny,nzm, 1.0) ! [#/l^3]
!  call hbuf_avg_put('NG',qng,1,nx,1,ny,nzm, 1.0) ! [#/l^3]
!  call hbuf_avg_put('NH',qnh,1,nx,1,ny,nzm, 1.0) ! [#/l^3]

end subroutine micro_statistics

!-----------------------------------------------------------------------
! Function that computes total water in a domain:
! Don't change this one.

real function total_water() 

  use vars, only : nstep,nprint,adz,dz,rho

  integer k,m

  total_water = 0.
  if(mod(nstep,nprint).ne.0) return

  do m=1,nmicro_fields

   if(flag_wmass(m).eq.1) then

    do k=1,nzm
      total_water = total_water + &
       sum(micro_field(1:nx,1:ny,k,m))*adz(k)*dz*rho(k)
    end do

   end if

  end do

end function total_water

! v6.9.4
function Get_reffc() ! liquid water
  real, dimension(nx,ny,nzm) :: Get_reffc
  Get_reffc = reffc
end function Get_reffc

function Get_reffi() ! ice
  real, dimension(nx,ny,nzm) :: Get_reffi
  Get_reffi = reffi
end function Get_reffi


end module microphysics



