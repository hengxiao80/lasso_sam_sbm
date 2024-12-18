module micro_prm

implicit none

!-------------------------------------------------------------------!
! SAM-SBM  by Jiwen Fan                                             !
!-------------------------------------------------------------------!

! number of bins for CCN
integer, parameter :: ncn=33, nkr=33
! add size bins for IN
integer, parameter :: nin=33

! number of bins for droplets
integer, parameter :: ncd=33

!total number of prognostic microphysical variables
integer, parameter :: nmicro_fields = 1+ncn+ncd*7+nin

! for bulk nucleation
integer, parameter :: BULKNUC=0

! COL=Ln2/3 (for xl(k+1)=xl(k)*2.)

real, parameter :: COL=0.23105

! Aerosols
!  add for diagnostic CCN. If diagCCN=.false., CCN is prognostic
! LOGICAL, PARAMETER :: diagCCN=.false.
LOGICAL, PARAMETER :: diagCCN=.True.

! Parameters are used for calculation of 

! RO_SOLUTE - density of aerosol+h20

!        PARAMETER (RO_SOLUTE=2.16)
real, parameter :: RO_SOLUTE=1.79

! RO_SOLUTE - density of ammo sulfate + organic
!        PARAMETER (RO_SOLUTE=(1.8D0*0.6D0+2.0D0*0.4D0))

! Molecular weight of dry aerosol
!        PARAMETER (mwaero = (22.9+35.5))
real, parameter :: mwaero = 115.0
! Number of ions
integer, parameter :: ions = 3
! TEMPERATURA IN SURFACE LAYER EQUAL 15 Celsius(288.15 K)
real, parameter :: tsfc_aero=288.15

! ROCCN0, g/cm^3 - density of aerosol

real, parameter :: ROCCN0=1.0E00, ROCCN02=0.1000E01, ROCCN03=0.1000E01

! ACCN in 1/cm^3, BCCN - coefficients for calculation of FCCNR(KR)
! FCCNR(KR)=1.5*ACCN*BCCN*S_KR**BCCN,
! where S_KR - water supersaturation, %

real, parameter :: ACCN=1.0000E02, BCCN=0.4620E00
real, parameter :: ACCN2=1.25000E03, BCCN2=0.3080E00
real, parameter :: ACCN3=1.0000E02, BCCN3=0.4620E00

! ICCN - set maximum mass of aerosol
! XL(ICCN) - maximum mass of aerosol, is determined as minimum
! mass drop in bin with number ICCN
                                                                                                             
integer, parameter :: ICCN=1

!c DELZCCN - changing with heihgt

!c if DELZCCN=0 then no vertical distribution concentration of aerosol

!c        PARAMETER (DELZCCN=0.0E0)
real, parameter :: DELZCCN=1.2E05

!c Z0CCN=1KM=1.E5CM

real, parameter :: Z0CCN=2.E5

real, parameter :: Z0IN=2.5E5

!c flags for turbulance
integer, parameter :: LIQTURB=1
integer, parameter :: ICETURB=1

!c parameter used for calculation of ice microphysics

!c if ICEPROCS=1 it is ice microphysics

integer, parameter :: ICEPROCS=1

! ICEFLAG=0 USE MEYERS ICE NUCLEATION SCHEME; ICEFLAG=1 USE Classical Theory
integer, parameter :: ICEFLAG=0

!c ICEMAX - number of ice crystal types

integer, parameter :: ICEMAX=3

!c FREEZING

integer, parameter :: NHYDR=5,NHYDRO=7                                 &
     &        ,ifreez_down1=0,ifreez_down2=1,ifreez_top=1              &
     &        ,K0_LL=8,KRMIN_LL=1,KRMAX_LL=19,L0_LL=6                  &
     &        , IEPS_400=1,IEPS_800=0,IEPS_1600=0                      &
     &        ,K0L_GL=16,K0G_GL=16                                     &
     &        ,KRMINL_GL=1,KRMAXL_GL=24                                &
     &        ,KRMING_GL=1,KRMAXG_GL=33                                &
     &        ,KRDROP=18,KRBREAK=17                                    &
     &        ,JMAX=33,NRG=2,JBREAK = 18 
                                                                                                             
!real, parameter :: C1_MEY=0.00033,C2_MEY=0.
real, parameter :: C1_MEY=0.00012,C2_MEY=0.

! For drop freezing (iceform=1) or drop evaporation nuclei (iceform = 2) mechanism 
integer, parameter :: iceform = 1 

! For the diagnostic ice formation from SHEBA intercomparion (fixice = 1 and iceprocs=1), totin is used in ice_nucl_constrain.f90
integer, parameter :: fixice = 0  ! this should be 0 if other ice nucleation mechanism is used
!real, parameter :: totin = 1.7e-3   ! 1.7 L-1
!real, parameter :: totin = 1.7e-3*3.
!real, parameter :: totin = 1.7e-3*0.1
!real, parameter :: totin = 0.4e-3    ! 0.4 L-1 for ISDAC flt31
!real, parameter :: totin = 0.5e-3    ! 0.5 L-1 for ISDAC flt31
!real, parameter :: totin = 0.    ! No ICE
real, parameter :: totin = 1.0e-3    ! 1.0 L-1 for sensitivity ISDAC flt31
!real, parameter :: totin = 4.0e-3    ! 4.0 L-1 for sensitivity
!real, parameter :: totin = 10.0e-3   ! 10.0 L-1 for sensitivity
!real, parameter :: totin = 1.5e-3    ! 1.5 L-1 for sensitivity

!c if ifreeze_top=1 homogeneous freezing at the top is included
!c if ifreeze_top=0 no homogeneous freezing at the top

!c	PARAMETER (IFREEZE_TOP=1)

!c Parameters used into algorithm of heterogenic freezing by Vali,1975:
!c an0_freeze=10/cm3(Vali,1975), gamma=4.4

real, parameter :: an0_freez=10.

!c        PARAMETER (gamma=4.4)

!c KRFREEZ set bin boundary between plates and hail with freezing

!        PARAMETER (KRFREEZ=21)

!c AFREEZE, BFREEZE, BFREEZEMAX - coefficients in formula of 
!c                                homogeneous freezing
!        PARAMETER (BFREEZMAX=0.660)

integer, PARAMETER :: KRFREEZ=21
REAL, PARAMETER :: BFREEZMAX=0.66E0
!c other parameters and thresholds 

real, parameter :: ANC1=1.0000E02, ANC2=1.0000E02, PARSP=0.9000E02, PEPL=0.6000E00,     &
     &  PEPI=0.6000E00,PERL=1.0000E-03,PERI=1.0000E-03,TAUR1=0.5000E00,TAUR2=0.8000E00, &
     &  EPS_R1=0.1500E09,EPS_R2=0.1500E09,AFREEZMY=0.3333E-04,         &
     &  BFREEZMY=0.6600E00,EPSF1=0.1000E-02,EPSF3=0.1000E-05,EPSF4=0.1000E-05

!c MELTING

!c DIFFUSIONAL GROWTH

!c Parameters are used in algorithm of diffusional growth

!c NCOND determine timestep (DT/NCOND) with diffusional growth

integer, parameter :: NCOND=1

!c Coeffients for diffusional growth
real, parameter :: A1_MYN=2.53,BB1_MYN=5.42, A2_MYN=3.41E1,BB2_MYN=6.13
real, parameter :: AA1_MY=2.53E12,BB1_MY=5.42E3, AA2_MY=3.41E13,BB2_MY=6.13E3 


!c COAGULATION

!c Kernels (collisions), depend on heights

!c P1=1000mb=1000000 dynes/cm^2

real, parameter :: p1=1000000.0,p2=750000.0,p3=500000.0

!c Parameters used into subroutine coal_bott (collisions):

!c if temperature less than TTCOAL then no calculations of
!c temperature dependence of collision kernels

real, parameter :: TTCOAL=233.15E00

!c TCRIT - temperature boundary between graupel and hail with
!c coagulation

real, parameter :: TCRIT=270.15E00

!c if i_ww=1 then calculates only for water-water=water

integer, parameter :: I_WW=0

!c alcr=1.0 g/m**3 :threshold for graupel formation snow-drop collisions

real, parameter :: ALCR=1.5

!c  threshold for hail formation graupel-drop collisions
integer, parameter :: alcr_hail=3

!c threshold bin number for hail formation from graupel-drop collisions
integer, parameter :: kp_hail=25


!c scal=1.0 is valid under condition of mass doubling : 
!c xl(k+1)=xl(k)*2 

real, parameter :: SCAL=1.0E00

!c Parameters used for ice multiplication :

!c if icempl=0 no ice multiplication; 
!c if icempl=1 - ice multiplication is included

integer, parameter :: ICEMPL=1
integer, parameter :: kr_icempl=9

!c 3-point remapping
integer, parameter :: I3POINT=1

!--- Constants specific to the parameterization follow:
!--- CLIMIT/CLIMIT1 are lower limits for treating accumulated precipitation
!
! ======================================================================
!--- Important tunable parameters that are exported to other modules
!  * RHgrd - threshold relative humidity for onset of condensation
!  * T_ICE - temperature (C) threshold at which all remaining liquid water
!            is glaciated to ice
!  * T_ICE_init - maximum temperature (C) at which ice nucleation occurs
!  * NLImax - maximum number concentrations (m**-3) of large ice (snow/graupel/sleet)
!  * NLImin - minimum number concentrations (m**-3) of large ice (snow/graupel/sleet)
!  * N0r0 - assumed intercept (m**-4) of rain drops if drop diameters are between 0.2 and 0.45 mm
!  * N0rmin - minimum intercept (m**-4) for rain drops
!  * NCW - number concentrations of cloud droplets (m**-3)
!  * FLARGE1, FLARGE2 - number fraction of large ice to total (large+snow) ice
!          at T>0C and in presence of sublimation (FLARGE1), otherwise in
!          presence of ice saturated/supersaturated conditions
! ======================================================================
real, parameter ::                                                      &
     &  RHgrd=1.0                                                       &
     & ,T_ICE=-30.                                                      &
     & ,T_ICE_init=-5.                                                  &
     & ,NLImax=20.E3                                                    &
     & ,NLImin=100.                                                     &
     & ,N0r0=8.E6                                                       &
     & ,N0rmin=1.E4                                                     &
     & ,NCW=200.E6                                                      &
     & ,FLARGE1=1.                                                      &
     & ,FLARGE2=.2 

! For new ice nucleation from J. Comstock (ICEFLAG =1)

!       ====================================================================
!       aerosol_prop_mks.inc
!       version
!               adapted from aerosol_prop_mks.h; from J.Comstock
!       purpose
!               -define properties for aerosol composition used in
!               heterogeneous nucleation (KC1999)
!               -define aerosol size distribution parameters
!       variables
!               total aerosol number concentration, prior to nucleation (m-3)
!               smallest aerosol radius in size distribution (m)
!               aerosol solubility in terms of volume fraction (Qv)
!               aerosol parameter describing composition
!               molecular weight of dry aerosol (g/mol) (NH4)2SO4 Ammonium Sulate KC1999
!               density of dry aerosol (kg/m3) (NH4)2SO4 Ammonium Sulate KC1999
!       description
                                                                                                           
!       NOTE MKS units
                                                                                                           
!       =========================================================================
                                                                                                           
!         PARAMETER (na=10.e3)     !total aerosol number concentration, prior to nucleation (m-3) org=500.e4
!real, parameter :: rmin= 0.03e-6  !smallest aerosol radius in size distribution (m)
real, parameter :: qvaero=0.5    !aerosol solubility in terms of volume fraction (Qv)
real, parameter :: betafr=0.5     !aerosol parameter describing composition
real, parameter :: Ms2=0.132      !molecular weight of dry aerosol (g/mol) (NH4)2SO4 Ammonium Sulfate KC1999
real, parameter :: rhos2=1770.0   !density of dry aerosol (kg/m3) (NH4)2SO4 Ammonium Sulate KC1999
real, parameter ::  wetcoef = 0.90 !wettability
real, parameter :: alf = 0.1e-5   ! relative area of active sites
real, parameter :: epsil = 2.5e-2 ! misfit strain parameter
real,parameter :: fracin = 0.008e-3 ! IN fraction
!         PARAMETER (mua=2.5)        !slope of junge aerosol size distribution
!real, parameter :: drdry=0.01e-6  !0.01 microns; width of aerosol size bin (dry)
                                                                                                           
!       ==========================================================================
!!! For ccn regeneration
real :: ccnreg
real :: inreg
!!!

end module micro_prm
