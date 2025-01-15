subroutine ice_nucl_constrain(nbin,TT,SUP2,ff2r,xi,ff3r,xs,ff4r,xg)

  use micro_prm

  implicit none
  integer  nbin

  real ff2r(nbin,3),ff3r(nbin),ff4r(nbin)
  real xi(nbin,3),xs(nbin),xg(nbin)

  real ni(nbin),ns(nbin),ng(nbin)

  reaL TT, SUP2,TPC,DEL2N, totice,newice,dx

  INTEGER ITYPE, ICE
  ! TYPE OF ICE WITH NUCLEATION (BEGIN)

  TPC=TT-273.15
  ITYPE=0

  IF((TPC.GT.-4.0) .OR. &
     (TPC.LE.-8.1 .AND. TPC.GT.-12.7) .OR.&
     (TPC.LE.-17.8 .AND. TPC.GT.-22.4)) THEN
    ITYPE=2
  ELSE
    IF((TPC.LE.-4.0 .AND. TPC.GT.-8.1) .OR. (TPC.LE.-22.4)) THEN
      ITYPE=1
    ELSE
      ITYPE=3
    ENDIF
  ENDIF

  ! NEW CRYSTAL SIZE DISTRIBUTION FUNCTION                      (BEGIN)

  ICE=ITYPE

  newice=0.0

  DEL2N=100.*SUP2
  IF(TPC.LT.0. .AND. TPC.GE.-35. .AND. DEL2N.GT.0.) THEN
    ni(:) = ff2r(:,1)*3.*xi(:,1)*col+ff2r(:,2)*3.*xi(:,2)*col &
          + ff2r(:,3)*3.*xi(:,3)*col 
    ns(:) = ff3r(:)*3.*xs(:)*col
    ng(:) = ff4r(:)*3.*xg(:)*col
    totice = sum(ni)+sum(ns)+sum(ng)

    ! totin = sum(fin)*col
    ! totin = 1.7e-3   ! 1.7 L-1
    ! totin = 1.7e-3*3.0   !  increase by a factor of 3
    ! totin = 1.7e-3*0.1   !  decrease by a factor of 10  
    if (totin .gt. totice .and. DEL2N .gt. 5.0) then
      newice = totin-totice
    else
      newice=0.0
    endif
       
    ! if (newice .gt. 1.7e-3) print*, 'newice', newice, ni,ns,ng

    ! add to the first ice bin
    dx = 3.*xi(1,ICE)*COL
    ff2r(1,ICE)=ff2r(1,ICE)+newice/dx
  endif

end subroutine ice_nucl_constrain 
