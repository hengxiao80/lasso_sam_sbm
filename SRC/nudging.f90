subroutine nudging()
	
use vars
use params
use microphysics, only: micro_field, index_water_vapor
implicit none

real coef, pii
integer i,j,k

real :: nudge_ramp(nzm)
real :: itau_transient, itau_transient_aloft, time_ramp
real :: itauz_t(nzm), itauz_q(nzm)
real :: zramp_transient(nzm), itau_transient_z(nzm)
real :: tmp_zinv_obs
real :: dist_zinv

real, external :: get_inversion_height

pii = acos(-1.)

tnudge = 0.
qnudge = 0.
unudge = 0.
vnudge = 0.

if(donudging_transient) then
  ! turn off nudging within 50m of observed inversion, so that model can adopt
  ! its own inversion structure if close to correct inversion height)
  if(transient_nudging_zinv .lt. 0.) then
    tmp_zinv_obs = get_inversion_height(nzm,z,pres,tg0+gamaz,tg0,qg0)
  else
    tmp_zinv_obs = transient_nudging_zinv
  end if

  if(masterproc .AND. (mod(nstep,nstat) .eq. 1) .AND. (icycle .eq. 1)) then
    write(*,935) transient_nudging_start, day, transient_nudging_end, tmp_zinv_obs
    935 format('transient nudging: starting day = ',F9.3,' current day = ',F9.3, &
                ' final day = ',F9.3,' obs_zinv = ',F9.2)
  end if

  !bloss: Define a vertically-uniform inverse nudging timescale
  !   that can be switched on/off during the simulation.
  !   The nudging of temperature and moisture will always 
  !   be at least as strong as itau_transient.
  if((day .gt. transient_nudging_start) .AND. (day .lt. transient_nudging_end)) then

    donudging_tq = .true.
    itau_transient = 1./tau_transient_nudging
    itau_transient_aloft = 1./600. ! ten minute nudging timescale above inversion
    
    !smooth start/finish to transient nudging
    time_ramp = &
      0.5 * (1. - cos(pi * &
      MAX(0., MIN(1., (day - transient_nudging_start)/transient_nudging_ramp)) &
      )) &
      * 0.5 * (1. - cos(pi * &
      MAX(0., MIN(1., (transient_nudging_end - day)/transient_nudging_ramp)) &
      )) 
    itau_transient = itau_transient * time_ramp

    zramp_transient(:) = 1.
    do k = 1,nzm
      dist_zinv = ABS(z(k)-tmp_zinv_obs) 
      IF (dist_zinv .lt. 50.) then
        zramp_transient(k) = 0.
      elseif (dist_zinv .lt. 100.) then
        zramp_transient(k) = 0.5*(1. - cos(pi*(dist_zinv-50.)/50.))
      end IF

      if(z(k) .lt. transient_nudging_zfloor) then
         zramp_transient(k) = zramp_transient(k) * &
           0.5 * (1. - cos(pi*z(k)/transient_nudging_zfloor))
      end if

      if(z(k) .lt. tmp_zinv_obs) then
        itau_transient_z(k) = itau_transient * zramp_transient(k)
      else
        ! different, shorter nudging time above inversion
        itau_transient_z(k) = itau_transient_aloft * zramp_transient(k)  
      end if
    end do

    if(icycle .eq. 1 .AND. mod(nstep,10) .eq. 1) then
      if(masterproc) write(*,*) 'Transient nudging on with itau(day) = ', itau_transient
    end if
  else
    itau_transient = 0.
    itau_transient_z(:) = 0.
  end if
end if

coef = 1./tauls
if(donudging_uv) then
    do k=1,nzm
      if(z(k).ge.nudging_uv_z1.and.z(k).le.nudging_uv_z2) then
        unudge(k)=unudge(k) - (u0(k)-ul0(k))*coef
        vnudge(k)=vnudge(k) - (v0(k)-vl0(k))*coef
        do j=1,ny
          do i=1,nx
            dudt(i,j,k,na)=dudt(i,j,k,na)-(u0(k)-ul0(k))*coef
            dvdt(i,j,k,na)=dvdt(i,j,k,na)-(v0(k)-vl0(k))*coef
          end do
        end do
      end if
    end do
endif

coef = 1./tautqls
if(donudging_tq .or. donudging_t .or. donudging_q) then
  ! Only nudge t, q between these two times
  if (time .gt. nudging_tq_t1 .and. time .lt. nudging_tq_t2) then

    ! set up nudging heights automatically by tracking inversion height
    if(dovariable_tauz) then
      nudging_t_z1 = get_inversion_height(nzm,z,pres,t0,tabs0,q0) &
          + variable_tauz_offset_above_inversion
      nudging_t_z1 = MAX(nudging_t_z1, variable_tauz_minimum_height)
      nudging_t_zramp = variable_tauz_thickness_of_onset

      nudging_q_z1 = nudging_t_z1
      nudging_q_zramp = nudging_t_zramp
    end if

    ! vertically-varying nudging amplitude for temperature
    nudge_ramp(:) = 0.
    do k=1,nzm
      if(z(k) .ge. nudging_t_z1 .and. z(k) .le. nudging_t_z2) then
        !nudging will be applied
        nudge_ramp(k) = 1.
        if(z(k) .lt. nudging_t_z1+nudging_t_zramp) then
          ! gradual onset of nudging between z1 and z1+zramp
          nudge_ramp(k) = 0.5*(1-cos(pii*(z(k)-nudging_t_z1)/nudging_t_zramp))
        elseif(z(k) .gt. nudging_t_z2 - nudging_t_zramp) then
          ! gradual falloff of nudging between z2 and z2-zramp
          nudge_ramp(k) = 0.5*(1-cos(pii*(nudging_t_z2 - z(k))/nudging_t_zramp))
        end if
      end if
    end do
    itauz_t(:) = coef*nudge_ramp(:)
    if(donudging_transient) then
      !bloss: Apply vertically-uniform nudging where stronger than standard nudging
      do k = 1,nzm
        itauz_t(k) = MAX(itauz_t(k), itau_transient_z(k))
      end do
    end if

    ! vertically-varying nudging amplitude for moisture
    nudge_ramp(:) = 0.
    do k=1,nzm
      if(z(k) .ge. nudging_q_z1 .and. z(k) .le. nudging_q_z2) then
        !nudging will be applied
        nudge_ramp(k) = 1.
        if(z(k) .lt. nudging_q_z1+nudging_q_zramp) then
          ! gradual onset of nudging between z1 and z1+zramp
          nudge_ramp(k) = 0.5*(1-cos(pii*(z(k)-nudging_q_z1)/nudging_q_zramp))
        elseif(z(k) .gt. nudging_q_z2 - nudging_q_zramp) then
          ! gradual falloff of nudging between z2 and z2-zramp
          nudge_ramp(k) = 0.5*(1-cos(pii*(nudging_q_z2 - z(k))/nudging_q_zramp))
        end if
      end if
    end do
    !bloss: Use itauz_q(1:nzm) to keep track of vertically-varying inverse nudging timescale
    itauz_q(:) = coef*nudge_ramp(:)
    if(donudging_transient) then
      !bloss: Apply vertically-uniform nudging where stronger than standard nudging
      do k = 1,nzm
        itauz_q(k) = MAX(itauz_q(k), itau_transient_z(k))
      end do
    end if

    if(donudging_tq .or. donudging_t) then
      do k=1,nzm
        !bloss: Apply nudging for all levels -- itauz_t(k) will be zero where no nudging is applied.
        tnudge(k)=tnudge(k) -(t0(k)-tg0(k)-gamaz(k))*itauz_t(k)
        do j=1,ny
          do i=1,nx
            t(i,j,k)=t(i,j,k)-(t0(k)-tg0(k)-gamaz(k))*dtn*itauz_t(k)
          end do
        end do
      end do
    endif

    if(donudging_tq .or. donudging_q) then
      do k=1,nzm
        !bloss: Apply nudging for all levels -- itauz_t(k) will be zero where no nudging is applied.
        qnudge(k)=qnudge(k) -(q0(k)-qg0(k))*itauz_q(k)
        do j=1,ny
          do i=1,nx
            micro_field(i,j,k,index_water_vapor) = micro_field(i,j,k,index_water_vapor) &
              - (q0(k) - qg0(k))*dtn*itauz_q(k)
          end do
        end do
      end do
    endif

  end if
endif

end subroutine nudging
