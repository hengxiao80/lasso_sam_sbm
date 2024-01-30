subroutine nudging()
	
use vars
use params
use microphysics, only: micro_field, index_water_vapor
implicit none

real coef, coef1, pii
integer i,j,k
real dqdt_qfloor
	
pii = acos(-1.)
coef = 1./tauls

!if(firststep) then
  tnudge = 0.
  qnudge = 0.
  unudge = 0.
  vnudge = 0.
!end if

if(donudging_uv) then

!  unudge = 0.
!  vnudge = 0.
    itauz(:) = 0.

    do k=1,nzm
      if(z(k).gt.tauzuv1) then
        if(z(k).gt.tauzuv2) then
          itauz(k) = 1./tauz_uv
        elseif(z(k).gt.tauzuv1.AND.z(k).le.tauzuv2) then
          itauz(k) = 1./tauz_uv*0.5*(1-cos(pii*(z(k)-tauzuv1)/(tauzuv2-tauzuv1)))
        end if

        unudge(k)=unudge(k) - (u0(k)-ug0(k))*itauz(k)
        vnudge(k)=vnudge(k) - (v0(k)-vg0(k))*itauz(k)
        do j=1,ny
          do i=1,nx
             dudt(i,j,k,na)=dudt(i,j,k,na)-(u0(k)-ug0(k))*itauz(k)
             dvdt(i,j,k,na)=dvdt(i,j,k,na)-(v0(k)-vg0(k))*itauz(k)
          end do
        end do
      end if
    end do

endif

coef = 1./tautqls

if(donudging_tq) then
    coef1 = dtn / tautqls
    do k=1,nzm
     if(z(k).ge.0.0 .and. z(k).le.tauz1) then
      tnudge(k)=tnudge(k)-(t0(k)-tg0(k)-gamaz(k))*coef
      do j=1,ny
       do i=1,nx
        t(i,j,k)=t(i,j,k)-(t0(k)-tg0(k)-gamaz(k))*coef1
       end do
      end do
     end if
    end do
end if

if(donudging_tq) then
    coef1 = dtn / tautqls
    do k=1,nzm
     if(z(k).ge.0.0 .and. z(k).le.tauz1) then
          qnudge(k)=qnudge(k) -(q0(k)-qg0(k))*coef
      do j=1,ny
       do i=1,nx
              micro_field(i,j,k,index_water_vapor) &
                   = micro_field(i,j,k,index_water_vapor) &
                    - (q0(k)-qg0(k))*coef1
       end do
      end do
     end if
    end do
end if

end subroutine nudging
