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

tnudge = 0.
qnudge = 0.
unudge = 0.
vnudge = 0.

if(donudging_uv) then

    itauz(:) = 0.

    do k=1,nzm
      if(z(k).gt.nudging_uv_z1) then
        if(z(k).gt.nudging_uv_z2) then
          itauz(k) = coef
        elseif(z(k).gt.nudging_uv_z1.AND.z(k).le.nudging_uv_z2) then
          itauz(k) = coef*0.5*(1-cos(pii*(z(k)-nudging_uv_z1)/(nudging_uv_z2-nudging_uv_z1)))
        end if

        unudge(k)=unudge(k) - (u0(k)-ul0(k))*itauz(k)
        vnudge(k)=vnudge(k) - (v0(k)-vl0(k))*itauz(k)
        do j=1,ny
          do i=1,nx
             dudt(i,j,k,na)=dudt(i,j,k,na)-(u0(k)-ul0(k))*itauz(k)
             dvdt(i,j,k,na)=dvdt(i,j,k,na)-(v0(k)-vl0(k))*itauz(k)
          end do
        end do
      end if
    end do

endif

coef = 1./tautqls

if(donudging_tq.or.donudging_t) then
  ! --- Heng Xiao, 02/19/2024  
  ! Only nudge t, q between these two times
  if (time .gt. nudging_tq_t1 .and. time .lt. nudging_tq_t2) then
  ! --- Heng Xiao, 02/19/2024  
    coef1 = dtn / tautqls
    do k=1,nzm
      if(z(k).ge.nudging_t_z1.and.z(k).le.nudging_t_z2) then
        tnudge(k)=tnudge(k) -(t0(k)-tg0(k)-gamaz(k))*coef
        do j=1,ny
          do i=1,nx
             t(i,j,k)=t(i,j,k)-(t0(k)-tg0(k)-gamaz(k))*coef1
          end do
        end do
      end if
    end do
  ! --- Heng Xiao, 02/19/2024  
  ! Only nudge t, q between these two times
  end if
  ! --- Heng Xiao, 02/19/2024  
endif

if(donudging_tq.or.donudging_q) then
  ! --- Heng Xiao, 02/19/2024  
  ! Only nudge t, q between these two times
  if (time .gt. nudging_tq_t1 .and. time .lt. nudging_tq_t2) then
  ! --- Heng Xiao, 02/19/2024  
    coef1 = dtn / tautqls
    do k=1,nzm
      if(z(k).ge.nudging_q_z1.and.z(k).le.nudging_q_z2) then
        qnudge(k)=qnudge(k) -(q0(k)-qg0(k))*coef
        do j=1,ny
          do i=1,nx
             micro_field(i,j,k,index_water_vapor)=micro_field(i,j,k,index_water_vapor)-(q0(k)-qg0(k))*coef1
          end do
        end do
      end if
    end do
  ! --- Heng Xiao, 02/19/2024  
  ! Only nudge t, q between these two times
  end if
  ! --- Heng Xiao, 02/19/2024  
endif

end subroutine nudging
