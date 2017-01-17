module mod_force_bond
  use global
  use mod_minimg
  contains

  subroutine bond_harmonic(step)
    implicit none
    integer :: step
    integer :: i
    integer :: m,i1,i2,itype,ifactor
    integer :: T1,T2,clock_rate,clock_max
    double precision:: delx,dely,delz,rsq,r,dr,rk,force
    
    
    call system_clock(T1,clock_rate,clock_max)

    
!    !$omp  parallel do reduction(+:e_bond) default(firstprivate),&
!    !$omp& shared(position,bondlist,bondcoeff,ff)
    do m = 1,numbond
       
       i1 = bondlist(1,m)
       i2 = bondlist(2,m)
       itype = bondlist(3,m)
    
       
       delx = position(i1)%x - position(i2)%x
       dely = position(i1)%y - position(i2)%y
       delz = position(i1)%z - position(i2)%z
       call minimg(delx,dely,delz)
       
       rsq = delx*delx + dely*dely + delz*delz
       r = sqrt(rsq)
       
       
       
       dr = r - bondcoeff(2,itype)
       rk = bondcoeff(1,itype) * dr
       e_bond = e_bond + rk*dr


       
       if (r.gt.0.0) then
          force = -2.0*rk/r
       else
          force = 0.0
       endif
       

       ff_bond(i1)%x = ff_bond(i1)%x + delx*force
       ff_bond(i1)%y = ff_bond(i1)%y + dely*force
       ff_bond(i1)%z = ff_bond(i1)%z + delz*force
       
       ff_bond(i2)%x = ff_bond(i2)%x - delx*force
       ff_bond(i2)%y = ff_bond(i2)%y - dely*force
       ff_bond(i2)%z = ff_bond(i2)%z - delz*force

       v_bondx = v_bondx  + delx*delx*force
       v_bondy = v_bondy  + dely*dely*force
       v_bondz = v_bondz  + delz*delz*force


       !ffh(i1)%x = ffh(i1)%x + delx*force
       !ffh(i1)%y = ffh(i1)%y + dely*force
       !ffh(i1)%z = ffh(i1)%z + delz*force
       
       !ffh(i2)%x = ffh(i2)%x - delx*force
       !ffh(i2)%y = ffh(i2)%y - dely*force
       !ffh(i2)%z = ffh(i2)%z - delz*force
    enddo
!    !$omp end parallel do
    call system_clock(T2,clock_rate,clock_max)
    time_bond = time_bond + real(T2-T1)/real(clock_rate)

  end subroutine bond_harmonic
  


end module mod_force_bond
