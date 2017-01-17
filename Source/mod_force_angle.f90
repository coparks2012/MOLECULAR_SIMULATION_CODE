module mod_force_angle
  use global
  use mod_minimg

  contains

  
  

  subroutine angle_harmonic(step)
    implicit none
    integer :: step
    integer :: m,i1,i2,i3,itype,ifactor
    integer :: i
    integer :: T1,T2,clock_rate,clock_max
    double precision :: small,delx1,dely1,delz1,rsq1,r1,delx2,dely2
    double precision ::  delz2,rsq2,r2,c,s,dtheta,tk,a,a11,a12,a22
    double precision ::  vx1,vx2,vy1,vy2,vz1,vz2

    small = 0.001

    call system_clock(T1,clock_rate,clock_max)
!    !$omp parallel do reduction(+:e_angle) default(firstprivate),&
!    !$omp& shared(anglelist,position,anglecoeff,ff)
   
    do m = 1,numangle
       
       i1 = anglelist(1,m)
       i2 = anglelist(2,m)
       i3 = anglelist(3,m)
       itype = anglelist(4,m)


       !---first bond
       delx1 = position(i1)%x - position(i2)%x
       dely1 = position(i1)%y - position(i2)%y
       delz1 = position(i1)%z - position(i2)%z
       call minimg(delx1,dely1,delz1)
       
       rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1
       r1 = sqrt(rsq1)
       
       !--- 2nd bond
       
       delx2 = position(i3)%x - position(i2)%x
       dely2 = position(i3)%y - position(i2)%y
       delz2 = position(i3)%z - position(i2)%z
       call minimg(delx2,dely2,delz2)
       
       rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2
       r2 = sqrt(rsq2)
       
       
       c = delx1*delx2 + dely1*dely2 + delz1*delz2
       c = c / (r1*r2)
       
       if (c.gt.1.0) c = 1.0
       if (c.lt.-1.0) c = -1.0
       
       s = sqrt(1.0 - c*c)
       if (s.lt.small) s = small
       s = 1.0/s
       
       
       
       dtheta = acos(c) - anglecoeff(2,itype)
       tk = anglecoeff(1,itype) * dtheta
       e_angle = e_angle + tk*dtheta
   

       a = 2.0 * tk * s
       
       a11 = a*c / rsq1
       a12 = -a / (r1*r2)
       a22 = a*c / rsq2
       
       vx1 = a11*delx1 + a12*delx2
       vx2 = a22*delx2 + a12*delx1
       vy1 = a11*dely1 + a12*dely2
       vy2 = a22*dely2 + a12*dely1
       vz1 = a11*delz1 + a12*delz2
       vz2 = a22*delz2 + a12*delz1
       



       ff_angle(i1)%x = ff_angle(i1)%x - vx1
       ff_angle(i1)%y = ff_angle(i1)%y - vy1
       ff_angle(i1)%z = ff_angle(i1)%z - vz1
       
       ff_angle(i2)%x = ff_angle(i2)%x + (vx2+vx1)
       ff_angle(i2)%y = ff_angle(i2)%y + (vy2+vy1)
       ff_angle(i2)%z = ff_angle(i2)%z + (vz2+vz1)
       
       ff_angle(i3)%x = ff_angle(i3)%x - vx2
       ff_angle(i3)%y = ff_angle(i3)%y - vy2
       ff_angle(i3)%z = ff_angle(i3)%z - vz2

       v_anglex = v_anglex - (delx1*vx1 + delx2*vx2)
       v_angley = v_angley - (dely1*vy1 + dely2*vy2)
       v_anglez = v_anglez - (delz1*vz1 + delz2*vz2)
   
       !ffh(i1)%x = ffh(i1)%x - vx1
       !ffh(i1)%y = ffh(i1)%y - vy1
       !ffh(i1)%z = ffh(i1)%z - vz1
       
       !ffh(i2)%x = ffh(i2)%x + (vx2+vx1)
       !ffh(i2)%y = ffh(i2)%y + (vy2+vy1)
       !ffh(i2)%z = ffh(i2)%z + (vz2+vz1)
       
       !ffh(i3)%x = ffh(i3)%x - vx2
       !ffh(i3)%y = ffh(i3)%y - vy2
       !ffh(i3)%z = ffh(i3)%z - vz2

    enddo
!    !$omp end parallel do

    call system_clock(T2,clock_rate,clock_max)
    time_angle = time_angle + real(T2-T1)/real(clock_rate)

    
  end subroutine angle_harmonic


end module mod_force_angle
