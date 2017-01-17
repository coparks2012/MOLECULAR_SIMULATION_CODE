module mod_force_impro_amber
  use global
  use mod_minimg
  implicit none
  
  contains

    subroutine impro_amber_force
      implicit none
      integer :: step
      integer :: i1,i2,i3,i4,itype
      integer :: i,m
      integer :: T1,T2,clock_rate,clock_max
      double precision :: tolerance,small,v1x,v1y,v1z,v2x,v2y,v2z,v3x
      double precision :: v3y,v3z,ss1,ss2,ss3,r1,r2,r3,c0,c1,c2,s1,s2
      double precision :: s12,c,s,domega,a,a11,a22,a33,a12,a13,a23,sx1
      double precision :: sx2,sx12,sy1,sy2,sy12,sz1,sz2,sz12
      double precision :: vb2xm,vb2ym,vb2zm
      double precision :: vb2x,vb2y,vb2z
      double precision :: vb1x,vb1y,vb1z
      double precision :: vb3x,vb3y,vb3z
      double precision :: rb1,rb3,rb2
      double precision :: b1mag2,b1mag, b2mag2,b2mag, b3mag2,b3mag
      double precision :: ctmp,r12c1,r12c2
      double precision :: sb1,sb2,sb3
      double precision :: sc1,sc2
      double precision :: pd,rc2
      double precision :: c1mag,c2mag
      double precision :: f1(3),f2(3),f3(3),f4(3)
      
      
      call system_clock(T1,clock_rate,clock_max)
      
      tolerance = 0.05
      small = 0.001
      do i = 1,numimpro
         
         i1 = improlist(1,i)
         i2 = improlist(2,i)
         i3 = improlist(3,i)
         i4 = improlist(4,i)
         itype = improlist(5,i)
         
         !---1st bond
         vb1x = position(i1)%x - position(i2)%x
         vb1y = position(i1)%y - position(i2)%y
         vb1z = position(i1)%z - position(i2)%z
         call minimg(vb1x,vb1y,vb1z)
         
         vb2x = position(i3)%x - position(i2)%x
         vb2y = position(i3)%y - position(i2)%y
         vb2z = position(i3)%z - position(i2)%z
         call minimg(vb2x,vb2y,vb2z)
         
         vb2xm = -vb2x
         vb2ym = -vb2y
         vb2zm = -vb2z
         call minimg(vb2xm,vb2ym,vb2zm)


         !---3rd bond
         vb3x = position(i4)%x - position(i3)%x
         vb3y = position(i4)%y - position(i3)%y
         vb3z = position(i4)%z - position(i3)%z
         call minimg(vb3x,vb3y,vb3z)
         
         !---c0 calculation
         
         sb1 = 1.0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z)
         sb2 = 1.0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z)
         sb3 = 1.0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z)
         
         rb1 = sqrt(sb1)
         rb3 = sqrt(sb3)
         
         c0 = (vb1x*vb3x + vb1y*vb3y + vb1z*vb3z) * rb1*rb3
         
         !---1st and 2nd angle
         
         b1mag2 = vb1x*vb1x + vb1y*vb1y + vb1z*vb1z
         b1mag = sqrt(b1mag2)
         b2mag2 = vb2x*vb2x + vb2y*vb2y + vb2z*vb2z
         b2mag = sqrt(b2mag2)
         b3mag2 = vb3x*vb3x + vb3y*vb3y + vb3z*vb3z
         b3mag = sqrt(b3mag2)
         
         ctmp = vb1x*vb2x + vb1y*vb2y + vb1z*vb2z
         r12c1 = 1.0 / (b1mag*b2mag)
         c1mag = ctmp * r12c1
         
         ctmp = vb2xm*vb3x + vb2ym*vb3y + vb2zm*vb3z
         r12c2 = 1.0 / (b2mag*b3mag)
         c2mag = ctmp * r12c2
         
         !---cos and sin of 2 angles and final c
         
         sc1 = sqrt(1.0 - c1mag*c1mag)
         if (sc1 < SMALL) sc1 = SMALL 
         sc1 = 1.0/sc1
         
         sc2 = sqrt(1.0 - c2mag*c2mag)
         if (sc2 < SMALL) sc2 = SMALL
         sc2 = 1.0/sc2
         
         s1 = sc1 * sc1
         s2 = sc2 * sc2
         s12 = sc1 * sc2
         c = (c0 + c1mag*c2mag) * s12
         
         !---error check
         if (c.gt.1.0+tolerance.or.c.lt.(-1.0-tolerance)) then
            print*, 'Improper:',i1,i2,i3,i4
            print*, ' Error1:',position(i1)%x,position(i1)%y,position(i1)%z
            print*, ' Error2:',position(i2)%x,position(i2)%y,position(i2)%z
            print*, ' Error3:',position(i3)%x,position(i3)%y,position(i3)%z
            print*, ' Error4:',position(i4)%x,position(i4)%y,position(i4)%z
            
         endif
         
         
         if (c.gt.1.0) then
            c = 1.0
         endif
         if (c.lt.-1.0)then
            c = -1.0
         endif

         
         m = multiplicity_impro(itype)
         
         if(m.eq.2)then
            p = 2.0*c*c
            pd = 2.0*c
         elseif(m.eq.3)then
            rc2 = c*c;
            p = (4.0*rc2-3.0)*c + 1.0;
            pd = 6.0*rc2 - 1.5;
         elseif(m.eq.4)then
            rc2 = c*c;
            p = 8.0*(rc2-1)*rc2 + 2.0;
            pd = (16.0*rc2-8.0)*c;
         elseif(m.eq.6)then 
            rc2 = c*c;
            p = ((32.0*rc2-48.0)*rc2 + 18.0)*rc2
            pd = (96.0*(rc2-1.0)*rc2 + 18.0)*c
         elseif(m.eq.1)then 
            p = c + 1.0
            pd = 0.5
         elseif(m.eq.5)then
            rc2 = c*c
            p = ((16.0*rc2-20.0)*rc2 + 5.0)*c + 1.0
            pd = (40.0*rc2-30.0)*rc2 + 2.5
         elseif(m.eq.0)then
            p = 2.0
            pd = 0.0
         endif
         
         if(sign_impro(itype).eq.-1)then
            p = 2.0 - p
            pd = -pd
         endif
         
         e_improper = e_improper + k_impro(itype)*p
         

         a = 2.0 * k_impro(itype) * pd
         c = c * a
         s12 = s12 * a
         a11 = c*sb1*s1
         a22 = -sb2*(2.0*c0*s12 - c*(s1+s2))
         a33 = c*sb3*s2
         a12 = -r12c1*(c1mag*c*s1 + c2mag*s12)
         a13 = -rb1*rb3*s12
         a23 = r12c2*(c2mag*c*s2 + c1mag*s12)
         
         sx2  = a12*vb1x + a22*vb2x + a23*vb3x
         sy2  = a12*vb1y + a22*vb2y + a23*vb3y
         sz2  = a12*vb1z + a22*vb2z + a23*vb3z
         
         f1(1) = a11*vb1x + a12*vb2x + a13*vb3x
         f1(2) = a11*vb1y + a12*vb2y + a13*vb3y
         f1(3) = a11*vb1z + a12*vb2z + a13*vb3z
         
         f2(1) = -sx2 - f1(1)
         f2(2) = -sy2 - f1(2)
         f2(3) = -sz2 - f1(3)
         
         f4(1) = a13*vb1x + a23*vb2x + a33*vb3x
         f4(2) = a13*vb1y + a23*vb2y + a33*vb3y
         f4(3) = a13*vb1z + a23*vb2z + a33*vb3z
         
         f3(1) = sx2 - f4(1)
         f3(2) = sy2 - f4(2)
         f3(3) = sz2 - f4(3)
         
         !---apply force to each of 4 atoms
         ff_impro(i1)%x = f1(1) + ff_impro(i1)%x
         ff_impro(i1)%y = f1(2) + ff_impro(i1)%y
         ff_impro(i1)%z = f1(3) + ff_impro(i1)%z
         
         
         ff_impro(i2)%x = f2(1) + ff_impro(i2)%x
         ff_impro(i2)%y = f2(2) + ff_impro(i2)%y
         ff_impro(i2)%z = f2(3) + ff_impro(i2)%z
         
         
         ff_impro(i3)%x = f3(1) + ff_impro(i3)%x
         ff_impro(i3)%y = f3(2) + ff_impro(i3)%y
         ff_impro(i3)%z = f3(3) + ff_impro(i3)%z
         
         
         ff_impro(i4)%x = f4(1) + ff_impro(i4)%x
         ff_impro(i4)%y = f4(2) + ff_impro(i4)%y
         ff_impro(i4)%z = f4(3) + ff_impro(i4)%z
         
         
         v_improx = v_improx + vb1x*f1(1) + vb2x*f3(1) +&
              (vb3x+vb2x)*f4(1)
         v_improy = v_improy+ vb1y*f1(2) + vb2y*f3(2) +&
              (vb3y+vb2y)*f4(2)
         v_improz = v_improz +  vb1z*f1(3) + vb2z*f3(3) +&
              (vb3z+vb2z)*f4(3)
         
        

      end do
      call system_clock(T2,clock_rate,clock_max)
      time_impro = time_impro + real(T2-T1)/real(clock_rate)
      
    end subroutine impro_amber_force
    
  end module mod_force_impro_amber
  
