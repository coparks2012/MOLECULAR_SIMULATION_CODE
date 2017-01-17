module mod_force_dihed_amber
  use global
  use mod_minimg
  implicit none
  contains

    subroutine dihed_PME_amber
      implicit none
      integer :: i1,i2,i3,i4,n,type
      integer :: itype,jtype
      integer :: step
      integer :: i
      integer :: offset
      integer :: T1,T2,clock_rate,clock_max
      double precision :: vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb2xm,vb2ym,vb2zm
      double precision :: f1(3),f2(3),f3(3),f4(3)
      double precision :: sb1,sb2,sb3,rb1,rb3,c0,b1mag2,b1mag,b2mag2
      double precision :: b2mag,b3mag2,b3mag,ctmp,r12c1,c1mag,r12c2
      double precision :: c2mag,sc1,sc2,s1,s12,c,p,pd,a,a11,a22
      double precision :: a33,a12,a13,a23,sx2,sy2,sz2
      double precision :: s2,cx,cy,cz,cmag,dx,phi,si,siinv,sin2
      double precision :: k1,k2,k3,k4
      double precision :: tolerance,SMALL,SMALLER
      double precision :: ddx,ddy,ddz
      double precision :: force,forcelj
      double precision :: dr2,dr2i,dr6i,dr12i
      double precision :: qtmp,r,prefactor,erfcc,erfcd,t
      double precision :: gbb,dtfx,dtfy,dtfz,dtgx,dtgy,dtgz,dthx
      double precision :: dthy,dthz,df
      double precision :: forcecoul
      
      !---new terms
      double precision :: ax,ay,az,bx,by,bz,rasq,rbsq,rgsq
      double precision :: rg,rginv,rb2inv,ra2inv
      double precision :: rabinv,s,m,ddf1,df1,fg,hg,fga,hgb,gaa
      
      double precision :: q1,q2
      
      call system_clock(T1,clock_rate,clock_max)
      
      e_dihedral = 0.0d0
      tolerance = 0.001
      SMALL     = 0.001
      SMALLER   = 0.0001  
      
      do n=1,numdihed
         
         i1 = dihedlist(1,n)
         i2 = dihedlist(2,n)
         i3 = dihedlist(3,n)
         i4 = dihedlist(4,n)
         type = dihedlist(5,n)
         
         !--- 1st bond
         vb1x = position(i1)%x - position(i2)%x
         vb1y = position(i1)%y - position(i2)%y
         vb1z = position(i1)%z - position(i2)%z
         call minimg(vb1x,vb1y,vb1z)
         
         
         !---2nd bond
         vb2x =  position(i3)%x - position(i2)%x
         vb2y =  position(i3)%y - position(i2)%y
         vb2z =  position(i3)%z - position(i2)%z
         call minimg(vb2x,vb2y,vb2z)
         
         vb2xm = -vb2x
         vb2ym = -vb2y
         vb2zm = -vb2z
         call minimg(vb2xm,vb2ym,vb2zm)
         
         
         ! // 3rd bond
         
         vb3x = position(i4)%x - position(i3)%x
         vb3y = position(i4)%y - position(i3)%y
         vb3z = position(i4)%z - position(i3)%z
         call minimg(vb3x,vb3y,vb3z)


         !--- 1-4 neighbors
         ddx = position(i1)%x-position(i4)%x
         ddy = position(i1)%y-position(i4)%y
         ddz = position(i1)%z-position(i4)%z
         call minimg(ddx,ddy,ddz)
         
         dr2 = ddx*ddx + ddy*ddy + ddz*ddz
         itype = position(i1)%type
         jtype = position(i4)%type


         
         ax = vb1y*vb2zm - vb1z*vb2ym
         ay = vb1z*vb2xm - vb1x*vb2zm
         az = vb1x*vb2ym - vb1y*vb2xm
         bx = vb3y*vb2zm - vb3z*vb2ym
         by = vb3z*vb2xm - vb3x*vb2zm
         bz = vb3x*vb2ym - vb3y*vb2xm
         
         rasq = ax*ax + ay*ay + az*az
         rbsq = bx*bx + by*by + bz*bz
         rgsq = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm
         rg = sqrt(rgsq)
         
         rginv = 0.0; ra2inv = 0.0; rb2inv = 0.0
         if (rg.gt.0)then
            rginv = 1.0/rg
         endif
         if (rasq.gt.0)then
            ra2inv = 1.0/rasq
         endif
         if (rbsq.gt.0)then 
            rb2inv = 1.0/rbsq
         endif

         rabinv = sqrt(ra2inv*rb2inv)
         
         c = (ax*bx + ay*by + az*bz)*rabinv
         s = rg*rabinv*(ax*vb3x + ay*vb3y + az*vb3z)
         
         !--- error check
         !    // error check
         if (c.gt.1.0 + TOLERANCE.or.c.lt.(-1.0 - TOLERANCE)) then
            print*,'c and tolerance:',c,TOLERANCE
            print*,'1ST ATOM:',position(i1)%x,position(i1)%y,position(i1)%z
            print*,'2ND ATOM:',position(i2)%x,position(i2)%y,position(i2)%z
            print*,'3RD ATOM:',position(i3)%x,position(i3)%y,position(i3)%z
            print*,'4TH ATOM:',position(i4)%x,position(i4)%y,position(i4)%z
         endif
         
         
         if (c.gt.1.0) c = 1.0
         if (c .lt.-1.0) c = -1.0
         
         !======
         !-----error starts here
         m = multiplicity(type)
         p = 1.0
         ddf1 = 0.0
         df1  = 0.0
         
         do i = 1,m
            ddf1 = p*c - df1*s
            df1 = p*s + df1*c
            p = ddf1
         enddo
         
         
         p = p*cos_shift(type)     + df1*sin_shift(type)
         df1 = df1*cos_shift(type) - ddf1*sin_shift(type)
         df1 = -m*df1
         p = p +1.0
         
         if (m.eq.0) then
            p = 1.0 + cos_shift(type)
            df1 = 0.0
         endif
         !==========
         !---error ends here
         
         e_dihedral = e_dihedral+ k_dihed(type) * p
         
         fg = vb1x*vb2xm + vb1y*vb2ym + vb1z*vb2zm
         hg = vb3x*vb2xm + vb3y*vb2ym + vb3z*vb2zm
         fga = fg*ra2inv*rginv
         hgb = hg*rb2inv*rginv
         gaa = -ra2inv*rg
         gbb = rb2inv*rg
         
         dtfx = gaa*ax
         dtfy = gaa*ay
         dtfz = gaa*az
         dtgx = fga*ax - hgb*bx
         dtgy = fga*ay - hgb*by
         dtgz = fga*az - hgb*bz
         dthx = gbb*bx
         dthy = gbb*by
         dthz = gbb*bz
         
         df = -k_dihed(type) * df1
         
         sx2 = df*dtgx
         sy2 = df*dtgy
         sz2 = df*dtgz
         
         f1(1) = df*dtfx
         f1(2) = df*dtfy
         f1(3) = df*dtfz
         
         f2(1) = sx2 - f1(1)
         f2(2) = sy2 - f1(2)
         f2(3) = sz2 - f1(3)
         
         f4(1) = df*dthx
         f4(2) = df*dthy
         f4(3) = df*dthz
         
         f3(1) = -sx2 - f4(1)
         f3(2) = -sy2 - f4(2)
         f3(3) = -sz2 - f4(3)        
         
              !========time to calculate 1-4 PME nonbonded forces
         force = 0.0d0        
         forcelj = 0.0d0
         if(dr2.lt.rcut2)then
            
            
            dr2i     =  1.0d0/dr2
            dr6i     =  dr2i*dr2i*dr2i 
            forcelj  =  dr6i*(lj14_1(itype,jtype)*dr6i-lj14_2(itype,jtype))             
            pot_14   =  pot_14 + dr6i*(dr6i*lj14_3(itype,jtype)-lj14_4(itype,jtype))
            force    = forcelj*dr2i
            
         end if
         
  
         v_dihedx = v_dihedx  + ddx*ddx*force
         v_dihedy = v_dihedy  + ddy*ddy*force
         v_dihedz = v_dihedz  + ddz*ddz*force
         
         
         
         !...// apply force to each of 4 atoms
         ff_dihed(i1)%x =  ff_dihed(i1)%x + f1(1) + ddx*force
         ff_dihed(i1)%y =  ff_dihed(i1)%y + f1(2) + ddy*force
         ff_dihed(i1)%z =  ff_dihed(i1)%z + f1(3) + ddz*force
         
         ff_dihed(i2)%x =  ff_dihed(i2)%x + f2(1)
         ff_dihed(i2)%y =  ff_dihed(i2)%y + f2(2)
         ff_dihed(i2)%z =  ff_dihed(i2)%z + f2(3)
         
         ff_dihed(i3)%x =  ff_dihed(i3)%x + f3(1)
         ff_dihed(i3)%y =  ff_dihed(i3)%y + f3(2)
         ff_dihed(i3)%z =  ff_dihed(i3)%z + f3(3)
         
         ff_dihed(i4)%x =  ff_dihed(i4)%x + f4(1) - ddx*force
         ff_dihed(i4)%y =  ff_dihed(i4)%y + f4(2) - ddy*force
         ff_dihed(i4)%z =  ff_dihed(i4)%z + f4(3) - ddz*force
         
         
         
         v_dihedx = v_dihedx + vb1x*f1(1) + vb2x*f3(1) + (vb3x+vb2x)*f4(1)
         v_dihedy = v_dihedy + vb1y*f1(2) + vb2y*f3(2) + (vb3y+vb2y)*f4(2)
         v_dihedz = v_dihedz + vb1z*f1(3) + vb2z*f3(3) + (vb3z+vb2z)*f4(3)
         
         
      end do
      
      call system_clock(T2,clock_rate,clock_max)
      time_dihed = time_dihed + real(T2-T1)/real(clock_rate)
   
    end subroutine dihed_PME_amber


  end module mod_force_dihed_amber
