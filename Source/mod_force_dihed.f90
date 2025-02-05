module mod_force_dihed
  use global
  use mod_minimg
  contains

 subroutine dihedral_opls_n2(step)
    implicit none
    integer :: i1,i2,i3,i4,n,type
    integer :: itype,jtype
    integer :: step
    integer :: offset
    integer :: i
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
    double precision :: lj11,lj12,lj13,lj21,lj22,lj23

    e_dihedral = 0.0
    tolerance = 0.001
    SMALL     = 0.001
    SMALLER   = 0.0001  

 
    
    do n=1,numdihed

       i1 = dihedlist(1,n)
       i2 = dihedlist(2,n)
       i3 = dihedlist(3,n)
       i4 = dihedlist(4,n)
       type = dihedlist(5,n)
      
      ! // 1st bond
       
       vb1x = position(i1)%x - position(i2)%x
       vb1y = position(i1)%y - position(i2)%y
       vb1z = position(i1)%z - position(i2)%z
       call minimg(vb1x,vb1y,vb1z)

       
      ! // 2nd bond

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

      
       
      !// c0 calculation
       
       sb1 = 1.0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z)
       sb2 = 1.0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z)
       sb3 = 1.0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z)
       
       rb1 = sqrt(sb1)
       rb3 = sqrt(sb3)
       
       c0 = (vb1x*vb3x + vb1y*vb3y + vb1z*vb3z) * rb1*rb3
       
       !// 1st and 2nd angle
       
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
       
       !// cos and sin of 2 angles and final c
       
       sin2 = MAX(1.0 - c1mag*c1mag,0.0)
       sc1 = sqrt(sin2)
       if (sc1 < SMALL) sc1 = SMALL
       sc1 = 1.0/sc1
       
       sin2 = MAX(1.0 - c2mag*c2mag,0.0)
       sc2 = sqrt(sin2)
       if (sc2 < SMALL) sc2 = SMALL
       sc2 = 1.0/sc2
       
       s1 = sc1 * sc1
       s2 = sc2 * sc2
       s12 = sc1 * sc2
       c = (c0 + c1mag*c2mag) * s12
       
  !     cx = vb1y*vb2z - vb1z*vb2y
  !     cy = vb1z*vb2x - vb1x*vb2z
  !     cz = vb1x*vb2y - vb1y*vb2x
  !     cmag = sqrt(cx*cx + cy*cy + cz*cz)
  !     dx = (cx*vb3x + cy*vb3y + cz*vb3z)/cmag/b3mag
       
   !    // error check
       if (c.gt.1.0 + TOLERANCE.or.c.lt.(-1.0 - TOLERANCE)) then
          print*,'c and tolerance:',c,TOLERANCE
          print*,'1ST ATOM:',position(i1)%x,position(i1)%y,position(i1)%z
          print*,'2ND ATOM:',position(i2)%x,position(i2)%y,position(i2)%z
          print*,'3RD ATOM:',position(i3)%x,position(i3)%y,position(i3)%z
          print*,'4TH ATOM:',position(i4)%x,position(i4)%y,position(i4)%z
       endif
          
       
       if (c.gt.1.0) c = 1.0
       if (c.lt.-1.0) c = -1.0
       
     
       phi = acos(c)
    !   if (dx.lt.0.0) phi = -1.0d0
       si = sin(phi)
       if (abs(si).lt.SMALLER) si = SMALLER
       siinv = 1.0/si

       k1 = dihedcoeff(1,type); k2 = dihedcoeff(2,type); k3 = dihedcoeff(3,type); k4 = dihedcoeff(4,type)
       
       p = k1*(1.0d0 + c) + k2*(1.0d0 - cos(2.0d0*phi)) + k3*(1.0d0 + cos(3.0d0*phi))+k4*(1.0d0 - cos(4.0d0*phi)) 

       pd = k1 - 2.0d0*k2*sin(2.0d0*phi)*siinv +3.0d0*k3*sin(3.0d0*phi)*siinv-4.0d0*k4*sin(4.0d0*phi)*siinv

       e_dihedral = e_dihedral+p
       
     
       a = pd
       c = c * a
       s12 = s12 * a
       a11 = c*sb1*s1
       a22 = -sb2 * (2.0*c0*s12 - c*(s1+s2))
       a33 = c*sb3*s2
       a12 = -r12c1 * (c1mag*c*s1 + c2mag*s12)
       a13 = -rb1*rb3*s12
       a23 = r12c2 * (c2mag*c*s2 + c1mag*s12)
       
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

     
       
       !...// apply force to each of 4 atoms
      
       ffh(i1)%x =  ffh(i1)%x + f1(1) 
       ffh(i1)%y =  ffh(i1)%y + f1(2) 
       ffh(i1)%z =  ffh(i1)%z + f1(3) 

       ffh(i2)%x =  ffh(i2)%x + f2(1)
       ffh(i2)%y =  ffh(i2)%y + f2(2)
       ffh(i2)%z =  ffh(i2)%z + f2(3)
       
       ffh(i3)%x =  ffh(i3)%x + f3(1)
       ffh(i3)%y =  ffh(i3)%y + f3(2)
       ffh(i3)%z =  ffh(i3)%z + f3(3)


       ffh(i4)%x =  ffh(i4)%x + f4(1) 
       ffh(i4)%y =  ffh(i4)%y + f4(2) 
       ffh(i4)%z =  ffh(i4)%z + f4(3) 


       
   !    virial(1) = virial(1) + vb1x*f1(1) + vb2x*f3(1) + (vb3x+vb2x)*f4(1)
   !    virial(2) = virial(1) + vb1y*f1(2) + vb2y*f3(2) + (vb3y+vb2y)*f4(2)
   !    virial(3) = virial(1) + vb1z*f1(3) + vb2z*f3(3) + (vb3z+vb2z)*f4(3)
   !    virial(4) = virial(1) + vb1x*f1(2) + vb2x*f3(2) + (vb3x+vb2x)*f4(2)
   !    virial(5) = virial(1) + vb1x*f1(3) + vb2x*f3(3) + (vb3x+vb2x)*f4(3)
   !    virial(6) = virial(1) + vb1y*f1(3) + vb2y*f3(3) + (vb3y+vb2y)*f4(3)







  !     open(unit = 9003, file = 'dihed.dat')
  !     write(9003,*)i1,i2,i3,i4
  !     write(9003,*)p*reps
  !     write(9003,*)f1(1),f1(2),f1(3)
  !     write(9003,*)f2(1),f2(2),f2(3)
  !     write(9003,*)f3(1),f3(2),f3(3)
  !     write(9003,*)f4(1),f4(2),f4(3)
  !     write(9003,*)
     
    enddo
  !  close(9003)

 

  end subroutine dihedral_opls_n2


  subroutine dihedral_opls(step)
    implicit none
    integer :: i1,i2,i3,i4,n,type
    integer :: itype,jtype
    integer :: step
    integer :: offset
    integer :: i
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
    double precision :: lj11,lj12,lj13,lj21,lj22,lj23

    e_dihedral = 0.0
    tolerance = 0.001
    SMALL     = 0.001
    SMALLER   = 0.0001  

 
    
    do n=1,numdihed

       i1 = dihedlist(1,n)
       i2 = dihedlist(2,n)
       i3 = dihedlist(3,n)
       i4 = dihedlist(4,n)
       type = dihedlist(5,n)
      
      ! // 1st bond
       
       vb1x = position(i1)%x - position(i2)%x
       vb1y = position(i1)%y - position(i2)%y
       vb1z = position(i1)%z - position(i2)%z
       call minimg(vb1x,vb1y,vb1z)

       
      ! // 2nd bond

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

      
       
      !// c0 calculation
       
       sb1 = 1.0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z)
       sb2 = 1.0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z)
       sb3 = 1.0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z)
       
       rb1 = sqrt(sb1)
       rb3 = sqrt(sb3)
       
       c0 = (vb1x*vb3x + vb1y*vb3y + vb1z*vb3z) * rb1*rb3
       
       !// 1st and 2nd angle
       
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
       
       !// cos and sin of 2 angles and final c
       
       sin2 = MAX(1.0 - c1mag*c1mag,0.0)
       sc1 = sqrt(sin2)
       if (sc1 < SMALL) sc1 = SMALL
       sc1 = 1.0/sc1
       
       sin2 = MAX(1.0 - c2mag*c2mag,0.0)
       sc2 = sqrt(sin2)
       if (sc2 < SMALL) sc2 = SMALL
       sc2 = 1.0/sc2
       
       s1 = sc1 * sc1
       s2 = sc2 * sc2
       s12 = sc1 * sc2
       c = (c0 + c1mag*c2mag) * s12
       
   !    cx = vb1y*vb2z - vb1z*vb2y
   !    cy = vb1z*vb2x - vb1x*vb2z
   !    cz = vb1x*vb2y - vb1y*vb2x
   !    cmag = sqrt(cx*cx + cy*cy + cz*cz)
   !    dx = (cx*vb3x + cy*vb3y + cz*vb3z)/cmag/b3mag
       
   !    // error check
       if (c.gt.1.0 + TOLERANCE.or.c.lt.(-1.0 - TOLERANCE)) then
          print*,'c and tolerance:',c,TOLERANCE
          print*,'1ST ATOM:',position(i1)%x,position(i1)%y,position(i1)%z
          print*,'2ND ATOM:',position(i2)%x,position(i2)%y,position(i2)%z
          print*,'3RD ATOM:',position(i3)%x,position(i3)%y,position(i3)%z
          print*,'4TH ATOM:',position(i4)%x,position(i4)%y,position(i4)%z
       endif
          
       
       if (c.gt.1.0) c = 1.0
       if (c.lt.-1.0) c = -1.0
       
     
       phi = acos(c)
     !  if (dx.lt.0.0) phi = -1.0d0
       si = sin(phi)
       if (abs(si).lt.SMALLER) si = SMALLER
       siinv = 1.0/si

       k1 = dihedcoeff(1,type); k2 = dihedcoeff(2,type); k3 = dihedcoeff(3,type); k4 = dihedcoeff(4,type)
       
       p = k1*(1.0d0 + c) + k2*(1.0d0 - cos(2.0d0*phi)) + k3*(1.0d0 + cos(3.0d0*phi))+k4*(1.0d0 - cos(4.0d0*phi)) 

       pd = k1 - 2.0d0*k2*sin(2.0d0*phi)*siinv +3.0d0*k3*sin(3.0d0*phi)*siinv-4.0d0*k4*sin(4.0d0*phi)*siinv

       e_dihedral = e_dihedral+p
       
       
       a = pd
       c = c * a
       s12 = s12 * a
       a11 = c*sb1*s1
       a22 = -sb2 * (2.0*c0*s12 - c*(s1+s2))
       a33 = c*sb3*s2
       a12 = -r12c1 * (c1mag*c*s1 + c2mag*s12)
       a13 = -rb1*rb3*s12
       a23 = r12c2 * (c2mag*c*s2 + c1mag*s12)
       
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

       force = 0.0d0
       if(dr2.lt.rcut2)then
          
          dr2i    = 1.0d0/dr2
          dr6i    = dr2i*dr2i*dr2i
          
          offset     = numAtomType*(itype-1)+jtype
          forcelj    = dr6i*(lj14_1(itype,jtype)*dr6i-lj14_2(itype,jtype))
          force      = forcelj*dr2i

      
          pot_14     = pot_14 + dr6i*(dr6i*lj14_3(itype,jtype)-lj14_4(itype,jtype))

          potential = potential + dr6i*(dr6i*lj14_3(itype,jtype)-lj14_4(itype,jtype))
          lj11 = ddx*force
          lj12 = ddy*force
          lj13 = ddz*force
          
          lj21 =  -ddx*force
          lj22 =  -ddy*force
          lj23 =  -ddz*force

       end if
       
       !...// apply force to each of 4 atoms
       
       ffh(i1)%x =  ffh(i1)%x + f1(1) + lj11
       ffh(i1)%y =  ffh(i1)%y + f1(2) + lj12
       ffh(i1)%z =  ffh(i1)%z + f1(3) + lj13


       ffh(i2)%x =  ffh(i2)%x + f2(1)
       ffh(i2)%y =  ffh(i2)%y + f2(2)
       ffh(i2)%z =  ffh(i2)%z + f2(3)
       
       ffh(i3)%x =  ffh(i3)%x + f3(1)
       ffh(i3)%y =  ffh(i3)%y + f3(2)
       ffh(i3)%z =  ffh(i3)%z + f3(3)

       ffh(i4)%x =  ffh(i4)%x + f4(1) + lj21
       ffh(i4)%y =  ffh(i4)%y + f4(2) + lj22
       ffh(i4)%z =  ffh(i4)%z + f4(3) + lj23

       
   !    virial(1) = virial(1) + vb1x*f1(1) + vb2x*f3(1) + (vb3x+vb2x)*f4(1)&
   !         +ddx*ddx*force
          
   !    virial(2) = virial(1) + vb1y*f1(2) + vb2y*f3(2) + (vb3y+vb2y)*f4(2)&
   !         +ddy*ddy*force
   !    virial(3) = virial(1) + vb1z*f1(3) + vb2z*f3(3) + (vb3z+vb2z)*f4(3)&
   !         +ddz*ddz*force
   !    virial(4) = virial(1) + vb1x*f1(2) + vb2x*f3(2) + (vb3x+vb2x)*f4(2)&
   !         +ddx*ddy*force
   !    virial(5) = virial(1) + vb1x*f1(3) + vb2x*f3(3) + (vb3x+vb2x)*f4(3)&
   !         +ddx*ddz*force
   !    virial(6) = virial(1) + vb1y*f1(3) + vb2y*f3(3) + (vb3y+vb2y)*f4(3)&
   !         +ddy*ddz*force


      ! open(unit = 9002, file = 'dihed.dat')
      ! write(9002,*)i1,i2,i3,i4
      ! write(9002,*)f1(1),f1(2),f1(3)
      ! write(9002,*)f2(1),f2(2),f2(3)
      ! write(9002,*)f3(1),f3(2),f3(3)
      ! write(9002,*)f4(1),f4(2),f4(3)
       


    enddo
      
 

  end subroutine dihedral_opls


  subroutine dihedral_opls_PME(step)
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
    double precision :: forcecoul
    double precision :: q1,q2,rrupdate
    double precision :: grij,erfc,expm2


    call system_clock(T1,clock_rate,clock_max)

    e_dihedral = 0.0d0
    tolerance = 0.001
    SMALL     = 0.001
    SMALLER   = 0.0001  

 !   !$omp parallel do schedule(dynamic) reduction(+:potential,e_dihedral,pot_14,e_coul_14) default(firstprivate),&
 !   !$omp& shared(dihedlist,dihedcoeff,position,ff)
    do n=1,numdihed

       i1 = dihedlist(1,n)
       i2 = dihedlist(2,n)
       i3 = dihedlist(3,n)
       i4 = dihedlist(4,n)
       type = dihedlist(5,n)
      
      ! // 1st bond
       
       vb1x = position(i1)%x - position(i2)%x
       vb1y = position(i1)%y - position(i2)%y
       vb1z = position(i1)%z - position(i2)%z
       call minimg(vb1x,vb1y,vb1z)

       
      ! // 2nd bond

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

       
      !// c0 calculation
       
       sb1 = 1.0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z)
       sb2 = 1.0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z)
       sb3 = 1.0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z)
       
       rb1 = sqrt(sb1)
       rb3 = sqrt(sb3)
       
       c0 = (vb1x*vb3x + vb1y*vb3y + vb1z*vb3z) * rb1*rb3
       
       !// 1st and 2nd angle
       
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
       
       !// cos and sin of 2 angles and final c
       
       sin2 = MAX(1.0 - c1mag*c1mag,0.0)
       sc1 = sqrt(sin2)
       if (sc1 < SMALL) sc1 = SMALL
       sc1 = 1.0/sc1
       
       sin2 = MAX(1.0 - c2mag*c2mag,0.0)
       sc2 = sqrt(sin2)
       if (sc2 < SMALL) sc2 = SMALL
       sc2 = 1.0/sc2
       
       s1 = sc1 * sc1
       s2 = sc2 * sc2
       s12 = sc1 * sc2
       c = (c0 + c1mag*c2mag) * s12

   !    // error check
       if (c.gt.1.0 + TOLERANCE.or.c.lt.(-1.0 - TOLERANCE)) then
          print*,'c and tolerance:',c,TOLERANCE
          print*,'1ST ATOM:',position(i1)%x,position(i1)%y,position(i1)%z
          print*,'2ND ATOM:',position(i2)%x,position(i2)%y,position(i2)%z
          print*,'3RD ATOM:',position(i3)%x,position(i3)%y,position(i3)%z
          print*,'4TH ATOM:',position(i4)%x,position(i4)%y,position(i4)%z
       endif
          
       
       if (c.gt.1.0) c = 1.0
       if (c.lt.-1.0) c = -1.0
       
     
       phi = acos(c)

   

       si = sin(phi)
       if (abs(si).lt.SMALLER) si = SMALLER
       siinv = 1.0/si

       k1 = dihedcoeff(1,type); k2 = dihedcoeff(2,type); k3 = dihedcoeff(3,type); k4 = dihedcoeff(4,type)
       
       p = k1*(1.0d0 + c) + k2*(1.0d0 - cos(2.0d0*phi)) + k3*(1.0d0 + cos(3.0d0*phi))+k4*(1.0d0 - cos(4.0d0*phi)) 

       pd = k1 - 2.0d0*k2*sin(2.0d0*phi)*siinv +3.0d0*k3*sin(3.0d0*phi)*siinv-4.0d0*k4*sin(4.0d0*phi)*siinv

       e_dihedral = e_dihedral+p
       

    
       a = pd
       c = c * a
       s12 = s12 * a
       a11 = c*sb1*s1
       a22 = -sb2 * (2.0*c0*s12 - c*(s1+s2))
       a33 = c*sb3*s2
       a12 = -r12c1 * (c1mag*c*s1 + c2mag*s12)
       a13 = -rb1*rb3*s12
       a23 = r12c2 * (c2mag*c*s2 + c1mag*s12)
       
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

       
       !========time to calculate 1-4 PME nonbonded forces
       force = 0.0d0
       if(dr2.lt.cut_coulsq)then
        
          
         ! r = sqrt(dr2)
         ! prefactor =  coulpre*q(i1)*q(i4)/r
         ! grij = gewald*r
         ! expm2 = exp(-grij*grij)
         ! t = 1.0/(1.0 + EWALD_P*grij)
         ! erfc = t*(A1+t*(A2+t*(A3+t*(A4+t*A5))))*expm2
         ! forcecoul = 0.50*prefactor * (erfc + EWALD_F*grij*expm2)
         ! e_coul_14 = e_coul_14 + 0.50*prefactor*erfc   


         ! r = sqrt(dr2)
         ! prefactor =  coulpre*q(i1)*q(i4)/r
         ! grij = gewald*r
         ! expm2 = exp(-grij*grij)
         ! t = 1.0/(1.0 + EWALD_P*grij)
         ! erfc = t*(A1+t*(A2+t*(A3+t*(A4+t*A5))))*expm2
         ! if(dr2.gt.cut_coulsq)prefactor = 0.0d0
         ! forcecoul = 0.50*prefactor * (erfc + EWALD_F*grij*expm2)
         ! e_coul_14 = e_coul_14 + 0.50*prefactor*erfc 


          forcelj = 0.0d0
          if(dr2.lt.rcut2)then
             

             dr2i     =  1.0d0/dr2
             dr6i     =  dr2i*dr2i*dr2i 
             forcelj  =  dr6i*(lj14_1(itype,jtype)*dr6i-lj14_2(itype,jtype))             
             pot_14   =  pot_14 + dr6i*(dr6i*lj14_3(itype,jtype)-lj14_4(itype,jtype))

          end if

!          force = (forcecoul+forcelj)*dr2i
          force = forcelj*dr2i
       endif
   

       !f1(1) = f1(1) + ddx*force
       !f1(2) = f1(2) + ddy*force
       !f1(3) = f1(3) + ddz*force
       
       !f4(1) = f4(1) - ddx*force
       !f4(2) = f4(2) - ddy*force
       !f4(3) = f4(3) - ddz*force

       v_dihedx = v_dihedx  + ddx*ddx*force
       v_dihedy = v_dihedy  + ddy*ddy*force
       v_dihedz = v_dihedz  + ddz*ddz*force
       

       
       !...// apply force to each of 4 atoms
       
     !  !$omp critical
       !ffh(i1)%x =  ffh(i1)%x + f1(1) + ddx*force
       !ffh(i1)%y =  ffh(i1)%y + f1(2) + ddy*force
       !ffh(i1)%z =  ffh(i1)%z + f1(3) + ddz*force

       !ffh(i2)%x =  ffh(i2)%x + f2(1)
       !ffh(i2)%y =  ffh(i2)%y + f2(2)
       !ffh(i2)%z =  ffh(i2)%z + f2(3)
       
       !ffh(i3)%x =  ffh(i3)%x + f3(1)
       !ffh(i3)%y =  ffh(i3)%y + f3(2)
       !ffh(i3)%z =  ffh(i3)%z + f3(3)

       !ffh(i4)%x =  ffh(i4)%x + f4(1) - ddx*force
       !ffh(i4)%y =  ffh(i4)%y + f4(2) - ddy*force
       !ffh(i4)%z =  ffh(i4)%z + f4(3) - ddz*force



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
       

     !  !$omp end critical

    enddo

 !   !$omp end parallel do

    call system_clock(T2,clock_rate,clock_max)
    time_dihed = time_dihed + real(T2-T1)/real(clock_rate)

 !   print*,'time in dihed',real(T2-T1)/real(clock_rate)
      
  end subroutine dihedral_opls_PME


  subroutine dihedral_opls_DSF(step)
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
    double precision :: forcecoul
    double precision :: q1,q2

    call system_clock(T1,clock_rate,clock_max)

    e_dihedral = 0.0d0
    tolerance = 0.001
    SMALL     = 0.001
    SMALLER   = 0.0001  
    
  !  !$omp parallel do reduction(+:potential,e_dihedral,pot_14,e_coul_14) default(firstprivate),&
  !  !$omp& shared(dihedlist,dihedcoeff,position,ff)
    do n=1,numdihed

       i1 = dihedlist(1,n)
       i2 = dihedlist(2,n)
       i3 = dihedlist(3,n)
       i4 = dihedlist(4,n)
       type = dihedlist(5,n)
      
      ! // 1st bond
       
       vb1x = position(i1)%x - position(i2)%x
       vb1y = position(i1)%y - position(i2)%y
       vb1z = position(i1)%z - position(i2)%z
       call minimg(vb1x,vb1y,vb1z)

       
      ! // 2nd bond

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

      
       
      !// c0 calculation
       
       sb1 = 1.0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z)
       sb2 = 1.0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z)
       sb3 = 1.0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z)
       
       rb1 = sqrt(sb1)
       rb3 = sqrt(sb3)
       
       c0 = (vb1x*vb3x + vb1y*vb3y + vb1z*vb3z) * rb1*rb3
       
       !// 1st and 2nd angle
       
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
       
       !// cos and sin of 2 angles and final c
       
       sin2 = MAX(1.0 - c1mag*c1mag,0.0)
       sc1 = sqrt(sin2)
       if (sc1 < SMALL) sc1 = SMALL
       sc1 = 1.0/sc1
       
       sin2 = MAX(1.0 - c2mag*c2mag,0.0)
       sc2 = sqrt(sin2)
       if (sc2 < SMALL) sc2 = SMALL
       sc2 = 1.0/sc2
       
       s1 = sc1 * sc1
       s2 = sc2 * sc2
       s12 = sc1 * sc2
       c = (c0 + c1mag*c2mag) * s12

   !    // error check
       if (c.gt.1.0 + TOLERANCE.or.c.lt.(-1.0 - TOLERANCE)) then
          print*,'c and tolerance:',c,TOLERANCE
          print*,'1ST ATOM:',position(i1)%x,position(i1)%y,position(i1)%z
          print*,'2ND ATOM:',position(i2)%x,position(i2)%y,position(i2)%z
          print*,'3RD ATOM:',position(i3)%x,position(i3)%y,position(i3)%z
          print*,'4TH ATOM:',position(i4)%x,position(i4)%y,position(i4)%z
       endif
          
       
       if (c.gt.1.0) c = 1.0
       if (c.lt.-1.0) c = -1.0
       
     
       phi = acos(c)
       si = sin(phi)
       if (abs(si).lt.SMALLER) si = SMALLER
       siinv = 1.0/si

       k1 = dihedcoeff(1,type); k2 = dihedcoeff(2,type); k3 = dihedcoeff(3,type); k4 = dihedcoeff(4,type)
       
       p = k1*(1.0d0 + c) + k2*(1.0d0 - cos(2.0d0*phi)) + k3*(1.0d0 + cos(3.0d0*phi))+k4*(1.0d0 - cos(4.0d0*phi)) 

       pd = k1 - 2.0d0*k2*sin(2.0d0*phi)*siinv +3.0d0*k3*sin(3.0d0*phi)*siinv-4.0d0*k4*sin(4.0d0*phi)*siinv

       e_dihedral = e_dihedral+p
       
       
       a = pd
       c = c * a
       s12 = s12 * a
       a11 = c*sb1*s1
       a22 = -sb2 * (2.0*c0*s12 - c*(s1+s2))
       a33 = c*sb3*s2
       a12 = -r12c1 * (c1mag*c*s1 + c2mag*s12)
       a13 = -rb1*rb3*s12
       a23 = r12c2 * (c2mag*c*s2 + c1mag*s12)
       
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

       
       !========time to calculate 1-4 nonbonded forces
       force = 0.0d0
       if(dr2.lt.cut_coulsq)then
          r = sqrt(dr2)
          dr2i    = 1.0d0/dr2

          !prefactor =  coulpre*q(i1)*q(i4)/r
          prefactor = coulpre*q(i1)*q(i4)/r
          erfcd = exp(-alpha*alpha*r*r)
          t = 1.0 / (1.0 + EWALD_P*alpha*r)
          erfcc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * erfcd

          !---scale 1-4 interactions by 0.5
          forcecoul = 0.50d0*prefactor * (erfcc/r + 2.0*alpha/MY_PIS * erfcd +&
               r*f_shift) * r
          
          e_coul_14 = e_coul_14 + 0.50d0*prefactor*(erfcc-r*e_shift-dr2*f_shift)


          forcelj = 0.0d0
          if(dr2.lt.rcut2)then
             
             dr6i       = dr2i*dr2i*dr2i
             offset     = numAtomType*(itype-1)+jtype
             forcelj    = dr6i*(lj14_1(itype,jtype)*dr6i-lj14_2(itype,jtype))
             !potential  = potential + dr6i*(dr6i*lj14_3(itype,jtype)-lj14_4(itype,jtype))
             
             pot_14  = pot_14 + dr6i*(dr6i*lj14_3(itype,jtype)-lj14_4(itype,jtype))

          end if


          force = (forcecoul + forcelj)*dr2i       


       endif

   

       f1(1) = f1(1) + ddx*force
       f1(2) = f1(2) + ddy*force
       f1(3) = f1(3) + ddz*force
       
       f4(1) = f4(1) - ddx*force
       f4(2) = f4(2) - ddy*force
       f4(3) = f4(3) - ddz*force


       
       !...// apply force to each of 4 atoms
       
       ffh(i1)%x =  ffh(i1)%x + f1(1)
       ffh(i1)%y =  ffh(i1)%y + f1(2)
       ffh(i1)%z =  ffh(i1)%z + f1(3)

       ffh(i2)%x =  ffh(i2)%x + f2(1)
       ffh(i2)%y =  ffh(i2)%y + f2(2)
       ffh(i2)%z =  ffh(i2)%z + f2(3)
       
       ffh(i3)%x =  ffh(i3)%x + f3(1)
       ffh(i3)%y =  ffh(i3)%y + f3(2)
       ffh(i3)%z =  ffh(i3)%z + f3(3)

       ffh(i4)%x =  ffh(i4)%x + f4(1)
       ffh(i4)%y =  ffh(i4)%y + f4(2)
       ffh(i4)%z =  ffh(i4)%z + f4(3)




   !    virial(1) = virial(1) + vb1x*f1(1) + vb2x*f3(1) + (vb3x+vb2x)*f4(1)&
   !         +ddx*ddx*force
       
   !    virial(2) = virial(1) + vb1y*f1(2) + vb2y*f3(2) + (vb3y+vb2y)*f4(2)&
   !         +ddy*ddy*force
   !    virial(3) = virial(1) + vb1z*f1(3) + vb2z*f3(3) + (vb3z+vb2z)*f4(3)&
   !         +ddz*ddz*force
   !    virial(4) = virial(1) + vb1x*f1(2) + vb2x*f3(2) + (vb3x+vb2x)*f4(2)&
   !         +ddx*ddy*force
   !    virial(5) = virial(1) + vb1x*f1(3) + vb2x*f3(3) + (vb3x+vb2x)*f4(3)&
   !         +ddx*ddz*force
   !    virial(6) = virial(1) + vb1y*f1(3) + vb2y*f3(3) + (vb3y+vb2y)*f4(3)&
   !         +ddy*ddz*force
    enddo
!    !$omp end parallel do

    call system_clock(T2,clock_rate,clock_max)
    time_dihed = time_dihed + real(T2-T1)/real(clock_rate)
      
  end subroutine dihedral_opls_DSF


! -------------------------------------------------------------------------
! CHARMM dihedral torsion, includes LJ & coul 1-4 terms

      subroutine dihedral_charmm(step)
      implicit none

      ! argument variables

      integer :: step,iflag

      !- local variables

      integer :: m,i1,i2,i3,i4,itype,n,ifactor
      integer :: itype1, itype4
      double precision :: rsq14,forcelj,force,r2inv,r6inv,delx,dely,delz
      double precision  :: philj,phicoul,forcecoul
      double precision  :: tolerance,small,vb1x,vb1y,vb1z,vb2x,vb2y
      double precision  :: vb2z,vb2xm,vb2ym,vb2zm,vb3x,vb3y,vb3z,sb1
      double precision  :: sb2,sb3,rb1,rb2,rb3,c0,b1mag2,b1mag,b2mag2
      double precision  :: b2mag,b3mag2,b3mag,ctmp,r12c1,c1mag,r12c2
      double precision  :: c2mag,sc1,sc2,s1,s12,c,p,pd,rc2,a,a11,a22
      double precision  :: a33,a12,a13,a23,sx1,sx2,sx12,sy1,sy2,sy12
      double precision  :: sz1,sz2,sz12,s2
      double precision  :: tmp


      e_dihedral = 0.0
      tolerance = 0.05
      small = 0.001
      iflag = 4

      
      do m = 1,numdihed
        
        i1 = dihedlist(1,m)
        i2 = dihedlist(2,m)
        i3 = dihedlist(3,m)
        i4 = dihedlist(4,m)
        itype = dihedlist(5,m)
        
! 1st bond in dihedral
        
        vb1x = position(i1)%x - position(i2)%x
        vb1y = position(i1)%y - position(i2)%y
        vb1z = position(i1)%z - position(i2)%z
        call minimg(vb1x,vb1y,vb1z)
        
! 2nd bond in dihedral
        
        vb2x = position(i3)%x - position(i2)%x
        vb2y = position(i3)%y - position(i2)%y
        vb2z = position(i3)%z - position(i2)%z
        call minimg(vb2x,vb2y,vb2z)

        vb2xm = -vb2x
        vb2ym = -vb2y
        vb2zm = -vb2z
        call minimg(vb2xm,vb2ym,vb2zm)

! 3rd bond in dihedral
        
        vb3x = position(i4)%x - position(i3)%x
        vb3y = position(i4)%y - position(i3)%y
        vb3z = position(i4)%z - position(i3)%z
        call minimg(vb3x,vb3y,vb3z)

! 1-4 interaction
        
        delx = position(i1)%x - position(i4)%x
        dely = position(i1)%y - position(i4)%y
        delz = position(i1)%z - position(i4)%z
        call minimg(delx,dely,delz)
        rsq14 = delx*delx + dely*dely + delz*delz
        itype1 = type(i1)
        itype4 = type(i4)
        
! c0 calculation
        
        sb1 = 1.0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z)
        sb2 = 1.0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z)
        sb3 = 1.0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z)
        
        rb1 = sqrt(sb1)
        rb2 = sqrt(sb2)
        rb3 = sqrt(sb3)
        
        c0 = (vb1x*vb3x + vb1y*vb3y + vb1z*vb3z) * rb1*rb3
        
! 1st and 2nd angle
        
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
        
! cos and sin of 2 angles and final c
        
        sc1 = sqrt(1.0 - c1mag*c1mag)
        if (sc1.lt.small) sc1 = small
        sc1 = 1.0/sc1
        
        sc2 = sqrt(1.0 - c2mag*c2mag)
        if (sc2.lt.small) sc2 = small
        sc2 = 1.0/sc2
        
        s1 = sc1 * sc1
        s2 = sc2 * sc2
        s12 = sc1 * sc2
        c = (c0 + c1mag*c2mag) * s12
        
! error check

        if (c.gt.1.0+tolerance.or.c.lt.(-1.0-tolerance)) then
          write (6,*) 'Dihedral:',step,&
             i1,i2,i3,i4
          write (6,*) ' Error1:',position(i1)%x,position(i1)%y,position(i1)%z
          write (6,*) ' Error2:',position(i2)%x,position(i2)%y,position(i2)%z
          write (6,*) ' Error3:',position(i3)%x,position(i3)%y,position(i3)%z
          write (6,*) ' Error4:',position(i4)%x,position(i4)%y,position(i4)%z
        endif

        if (c.gt.1.0) c = 1.0
        if (c.lt.-1.0) c = -1.0
        
! energy
!  p = 1 + cos(n*phi) for d = 1
!  p = 1 - cos(n*phi) for d = -1
!  pd = dp/dc / 2
        
        n = nint(dihedcoeff(2,itype))

        if (n.eq.2) then
          p = 2.0*c*c
          pd = 2.0*c
        else if (n.eq.3) then
          rc2 = c*c
          p = (4.0*rc2-3.0)*c + 1.0
          pd = 6.0*rc2 - 1.5
        else if (n.eq.4) then
          rc2 = c*c
          p = 8.0*(rc2-1)*rc2 + 2.0
          pd = (16.0*rc2-8.0)*c
        else if (n.eq.6) then
          rc2 = c*c
          p = ((32.0*rc2-48.0)*rc2 + 18.0)*rc2
          pd = (96.0*(rc2-1.0)*rc2 + 18.0)*c
        else if (n.eq.1) then
          p = c + 1.0
          pd = 0.5
        else if (n.eq.5) then
          rc2 = c*c
          p = ((16.0*rc2-20.0)*rc2 + 5.0)*c + 1.0
          pd = (40.0*rc2-30.0)*rc2 + 2.5
        else if (n.eq.0) then
          p = 2.0
          pd = 0.0
        endif

        if (dihedcoeff(3,itype).eq.180.) then
          p = 2.0 - p
          pd = -pd
        endif


        e_dihedral = e_dihedral + ifactor/4.0 *&
            dihedcoeff(1,itype)*p

        print*,'what the fuck is ifactor'

     !- force on all 4 atoms
        
        a = 2.0 * dihedcoeff(1,itype) * pd
        c = c * a
        s12 = s12 * a
        a11 = (-c*sb1*s1)
        a22 = sb2*(2.0*c0*s12 - c*(s1+s2))
        a33 = (-c*sb3*s2)
        a12 = r12c1*(c1mag*c*s1 + c2mag*s12)
        a13 = rb1*rb3*s12
        a23 = r12c2*(-c2mag*c*s2 - c1mag*s12)
        
        sx1  = a11*vb1x + a12*vb2x + a13*vb3x
        sx2  = a12*vb1x + a22*vb2x + a23*vb3x
        sx12 = a13*vb1x + a23*vb2x + a33*vb3x
        sy1  = a11*vb1y + a12*vb2y + a13*vb3y
        sy2  = a12*vb1y + a22*vb2y + a23*vb3y
        sy12 = a13*vb1y + a23*vb2y + a33*vb3y
        sz1  = a11*vb1z + a12*vb2z + a13*vb3z
        sz2  = a12*vb1z + a22*vb2z + a23*vb3z
        sz12 = a13*vb1z + a23*vb2z + a33*vb3z

        ff(i1)%x = ff(i1)%x - sx1
        ff(i1)%y = ff(i1)%y - sy1
        ff(i1)%z = ff(i1)%z - sz1
        
        ff(i2)%x = ff(i2)%x + (sx2 + sx1)
        ff(i2)%y = ff(i2)%y + (sy2 + sy1)
        ff(i2)%z = ff(i2)%z + (sz2 + sz1)
        
        ff(i3)%x = ff(i3)%x + sx12 - sx2
        ff(i3)%y = ff(i3)%y + sy12 - sy2
        ff(i3)%z = ff(i3)%z + sz12 - sz2
        
        ff(i4)%x = ff(i4)%x - sx12
        ff(i4)%y = ff(i4)%y - sy12
        ff(i4)%z = ff(i4)%z - sz12

        !--- virial contribution: Sum r_ij * F_ij

        if (iflag >= 1) then
          virial(1) = virial(1) - ifactor/4.0 *&
              (vb1x*sx1 + vb2x*sx2 + vb3x*sx12)
          virial(2) = virial(2) - ifactor/4.0 *&
              (vb1y*sy1 + vb2y*sy2 + vb3y*sy12)
         virial(3) = virial(3) - ifactor/4.0 *&
              (vb1z*sz1 + vb2z*sz2 + vb3z*sz12)
         virial(4) = virial(4) - ifactor/4.0 *&
              (vb1x*sy1 + vb2x*sy2 + vb3x*sy12)
          virial(5) = virial(5) - ifactor/4.0 *&
              (vb1x*sz1 + vb2x*sz2 + vb3x*sz12)
         virial(6) = virial(6) - ifactor/4.0 *&
              (vb1y*sz1 + vb2y*sz2 + vb3y*sz12)
        endif

        ! 1-4 LJ & Coulomb interactions

        !    Warning: are assuming here that all nonbond 1-4 interactions (generally 
        !   less that 6 angstroms apart) are within the inner cutoffs,
        !   cutljinterior and cutcoulint (generally at least 8 angstroms)

        if (dihedcoeff(4,itype) /= 0) then

          r2inv = 1.0/rsq14
          r6inv = r2inv*r2inv*r2inv

         ! forcecoul = coulpre*q(i1)*q(i4)*sqrt(r2inv)
          forcelj = r6inv*(lj14_1(itype1,itype4)*r6inv-&
                          lj14_2(itype1,itype4))
          force = dihedcoeff(4,itype)*(forcelj + forcecoul)*r2inv

          ! compute energy if this is thermodynamic timestep

            phicoul = forcecoul
            phicoul = ifactor/4.0 * phicoul * dihedcoeff(4,itype)

            philj = r6inv *&
               (lj14_3(itype1,itype4)*r6inv - lj14_4(itype1,itype4))
            philj = ifactor/4.0 * philj * dihedcoeff(4,itype)

          !  potential = potential + phicoul
          !  potential = potential + philj


            ff(i1)%x = ff(i1)%x + delx*force
            ff(i1)%y = ff(i1)%y + dely*force
            ff(i1)%z = ff(i1)%z + delz*force

            ff(i4)%x = ff(i4)%x - delx*force
            ff(i4)%y = ff(i4)%y - dely*force
            ff(i4)%z = ff(i4)%z - delz*force
          

            force = force*ifactor/4.0

            virial(1) = virial(1) + delx*delx*force
            virial(2) = virial(2) + dely*dely*force
            virial(3) = virial(3) + delz*delz*force
            virial(4) = virial(4) + delx*dely*force
            virial(5) = virial(5) + delx*delz*force
            virial(6) = virial(6) + dely*delz*force

         endif

      enddo

    end subroutine dihedral_charmm

 
   subroutine improper_harmonic(step)
      implicit none
      integer :: step
      integer :: m,i1,i2,i3,i4,itype
      integer :: i
      integer :: T1,T2,clock_rate,clock_max
      double precision :: tolerance,small,v1x,v1y,v1z,v2x,v2y,v2z,v3x
      double precision :: v3y,v3z,ss1,ss2,ss3,r1,r2,r3,c0,c1,c2,s1,s2
      double precision :: s12,c,s,domega,a,a11,a22,a33,a12,a13,a23,sx1
      double precision :: sx2,sx12,sy1,sy2,sy12,sz1,sz2,sz12


      call system_clock(T1,clock_rate,clock_max)

      tolerance = 0.05
      small = 0.001
 !     !$omp parallel do reduction(+:e_improper)  default(firstprivate),&
 !     !$omp& shared(improlist,improcoeff,position,ff)
      do m = 1,numimpro
        
        i1 = improlist(1,m)
        i2 = improlist(2,m)
        i3 = improlist(3,m)
        i4 = improlist(4,m)
        itype = improlist(5,m)
                

        !--- geometry of 4-body
        
        v1x = position(i2)%x - position(i1)%x
        v1y = position(i2)%y - position(i1)%y
        v1z = position(i2)%z - position(i1)%z
        call minimg(v1x,v1y,v1z)
        
        v2x = position(i3)%x - position(i2)%x
        v2y = position(i3)%y - position(i2)%y
        v2z = position(i3)%z - position(i2)%z
        call minimg(v2x,v2y,v2z)
        
        v3x = position(i4)%x - position(i3)%x
        v3y = position(i4)%y - position(i3)%y
        v3z = position(i4)%z - position(i3)%z
        call minimg(v3x,v3y,v3z)
        
        ss1 = 1.0 / (v1x*v1x + v1y*v1y + v1z*v1z)
        ss2 = 1.0 / (v2x*v2x + v2y*v2y + v2z*v2z)
        ss3 = 1.0 / (v3x*v3x + v3y*v3y + v3z*v3z)
        
        r1 = sqrt(ss1)
        r2 = sqrt(ss2)
        r3 = sqrt(ss3)
        
        !--- sin and cos of angle
        
        c0 = -(v1x * v3x + v1y * v3y + v1z * v3z) * r1 * r3
        c1 = -(v1x * v2x + v1y * v2y + v1z * v2z) * r1 * r2
        c2 = -(v3x * v2x + v3y * v2y + v3z * v2z) * r3 * r2
        
        s1 = 1.0 - c1*c1
        if (s1.lt.small) s1 = small
        s1 = 1.0 / s1
        
        s2 = 1.0 - c2*c2
        if (s2.lt.small) s2 = small
        s2 = 1.0 / s2
        
        s12 = sqrt(s1*s2)
        c = (c1*c2 + c0) * s12
        
        !--- error check
        
        if (c.gt.1.0+tolerance.or.c.lt.(-1.0-tolerance)) then
          write (6,*) 'Improper:',step,&
             (i1),(i2),(i3),(i4)
          write (6,*) ' Error1:',position(i1)%x,position(i1)%y,position(i1)%z
          write (6,*) ' Error2:',position(i2)%x,position(i2)%y,position(i2)%z
          write (6,*) ' Error3:',position(i3)%x,position(i3)%y,position(i3)%z
          write (6,*) ' Error4:',position(i4)%x,position(i4)%y,position(i4)%z
        endif

        if (c.gt.1.0) c = 1.0
        if (c.lt.-1.0) c = -1.0
        
        s = sqrt(1.0 - c*c)
        if (s.lt.small) s = small
        
        !--- energy
        domega = acos(c) - improcoeff(2,itype)
        a = improcoeff(1,itype) * domega

        e_improper = e_improper +  a*domega

        
        !--- force on all 4 atoms
        
        a = -a * 2.0/s
        c = c * a
        
        s12 = s12 * a
        a11 = (-c * ss1 * s1)
        a22 = ss2 * (2.0 * c0 * s12 - c * (s1 + s2))
        a33 = (-c * ss3 * s2)
        a12 = r1 * r2 * (c1 * c * s1 + c2 * s12)
        a13 = r1 * r3 * s12
        a23 = r2 * r3 * (-c2 * c * s2 - c1 * s12)
        
        sx1  = a12*v2x + a13*v3x - a11*v1x
        sx2  = a22*v2x + a23*v3x - a12*v1x
        sx12 = a23*v2x + a33*v3x - a13*v1x
        sy1  = a12*v2y + a13*v3y - a11*v1y
        sy2  = a22*v2y + a23*v3y - a12*v1y
        sy12 = a23*v2y + a33*v3y - a13*v1y
        sz1  = a12*v2z + a13*v3z - a11*v1z
        sz2  = a22*v2z + a23*v3z - a12*v1z
        sz12 = a23*v2z + a33*v3z - a13*v1z
        
        ffh(i1)%x = ffh(i1)%x - sx1
        ffh(i1)%y = ffh(i1)%y - sy1
        ffh(i1)%z = ffh(i1)%z - sz1
        
        ffh(i2)%x = ffh(i2)%x + (sx2 + sx1)
        ffh(i2)%y = ffh(i2)%y + (sy2 + sy1)
        ffh(i2)%z = ffh(i2)%z + (sz2 + sz1)
        
        ffh(i3)%x = ffh(i3)%x + sx12 - sx2
        ffh(i3)%y = ffh(i3)%y + sy12 - sy2
        ffh(i3)%z = ffh(i3)%z + sz12 - sz2
        
        ffh(i4)%x = ffh(i4)%x - sx12
        ffh(i4)%y = ffh(i4)%y - sy12
        ffh(i4)%z = ffh(i4)%z - sz12



    !    open(unit = 9005, file = 'impro.dat')
    !    write(9005,*)i1,i2,i3,i4
    !    write(9005,*)a*domega*reps
    !    write(9005,*)i1,-sx1,-sy1,-sz1
    !    write(9005,*)i2,sx2+sx1,sy2+sy1,sz2+sz1
    !    write(9005,*)i3,sx12-sx2,sy12-sy2,sz12-sz2
    !    write(9005,*)i4,-sx12,-sy12,-sz12
    !    write(9005,*)

        
        !--- virial contribution: Sum r_ij * F_ij
   !     virial(1) = virial(1) + (v1x*sx1 - v2x*sx2 - v3x*sx12)
   !     virial(2) = virial(2) + (v1y*sy1 - v2y*sy2 - v3y*sy12)
   !     virial(3) = virial(3) + (v1z*sz1 - v2z*sz2 - v3z*sz12)
   !     virial(4) = virial(4) + (v1x*sy1 - v2x*sy2 - v3x*sy12)
   !     virial(5) = virial(5) + (v1x*sz1 - v2x*sz2 - v3x*sz12)
   !     virial(6) = virial(6) + (v1y*sz1 - v2y*sz2 - v3y*sz12)
      enddo
   !   !$omp end parallel do

      call system_clock(T2,clock_rate,clock_max)
   !   print*,'time in impropoer',real(T2-T1)/real(clock_rate)
      
      time_impro = time_impro + real(T2-T1)/real(clock_rate)

    end subroutine improper_harmonic



     subroutine improper_harmonic_lammps14(step)
      implicit none
      integer :: step
      integer :: m,i1,i2,i3,i4,itype
      integer :: i
      integer :: T1,T2,clock_rate,clock_max
      double precision :: tolerance,small,v1x,v1y,v1z,v2x,v2y,v2z,v3x
      double precision :: v3y,v3z,ss1,ss2,ss3,r1,r2,r3,c0,c1,c2,s1,s2
      double precision :: s12,c,s,domega,a,a11,a22,a33,a12,a13,a23,sx1
      double precision :: sx2,sx12,sy1,sy2,sy12,sz1,sz2,sz12
      double precision :: f1(3),f2(3),f3(3),f4(3)
      double precision :: vb1x,vb3x,vb3y,vb1z,vb3z,vb2x,vb2y,vb2z,vb1y

      call system_clock(T1,clock_rate,clock_max)


      tolerance = 0.05
      small = 0.001
      do m = 1,numimpro
        
        i1 = improlist(1,m)
        i2 = improlist(2,m)
        i3 = improlist(3,m)
        i4 = improlist(4,m)
        itype = improlist(5,m)
                
 

        v1x = position(i1)%x - position(i2)%x
        v1y = position(i1)%y - position(i2)%y
        v1z = position(i1)%z - position(i2)%z
        call minimg(v1x,v1y,v1z)
        
        v2x = position(i3)%x - position(i2)%x
        v2y = position(i3)%y - position(i2)%y
        v2z = position(i3)%z - position(i2)%z
        call minimg(v2x,v2y,v2z)
        
        v3x = position(i4)%x - position(i3)%x
        v3y = position(i4)%y - position(i3)%y
        v3z = position(i4)%z - position(i3)%z
        call minimg(v3x,v3y,v3z)


        
        ss1 = 1.0 / (v1x*v1x + v1y*v1y + v1z*v1z)
        ss2 = 1.0 / (v2x*v2x + v2y*v2y + v2z*v2z)
        ss3 = 1.0 / (v3x*v3x + v3y*v3y + v3z*v3z)

 
        r1 = sqrt(ss1)
        r2 = sqrt(ss2)
        r3 = sqrt(ss3)
        
        
        ! sin and cos of angle
        
        c0 = (v1x * v3x + v1y * v3y + v1z * v3z) * r1 * r3
        c1 = (v1x * v2x + v1y * v2y + v1z * v2z) * r1 * r2
        c2 = -(v3x * v2x + v3y * v2y + v3z * v2z) * r3 * r2
    
        
        s1 = 1.0 - c1*c1
        if (s1 < SMALL) s1 = SMALL
        s1 = 1.0 / s1
        
        s2 = 1.0 - c2*c2
        if (s2 < SMALL) s2 = SMALL
        s2 = 1.0 / s2

        s12 = sqrt(s1*s2)
        c = (c1*c2 + c0) * s12
        
        ! error check
        
        if (c .gt. 1.0 + TOLERANCE.or.c .lt. (-1.0 - TOLERANCE)) then
           print*,'error in improper harmonic'
        endif
        
        if (c > 1.0) c = 1.0
        if (c < -1.0) c = -1.0
        
        s = sqrt(1.0 - c*c)
        if (s < SMALL) s = SMALL
        
        ! force & energy
        
        domega = acos(c) - improcoeff(2,itype)
        a = improcoeff(1,itype)* domega
        
        e_improper = e_improper +  a*domega
        

        a = -a * 2.0/s
        c = c * a
        s12 = s12 * a
        a11 = c*ss1*s1
        a22 = -ss2 * (2.0*c0*s12 - c*(s1+s2))
        a33 = c*ss3*s2
        a12 = -r1*r2*(c1*c*s1 + c2*s12)
        a13 = -r1*r3*s12
        a23 = r2*r3*(c2*c*s2 + c1*s12)
        
        sx2  = a22*v2x + a23*v3x + a12*v1x
        sy2  = a22*v2y + a23*v3y + a12*v1y
        sz2  = a22*v2z + a23*v3z + a12*v1z
        
        f1(1) = a12*v2x + a13*v3x + a11*v1x
        f1(2) = a12*v2y + a13*v3y + a11*v1y
        f1(3) = a12*v2z + a13*v3z + a11*v1z
        
        f2(1) = -sx2 - f1(1)
        f2(2) = -sy2 - f1(2)
        f2(3) = -sz2 - f1(3)
        
        f4(1) = a23*v2x + a33*v3x + a13*v1x
        f4(2) = a23*v2y + a33*v3y + a13*v1y
        f4(3) = a23*v2z + a33*v3z + a13*v1z
        
        f3(1) = sx2 - f4(1)
        f3(2) = sy2 - f4(2)
        f3(3) = sz2 - f4(3)
        
        ! apply force to each of 4 atoms
        
        !ffh(i1)%x = f1(1) + ffh(i1)%x
        !ffh(i1)%y = f1(2) + ffh(i1)%y
        !ffh(i1)%z = f1(3) + ffh(i1)%z
        
        
        !ffh(i2)%x = f2(1) + ffh(i2)%x
        !ffh(i2)%y = f2(2) + ffh(i2)%y
        !ffh(i2)%z = f2(3) + ffh(i2)%z
        
        
        !ffh(i3)%x = f3(1) + ffh(i3)%x
        !ffh(i3)%y = f3(2) + ffh(i3)%y
        !ffh(i3)%z = f3(3) + ffh(i3)%z
        
        
        !ffh(i4)%x = f4(1) + ffh(i4)%x
        !ffh(i4)%y = f4(2) + ffh(i4)%y
        !ffh(i4)%z = f4(3) + ffh(i4)%z



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




        v_improx = v_improx + v1x*f1(1) + v2x*f3(1) +&
             (v3x+v2x)*f4(1)
        v_improy = v_improy+ v1y*f1(2) + v2y*f3(2) +&
             (v3y+v2y)*f4(2)
        v_improz = v_improz +  v1z*f1(3) + v2z*f3(3) +&
             (v3z+v2z)*f4(3)
        
        
        
     enddo
   end subroutine improper_harmonic_lammps14




  end module mod_force_dihed
