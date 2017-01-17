module mod_prune
  use global
  use mod_minimg
  implicit none
  contains

    subroutine prune
      implicit none
      integer :: i,j,atom1,atom2,atom3,atom4
      double precision :: vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb2xm,vb2ym,vb2zm
      double precision :: sb1,sb2,sb3,rb1,rb3,c0,b1mag2,b1mag,b2mag2
      double precision :: b2mag,b3mag2,b3mag,ctmp,r12c1,c1mag,r12c2
      double precision :: c2mag,sc1,sc2,s1,s12,c,p,pd,a,a11,a22
      double precision :: s2,cx,cy,cz,cmag,dx,phi,si,siinv,sin2
      double precision :: small, smaller
      double precision :: dihedangle(2)
      double precision :: angle1,angle2
      double precision :: signt(2)
      double precision :: abs1,abs2
      double precision :: u(3),v(3),w(3),n(3),umag
      double precision :: aphi,dot

      integer :: spot1,spot2,spot3,spot4,spot5
      integer :: spot11,spot12,spot13,spot14,spot15

      

      SMALL     = 0.001
      SMALLER   = 0.0001

     !open(unit = 5032, file = 'check.dat')

      do i = 1,numMolArray(1)

         spot1 = (i-1)*numatns +1
         spot2 = (i-1)*numatns +2
         spot3 = (i-1)*numatns +3
         spot4 = (i-1)*numatns +4
         spot5 = (i-1)*numatns +5

         spot11 = (i-1)*numatompertemp +1
         spot12 = (i-1)*numatompertemp +2
         spot13 = (i-1)*numatompertemp +3
         spot14 = (i-1)*numatompertemp +4
         
         pruned(spot11)%x = ptta(spot1)%x
         pruned(spot11)%y = ptta(spot1)%y
         pruned(spot11)%z = ptta(spot1)%z


         pruned(spot12)%x = ptta(spot2)%x
         pruned(spot12)%y = ptta(spot2)%y
         pruned(spot12)%z = ptta(spot2)%z
         
         pruned(spot13)%x = ptta(spot3)%x
         pruned(spot13)%y = ptta(spot3)%y
         pruned(spot13)%z = ptta(spot3)%z


      !   pruned(spot14)%x = ptta(spot4)%x
      !   pruned(spot14)%y = ptta(spot4)%y
      !   pruned(spot14)%z = ptta(spot4)%z


         !---now obtain dihedrals
         do j = 1,2
            atom1 = (i-1)*numatns + 1
            atom2 = (i-1)*numatns + 2
            atom3 = (i-1)*numatns + 3
            atom4 = (i-1)*numatns + 4 + j-1

 
            ! // 1st bond
            vb1x = ptta(atom3)%x  -ptta(atom2)%x
            vb1y = ptta(atom3)%y - ptta(atom2)%y
            vb1z = ptta(atom3)%z - ptta(atom2)%z
            call minimg(vb1x,vb1y,vb1z)
            

            
            ! // 2nd bond
            vb2x =  ptta(atom1)%x - ptta(atom2)%x
            vb2y =  ptta(atom1)%y - ptta(atom2)%y
            vb2z =  ptta(atom1)%z - ptta(atom2)%z
            call minimg(vb2x,vb2y,vb2z)

            
            vb2xm = -vb2x
            vb2ym = -vb2y
            vb2zm = -vb2z
            call minimg(vb2xm,vb2ym,vb2zm)
            
           
            ! // 3rd bond
            vb3x = ptta(atom4)%x - ptta(atom1)%x
            vb3y = ptta(atom4)%y - ptta(atom1)%y
            vb3z = ptta(atom4)%z - ptta(atom1)%z
            call minimg(vb3x,vb3y,vb3z)
            
          
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
            
            
            
            if (c.gt.1.0) c = 1.0
            if (c.lt.-1.0) c = -1.0
            
            !---returns angles in radians
            phi = acos(c)

            !---assumes carbonyl o is (1,i)
            !---assumes OH o is (2,i)
            !dihedprune(j,i) = phi

            dihedangle(j) = phi


            if(phi.gt.0.0d0)then
               signt(j) = 1.0d0
            elseif(signt(j).lt.0.0d0)then
               signt(j) = -1.0d0
            endif

         enddo

         abs1 = abs(dihedangle(1))
         abs2 = abs(dihedangle(2))
         
         if(abs1.lt.abs2)then
            angle1 = abs1*signt(1)
            angle2 = abs2*signt(2) - signt(2)*acos(-1.0d0)
         else
            angle1 = abs2*signt(2)
            angle2 = abs1*signt(1) - signt(1)*acos(-1.0d0)
         endif


         !---now build u vector
         u(:) = (/vb2x,vb2y,vb2z/)
         umag = sqrt(vb2x*vb2x + vb2y*vb2y + vb2z*vb2z)
         u(:) = u(:)/umag
         
         !---now build n
         n(:) = (/vb1x,vb1y,vb1z/)
!         n(:) = (/vb3x,vb3y,vb3z/)

         dot  = n(1)*u(1) + n(2)*u(2) + n(3)*u(3)
         n(:) = n(:)-dot*u(:)
         !---get w from n
         umag = sqrt(n(1)*n(1) + n(2)*n(2) + n(3)*n(3))
         w(:) = n(:)/umag
         
         !---now build v vector
         v(1) = w(2)*u(3)-w(3)*u(2)
         v(2) = w(3)*u(1)-w(1)*u(3)
         v(3) = w(1)*u(2) -w(2)*u(1)
         
         aphi = 0.50d0*(angle1+angle2)
         
         pruned(spot14)%x = CObond*w(1)*(1.0d0+cos(2.0d0*aphi)) + CObond*v(1)*sin(2.0d0*aphi)& 
              + pruned(spot11)%x
         
         pruned(spot14)%y = CObond*w(2)*(1.0d0+cos(2.0d0*aphi)) + CObond*v(2)*sin(2.0d0*aphi)& 
              + pruned(spot11)%y
         
         pruned(spot14)%z = CObond*w(3)*(1.0d0+cos(2.0d0*aphi)) + CObond*v(3)*sin(2.0d0*aphi)& 
              + pruned(spot11)%z
     
    !     write(5032,*)pruned(spot11)%x,pruned(spot11)%y,pruned(spot11)%z
    !     write(5032,*)pruned(spot12)%x,pruned(spot12)%y,pruned(spot12)%z
    !     write(5032,*)pruned(spot13)%x,pruned(spot13)%y,pruned(spot13)%z
    !     write(5032,*)pruned(spot14)%x,pruned(spot14)%y,pruned(spot14)%z
    !     write(5032,*)dihedangle(1),dihedangle(2)

    !     write(5032,*)'---------------------------------------------------'
      enddo
      
    end subroutine prune

    subroutine prune_template
      implicit none
      integer :: i,j,atom1,atom2,atom3,atom4,k
      double precision :: vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb2xm,vb2ym,vb2zm
      double precision :: sb1,sb2,sb3,rb1,rb3,c0,b1mag2,b1mag,b2mag2
      double precision :: b2mag,b3mag2,b3mag,ctmp,r12c1,c1mag,r12c2
      double precision :: c2mag,sc1,sc2,s1,s12,c,p,pd,a,a11,a22
      double precision :: s2,cx,cy,cz,cmag,dx,phi,si,siinv,sin2
      double precision :: small, smaller
      double precision :: dihedangle(2)
      double precision :: angle1,angle2
      double precision :: signt(2)
      double precision :: abs1,abs2
      double precision :: u(3),v(3),w(3),n(3),umag
      double precision :: aphi,dot

      integer :: spot1,spot2,spot3,spot4,spot5
      integer :: spot11,spot12,spot13,spot14,spot15

      

      SMALL     = 0.001
      SMALLER   = 0.0001  

      do j = 1,numtemplate
         do i = 1,nummoltemp
            
            spot1 = (i-1)*numatns +1
            spot2 = (i-1)*numatns +2
            spot3 = (i-1)*numatns +3
            
            spot11 = (i-1)*numatompertemp +1
            spot12 = (i-1)*numatompertemp +2
            spot13 = (i-1)*numatompertemp +3
            spot14 = (i-1)*numatompertemp +4
            
            tempa(j)%temp(1,spot11) = tempa_noprune(j)%temp(1,spot1)
            tempa(j)%temp(2,spot11) = tempa_noprune(j)%temp(2,spot1)
            tempa(j)%temp(3,spot11) = tempa_noprune(j)%temp(3,spot1)
            
            
            tempa(j)%temp(1,spot12) = tempa_noprune(j)%temp(1,spot2)
            tempa(j)%temp(2,spot12) = tempa_noprune(j)%temp(2,spot2)
            tempa(j)%temp(3,spot12) = tempa_noprune(j)%temp(3,spot2)
            
            tempa(j)%temp(1,spot13) = tempa_noprune(j)%temp(1,spot3)
            tempa(j)%temp(2,spot13) = tempa_noprune(j)%temp(2,spot3)
            tempa(j)%temp(3,spot13) = tempa_noprune(j)%temp(3,spot3)

            
            
            !---now obtain dihedrals
            do k = 1,2
               atom1 = (i-1)*numatns + 3
               atom2 = (i-1)*numatns + 2
               atom3 = (i-1)*numatns + 1
               atom4 = (i-1)*numatns + 4 + k-1
               
               ! // 1st bond
               vb1x = tempa_noprune(j)%temp(1,atom1) - tempa_noprune(j)%temp(1,atom2)
               vb1y = tempa_noprune(j)%temp(2,atom1) - tempa_noprune(j)%temp(2,atom2)
               vb1z = tempa_noprune(j)%temp(3,atom1) - tempa_noprune(j)%temp(3,atom2)
               call minimg(vb1x,vb1y,vb1z)
               
               
               ! // 2nd bond
               vb2x =  tempa_noprune(j)%temp(1,atom3) - tempa_noprune(j)%temp(1,atom2)
               vb2y =  tempa_noprune(j)%temp(2,atom3) - tempa_noprune(j)%temp(2,atom2)
               vb2z =  tempa_noprune(j)%temp(3,atom3) - tempa_noprune(j)%temp(3,atom2)
               call minimg(vb2x,vb2y,vb2z)
               
               
               vb2xm = -vb2x
               vb2ym = -vb2y
               vb2zm = -vb2z
               call minimg(vb2xm,vb2ym,vb2zm)
               
               
               ! // 3rd bond
               vb3x = tempa_noprune(j)%temp(1,atom4) - tempa_noprune(j)%temp(1,atom3)
               vb3y = tempa_noprune(j)%temp(2,atom4) - tempa_noprune(j)%temp(2,atom3)
               vb3z = tempa_noprune(j)%temp(3,atom4) - tempa_noprune(j)%temp(3,atom3)
               call minimg(vb3x,vb3y,vb3z)
               
               
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
               
               
               
               if (c.gt.1.0) c = 1.0
               if (c.lt.-1.0) c = -1.0
               
               !---returns angles in radians
               phi = acos(c)
               
               !---assumes carbonyl o is (1,i)
               !---assumes OH o is (2,i)
               !dihedprune(j,i) = phi
               
               dihedangle(k) = phi

               
               
               if(phi.gt.0.0d0)then
                  signt(k) = 1.0d0
               elseif(signt(k).lt.0.0d0)then
                  signt(k) = -1.0d0
               endif
               
            enddo
            
            abs1 = abs(dihedangle(1))
            abs2 = abs(dihedangle(2))
            
            if(abs1.lt.abs2)then
               angle1 = abs1*signt(1)
               angle2 = abs2*signt(2) - signt(2)*acos(-1.0d0)
            else
               angle1 = abs2*signt(2)
               angle2 = abs1*signt(1) - signt(1)*acos(-1.0d0)
            endif
            
            
            !---now build u vector
            u(:) = (/vb2x,vb2y,vb2z/)
            umag = sqrt(vb2x*vb2x + vb2y*vb2y + vb2z*vb2z)
            u(:) = u(:)/umag
            
            !---now build n
            n(:) = (/vb1x,vb1y,vb1z/)
            !n(:) = (/vb3x,vb3y,vb3z/)
            dot  = n(1)*u(1) + n(2)*u(2) + n(3)*u(3)
            n(:) = n(:)-dot*u(:)

            !---get w from n
            umag = sqrt(n(1)*n(1) + n(2)*n(2) + n(3)*n(3))
            w(:) = n(:)/umag
            
            !---now build v vector
            v(1) = w(2)*u(3) -w(3)*u(2)
            v(2) = w(3)*u(1) -w(1)*u(3)
            v(3) = w(1)*u(2) -w(2)*u(1)
            
            aphi = 0.50d0*(angle1+angle2)
            
            tempa(j)%temp(1,spot14) = CObond*w(1)*(1.0d0+cos(2.0d0*aphi)) + CObond*v(1)*sin(2.0d0*aphi)& 
                 + tempa(j)%temp(1,spot11)
           
            tempa(j)%temp(2,spot14) = CObond*w(2)*(1.0d0+cos(2.0d0*aphi)) + CObond*v(2)*sin(2.0d0*aphi)& 
                 +  tempa(j)%temp(2,spot11)
           
            tempa(j)%temp(3,spot14) = CObond*w(3)*(1.0d0+cos(2.0d0*aphi)) + CObond*v(3)*sin(2.0d0*aphi)& 
                 +  tempa(j)%temp(3,spot11)
            
           ! tempa(j)%temp(1,spot14) =  tempa(j)%temp(1,spot11)
           ! tempa(j)%temp(2,spot14) =  tempa(j)%temp(2,spot11)
           ! tempa(j)%temp(3,spot14) =  tempa(j)%temp(3,spot11)


         enddo
      enddo
    end subroutine prune_template


    subroutine prune_template_debug
         implicit none
         integer :: i,j,atom1,atom2,atom3,atom4,k
         double precision :: vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb2xm,vb2ym,vb2zm
         double precision :: sb1,sb2,sb3,rb1,rb3,c0,b1mag2,b1mag,b2mag2
         double precision :: b2mag,b3mag2,b3mag,ctmp,r12c1,c1mag,r12c2
         double precision :: c2mag,sc1,sc2,s1,s12,c,p,pd,a,a11,a22
         double precision :: s2,cx,cy,cz,cmag,dx,phi,si,siinv,sin2
         double precision :: small, smaller
         double precision :: dihedangle(2)
         double precision :: angle1,angle2
         double precision :: signt(2)
         double precision :: abs1,abs2
         double precision :: u(3),v(3),w(3),n(3),umag
         double precision :: aphi,dot
         
         integer :: spot1,spot2,spot3,spot4,spot5,count
         integer :: spot11,spot12,spot13,spot14,spot15

         tempa(1)%temp(:,:) = tempa_noprune(1)%temp(:,:)
         tempa(2)%temp(:,:) = tempa_noprune(2)%temp(:,:)
         tempa(3)%temp(:,:) = tempa_noprune(3)%temp(:,:)
         tempa(4)%temp(:,:) = tempa_noprune(4)%temp(:,:)
         tempa(5)%temp(:,:) = tempa_noprune(5)%temp(:,:)
         tempa(6)%temp(:,:) = tempa_noprune(6)%temp(:,:)

         
         do i =1,numTemPlate

            count = 0
            do j = 1,nummoltemp
               do k = 1,numatns

                     spot5  = (j-1)*10+k
                     count = count + 1
                  !   tempa(i)%temp(1,count) = tempa_noprune(i)%temp(1,spot5)
                  !   tempa(i)%temp(2,count) = tempa_noprune(i)%temp(2,spot5)
                  !   tempa(i)%temp(3,count) = tempa_noprune(i)%temp(3,spot5)

               enddo
            enddo
         enddo
       end subroutine prune_template_debug


    subroutine prune_debug
      implicit none
      integer :: i,j,spot

      do i = 1,numMolArray(1)
         do j = 1,numatns

            spot = (i-1)*numatns+j
            pruned(spot)%x = ptta(spot)%x
            pruned(spot)%y = ptta(spot)%y
            pruned(spot)%z = ptta(spot)%z
         enddo
      enddo

    end subroutine prune_debug

  end module mod_prune
