module mod_tip4p
  use global
  use mod_minimg
  contains

    subroutine calculate_m_posit
      implicit none
      double precision :: o_x1, o_y1, o_z1
      double precision :: h1_x1, h1_y1, h1_z1
      double precision :: h2_x1, h2_y1, h2_z1
      double precision :: a1,a2,a3
      double precision :: c1,c2,c3
      integer :: i
      integer :: atom_o, atom_h1,atom_h2,atom_m

      do i = 1,num_tip4p_pairs
         atom_o  = tip4p_list(1,i)
         atom_h1 = tip4p_list(2,i)
         atom_h2 = tip4p_list(3,i)
         atom_m  = tip4p_list(4,i)

         o_x1  = position(atom_o)%x;  o_y1  = position(atom_o)%y;  o_z1  = position(atom_o)%z
         h1_x1 = position(atom_h1)%x; h1_y1 = position(atom_h1)%y; h1_z1 = position(atom_h1)%z
         h2_x1 = position(atom_h2)%x; h2_y1 = position(atom_h2)%y; h2_z1 = position(atom_h2)%z

         a1 = h1_x1 - o_x1
         a1 = a1 - box*nint(a1*ibox)
         
         a2 = h1_y1 - o_y1
         a2 = a2 - box*nint(a2*ibox)
         
         a3 = h1_z1 - o_z1
         a3 = a3 - box*nint(a3*ibox)
         
         c1 = h2_x1 - o_x1
         c1 = c1 - box*nint(c1*ibox)
         
         c2 = h2_y1  - o_y1
         c2 = c2 - box*nint(c2*ibox)
         
         c3 = h2_z1 -  o_z1
         c3 = c3 - box*nint(c3*ibox)
         
         maga = a1*a1 + a2*a2 + a3*a3
         magc = c1*c1 + c2*c2 + c3*c3
         maga = sqrt(maga)
         magc = sqrt(magc)
         
         gam1 = maga*c1 + magc*a1; gam2 = maga*c2 + magc*a2; gam3 = maga*c3 + magc*a3
         magy = gam1*gam1 + gam2*gam2 + gam3*gam3
         magy = sqrt(magy)
         
         gam1 = gam1/magy; gam2 = gam2/magy; gam3 = gam3/magy
         x4   = gam1*r_om; y4 = gam2*r_om; z4 = gam3*r_om
     
         
         x4   = x4 + x1
         y4   = y4 + y1
         z4   = z4 + z1

         position(atom_m)%x = x4
         position(atom_m)%y = y4
         position(atom_m)%z = z4

      enddo


    subroutine redistrib_m_water_force
      implicit none
      double precision :: fad_ox(3), fad_h(3)
      double precision :: r_ox_m(3)
      double precision :: invlen_r_ox_m,len_r_ox_m
      double precision :: rad_factor(3)
      double precision :: f_rad(3)
      double precision :: dx,dy,dz,dr2
      double precision :: oxcomp,hydcomp
      integer :: i
      integer :: atom_o,atom_h1,atom_h2,atom_m

      
      do i = 1,num_tip4p_pairs
         fad_ox(:) = 0.0d0; fad_h(:) =0.0d0


         atom_o  = tip4p_list(1,i)
         atom_h1 = tip4p_list(2,i)
         atom_h2 = tip4p_list(3,i)
         atom_m  = tip4p_list(4,i)

         !---calculate radial component and add to oxygen
         dx = nopbc(atom_m)%x - nopbc(atom_o)%x
         dy = nopbc(atom_m)%y - nopbc(atom_o)%y
         dz = nopbc(atom_m)%z - nopbc(atom_o)%z
         
         
         r_ox_m(1) = dx; r_ox_m(2) = dy; r_ox_m(3) = dz

         len_r_ox_m    = sqrt(dx*dx + dy*dy + dz*dz)
         invlen_r_ox_m = 1.0d0/len_r_ox_m

         rad_factor(1) = ff(atom_4)%x*dx*invlen_r_ox_m
         rad_factor(2) = ff(atom_4)%y*dy*invlen_r_ox_m         
         rad_factor(3) = ff(atom_4)%z*dz*invlen_r_ox_m
         
         f_rad(1) = r_ox_m(1)*rad_factor(1)
         f_rad(2) = r_ox_m(2)*rad_factor(2)
         f_rad(3) = r_ox_m(3)*rad_factor(3)

         fad_ox(:) = fad_ox(:) + f_rad(:)

         !---angular component
         r_hcom_ox(1) = nopbc(atom_o)%x -&
              (nopbc(atom_h1)%x+nopbc(atom_h2)%x)*0.5

         r_hcom_ox(2) = nopbc(atom_o)%y -&
              (nopbc(atom_h1)%y+nopbc(atom_h2)%y)*0.5
         
         r_hcom_ox(3) = nopbc(atom_o)%z -&
              (nopbc(atom_h1)%z+nopbc(atom_h2)%z)*0.5

         len_r_hcom_ox = sqrt(r_hcom_ox(1)*r_hcom_ox(1) + r_hcom_ox(2)*r_hcom_ox(2) + r_hcom_ox(3)*r_hcom_ox(3))
         invlen_r_hcom_ox  = 1.0d0/len_r_hcom_ox


         r_h2_h1_2(1) = (nopbc(atom_h1)%x-nopbc(atom_h2)%x)*0.50
         r_h2_h1_2(2) = (nopbc(atom_h1)%y-nopbc(atom_h2)%y)*0.50
         r_h2_h1_2(3) = (nopbc(atom_h1)%z-nopbc(atom_h2)%z)*0.50

         
         f_ang(1) = ff(atom_m)%x - f_rad(1)
         f_ang(2) = ff(atom_m)%y - f_rad(2)
         f_ang(3) = ff(atom_m)%z - f_rad(3)

         
         oxcomp = (len_r_hcom_ox-len_r_ox_m)*invlen_r_hcom_ox
         hydcomp = 0.5*len_r_ox_m*invlen_r_hcom_ox
         
         fad_ox(1) = fad_ox(1) + f_ang(1)*oxcomp
         fad_ox(2) = fad_ox(2) + f_ang(2)*oxcomp
         fad_ox(3) = fad_ox(3) + f_ang(3)*oxcomp

         fad_h(1) = f_ang(1) * hydcomp
         fad_h(2) = f_ang(2) * hydcomp
         fad_h(3) = f_ang(3) * hydcomp

         ff(atom_o)%x = ff(atom_o)%x + fad_ox(1)
         ff(atom_o)%y = ff(atom_o)%y + fad_ox(2)
         ff(atom_o)%z = ff(atom_o)%z + fad_ox(3)

         ff(atom_h1)%x = ff(atom_h1)%x + fad_h(1)
         ff(atom_h1)%y = ff(atom_h1)%y + fad_h(2)
         ff(atom_h1)%z = ff(atom_h1)%z + fad_h(3)
         
         ff(atom_h2)%x = ff(atom_h2)%x + fad_h(1)
         ff(atom_h2)%y = ff(atom_h2)%y + fad_h(2)
         ff(atom_h2)%z = ff(atom_h2)%z + fad_h(3)

       
         ff(atom_m)%x = 0.0d0; ff(atom_m)%y = 0.0d0; ff(atom_m)%z =0.0d0

      end do

    end subroutine redistrib_m_water_force

  end module mod_tip4p
