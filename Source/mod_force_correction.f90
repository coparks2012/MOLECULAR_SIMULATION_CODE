module mod_force_correction
  use global
  contains


   subroutine correct_12_bonds
      implicit none
      integer :: m,i1,i2,ioffset
      integer :: T1,T2,clock_rate,clock_max
      double precision :: force,forcecoul
      double precision :: x1,y1,z1,x2,y2,z2
      double precision :: dx,dy,dz,dr,dr2,dr2i,dri
      double precision :: grij,expm2
      double precision :: r,prefactor,erfc,t
      double precision :: boxdx,boxdy,boxdz
      integer :: i,j,l,step
      integer :: itype,jtype,neigh
      integer :: tid,num,offset
      

      call system_clock(T1,clock_rate,clock_max)
      
      do m = 1,numbond
         
         i1 = bondlist(1,m)
         i2 = bondlist(2,m)
         itype = position(i1)%type
         jtype = position(i2)%type
         
         ioffset = (itype-1)*numAtomType

         dx = position(i1)%x - position(i2)%x
         dy = position(i1)%y - position(i2)%y
         dz = position(i1)%z - position(i2)%z
         
         boxdx = dx*ibox; boxdy = dy*ibox; boxdz = dz*ibox
         boxdx = (boxdx+sign(1/(epsilon(boxdx)),boxdx)) -sign(1/epsilon(boxdx),dx)
         boxdy = (boxdy+sign(1/(epsilon(boxdy)),boxdy)) -sign(1/epsilon(boxdy),dy)
         boxdz = (boxdz+sign(1/(epsilon(boxdz)),boxdz)) -sign(1/epsilon(boxdz),dz)
         
         dx = dx-box*boxdx; dy = dy-box*boxdy; dz = dz-box*boxdz
         dr2 = dx*dx + dy*dy + dz*dz
         dr2i = 1.0/dr2

         !---electrostatic calculations     
         r = sqrt(dr2)
         prefactor =  coulpre*q(i1)*q(i2)/r
         grij = gewald*r
         expm2 = exp(-grij*grij)
         t = 1.0/(1.0 + EWALD_P*grij)
         erfc = t*(A1+t*(A2+t*(A3+t*(A4+t*A5))))*expm2
         if(dr2.gt.cut_coulsq)prefactor = 0.0d0
         forcecoul = prefactor * (erfc + EWALD_F*grij*expm2)
         forcecoul = forcecoul - (1.0-fudge_12)*prefactor
         !e_coul_host = e_coul_host + prefactor*erfc    
         !e_coul_host = e_coul_host - (1.0-fudge_12)*prefactor
         e_coul_bond = e_coul_bond + prefactor*erfc    
         e_coul_bond = e_coul_bond - (1.0-fudge_12)*prefactor

         force = forcecoul*dr2i
         
         !ffh(i1)%x = ffh(i1)%x + dx*force
         !ffh(i1)%y = ffh(i1)%y + dy*force
         !ffh(i1)%z = ffh(i1)%z + dz*force
         
         !ffh(i2)%x = ffh(i2)%x - dx*force
         !ffh(i2)%y = ffh(i2)%y - dy*force
         !ffh(i2)%z = ffh(i2)%z - dz*force

         ff_bond(i1)%x = ff_bond(i1)%x + dx*force
         ff_bond(i1)%y = ff_bond(i1)%y + dy*force
         ff_bond(i1)%z = ff_bond(i1)%z + dz*force
         
         ff_bond(i2)%x = ff_bond(i2)%x - dx*force
         ff_bond(i2)%y = ff_bond(i2)%y - dy*force
         ff_bond(i2)%z = ff_bond(i2)%z - dz*force
         
         
         v_bondx = v_bondx + dx*dx*force
         v_bondy = v_bondy + dy*dy*force
         v_bondz = v_bondz + dz*dz*force




      end do
      call system_clock(T2,clock_rate,clock_max)
      time_12corr = time_12corr + real(T2-T1)/real(clock_rate)
     ! print*,'time in 12 correction',real(T2-T1)/real(clock_rate)

      

    end subroutine correct_12_bonds

    subroutine correct_13_bonds
      implicit none
      integer :: m,i1,i2,ioffset
      integer :: T1,T2,clock_rate,clock_max
      double precision :: force,forcecoul
      double precision :: x1,y1,z1,x2,y2,z2
      double precision :: dx,dy,dz,dr,dr2,dr2i,dr6i,dr12i,dri
      double precision :: grij,expm2
      double precision :: r,prefactor,erfc,t
      double precision :: boxdx,boxdy,boxdz
      integer :: i,j,l,step
      integer :: itype,jtype,neigh
      integer :: tid,num,offset
    
      call system_clock(T1,clock_rate,clock_max)
      
      do m = 1,numangle
         
         i1 = anglelist(1,m)
         i2 = anglelist(3,m)
         itype = position(i1)%type
         jtype = position(i2)%type
         
         ioffset = (itype-1)*numAtomType
         
         dx = position(i1)%x - position(i2)%x
         dy = position(i1)%y - position(i2)%y
         dz = position(i1)%z - position(i2)%z
         
         boxdx = dx*ibox; boxdy = dy*ibox; boxdz = dz*ibox
         boxdx = (boxdx+sign(1/(epsilon(boxdx)),boxdx)) -sign(1/epsilon(boxdx),dx)
         boxdy = (boxdy+sign(1/(epsilon(boxdy)),boxdy)) -sign(1/epsilon(boxdy),dy)
         boxdz = (boxdz+sign(1/(epsilon(boxdz)),boxdz)) -sign(1/epsilon(boxdz),dz)
         
         dx = dx-box*boxdx; dy = dy-box*boxdy; dz = dz-box*boxdz
         dr2 = dx*dx + dy*dy + dz*dz
         dr2i = 1.0/dr2
         
         !---electrostatic calculations     
         r = sqrt(dr2)
         prefactor =  coulpre*q(i1)*q(i2)/r
         grij = gewald*r
         expm2 = exp(-grij*grij)
         t = 1.0/(1.0 + EWALD_P*grij)
         erfc = t*(A1+t*(A2+t*(A3+t*(A4+t*A5))))*expm2
         if(dr2.gt.cut_coulsq)prefactor = 0.0d0
         forcecoul = prefactor * (erfc + EWALD_F*grij*expm2)
         forcecoul = forcecoul - (1.0-fudge_13)*prefactor

         e_coul_angle = e_coul_angle + prefactor*erfc    
         e_coul_angle = e_coul_angle - (1.0-fudge_13)*prefactor



         !e_coul_host = e_coul_host + prefactor*erfc    
         !e_coul_host = e_coul_host - (1.0-fudge_13)*prefactor
         
         force = forcecoul*dr2i
         
         !ffh(i1)%x = ffh(i1)%x + dx*force
         !ffh(i1)%y = ffh(i1)%y + dy*force
         !ffh(i1)%z = ffh(i1)%z + dz*force
         
         !ffh(i2)%x = ffh(i2)%x - dx*force
         !ffh(i2)%y = ffh(i2)%y - dy*force
         !ffh(i2)%z = ffh(i2)%z - dz*force


         ff_angle(i1)%x = ff_angle(i1)%x + dx*force
         ff_angle(i1)%y = ff_angle(i1)%y + dy*force
         ff_angle(i1)%z = ff_angle(i1)%z + dz*force
         
         ff_angle(i2)%x = ff_angle(i2)%x - dx*force
         ff_angle(i2)%y = ff_angle(i2)%y - dy*force
         ff_angle(i2)%z = ff_angle(i2)%z - dz*force

         v_anglex = v_anglex + dx*dx*force
         v_angley = v_angley + dy*dy*force
         v_anglez = v_anglez + dz*dz*force

         
      end do

      call system_clock(T2,clock_rate,clock_max)
      time_13corr = time_13corr + real(T2-T1)/real(clock_rate)

    !  print*,'ela[sed time 13 correc',real(T2-T1)/real(clock_rate)


    end subroutine correct_13_bonds
    
    
    subroutine correct_14_bonds
      implicit none
      integer :: m,i1,i2,ioffset
      integer :: T1,T2,clock_rate,clock_max
      double precision :: force,forcecoul
      double precision :: x1,y1,z1,x2,y2,z2
      double precision :: dx,dy,dz,dr,dr2,dr2i,dr6i,dr12i,dri
      double precision :: grij,expm2
      double precision :: qtmp,r,prefactor,erfc,t
      double precision :: boxdx,boxdy,boxdz
      integer :: i,j,l,step
      integer :: itype,jtype,neigh
      integer :: tid,num,offset
      

      call system_clock(T1,clock_rate,clock_max)

      do m = 1,numdihed
         
         i1 = dihedlist(1,m)
         i2 = dihedlist(4,m)
         itype = position(i1)%type
         jtype = position(i2)%type
         
         ioffset = (itype-1)*numAtomType
         
         dx = position(i1)%x - position(i2)%x
         dy = position(i1)%y - position(i2)%y
         dz = position(i1)%z - position(i2)%z
         
         boxdx = dx*ibox; boxdy = dy*ibox; boxdz = dz*ibox
         boxdx = (boxdx+sign(1/(epsilon(boxdx)),boxdx)) -sign(1/epsilon(boxdx),dx)
         boxdy = (boxdy+sign(1/(epsilon(boxdy)),boxdy)) -sign(1/epsilon(boxdy),dy)
         boxdz = (boxdz+sign(1/(epsilon(boxdz)),boxdz)) -sign(1/epsilon(boxdz),dz)
         
         dx = dx-box*boxdx; dy = dy-box*boxdy; dz = dz-box*boxdz
         dr2 = dx*dx + dy*dy + dz*dz
         dr2i = 1.0/dr2



         !---electrostatic calculations     
         r = sqrt(dr2)
         prefactor =  coulpre*q(i1)*q(i2)/r
         grij = gewald*r
         expm2 = exp(-grij*grij)
         t = 1.0/(1.0 + EWALD_P*grij)
         erfc = t*(A1+t*(A2+t*(A3+t*(A4+t*A5))))*expm2
         if(dr2.gt.cut_coulsq)prefactor = 0.0d0
         forcecoul = prefactor * (erfc + EWALD_F*grij*expm2)
         forcecoul = forcecoul - (1.0-fudge_14_coul)*prefactor
         

         e_coul_dihed = e_coul_dihed + prefactor*erfc    
         e_coul_dihed = e_coul_dihed - (1.0-fudge_14_coul)*prefactor
         !e_coul_host = e_coul_host + prefactor*erfc    
         !e_coul_host = e_coul_host - (1.0-fudge_14_coul)*prefactor
         
         
         force = forcecoul*dr2i
         
         
         !ffh(i1)%x = ffh(i1)%x + dx*force
         !ffh(i1)%y = ffh(i1)%y + dy*force
         !ffh(i1)%z = ffh(i1)%z + dz*force
         
         !ffh(i2)%x = ffh(i2)%x - dx*force
         !ffh(i2)%y = ffh(i2)%y - dy*force
         !ffh(i2)%z = ffh(i2)%z - dz*force

         ff_dihed(i1)%x = ff_dihed(i1)%x + dx*force
         ff_dihed(i1)%y = ff_dihed(i1)%y + dy*force
         ff_dihed(i1)%z = ff_dihed(i1)%z + dz*force
         
         ff_dihed(i2)%x = ff_dihed(i2)%x - dx*force
         ff_dihed(i2)%y = ff_dihed(i2)%y - dy*force
         ff_dihed(i2)%z = ff_dihed(i2)%z - dz*force

         v_dihedx = v_dihedx + dx*dx*force
         v_dihedy = v_dihedy + dy*dy*force
         v_dihedz = v_dihedz + dz*dz*force

         
      end do
      call system_clock(T2,clock_rate,clock_max)
      time_14corr = time_14corr + real(T2-T1)/real(clock_rate)

    !  print*,'elapsed time 14 corr',real(T2-T1)/real(clock_rate)



    end subroutine correct_14_bonds
    

end module mod_force_correction
