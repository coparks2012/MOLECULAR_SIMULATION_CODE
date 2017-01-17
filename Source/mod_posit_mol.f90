module mod_posit_mol

  use global
  use mod_neighbor
  use mod_pbc
  use mod_buildbond_lists
  use mod_minimg
  implicit none

contains


  subroutine random_mol
    implicit none
    double precision :: x1,y1,z1,x,y,z
    double precision :: r11,r12,r13,r21,r22,r23,r31,r32,r33
    double precision :: dx,dy,dz
    integer :: i,j,k
    integer :: particle
    
 
     do i = 1,numMolType
       do j = 1, numMolArray(i)

          !---pick random displacement vector within box
          x1 = box*ran(); y1 = box*ran(); z1 = box*ran()

          
          !---put first particle of the molecule at displacement
          particle = particle + 1
          position(particle)%x = x1; position(particle)%y = y1; position(particle)%z = z1
                    
          
          call build_rotation_matrix(r11,r12,r13,r21,r22,r23,r31,r32,r33)
          
          
          !---rotate other molecules
          do k = 2,numMolAtom(i)
             
             particle = particle +1
             x = FC(i)%coords(k)%x; y = FC(i)%coords(k)%y; z = FC(i)%coords(k)%z
             
         !    position(particle)%x = x*r11 +y*r12 +z*r13 +x1
         !    position(particle)%y = x*r21 +y*r22 +z*r23 +y1
         !    position(particle)%z = x*r31 +y*r32 +z*r33 +z1   


             position(particle)%x = x +x1
             position(particle)%y = y +y1
             position(particle)%z = z +z1 

          enddo

       enddo
    enddo


  end subroutine random_mol


  subroutine posit_LJ
    implicit none
    double precision :: x1,y1,z1
    double precision :: pot0
    double precision :: delta,pref,rand1
    integer :: i,j,k
    integer :: particle,particleold
    integer :: accept
    integer :: T1,T2,clock_rate,clock_max

    call system_clock(T1,clock_rate,clock_max)
    potential = 0.0d0

    position(1)%x = box*ran(); position(1)%y = box*ran(); position(1)%z = box*ran()

    !--- global particle tag
    particle = 1
    do i = 1,np
       accept = 0

       do while(accept.eq.0.and.i.le.np)
          
          particleold = particle
          pot0 = potential
                
          
          position(particle)%x = box*ran(); position(particle)%y = box*ran(); position(particle)%z = box*ran()

          !position(particle)%x = position(particle-1)%x + 2.0d0*(ran()-0.50d0)
          !position(particle)%y = position(particle-1)%y + 2.0d0*(ran()-0.50d0)
          !position(particle)%z = position(particle-1)%z + 2.0d0*(ran()-0.50d0)



          call pbc_particle(particle)
          call bin(particle)
          call linked_cell_energy

          delta = potential - pot0

          pref = exp(-delta/(kboltz*temp))
          
          rand1 = ran()

          if(delta.lt.0.0d0)then
             potential = potential + delta
             accept    = 1

          elseif(pref.gt.rand1)then
             potential = potential + delta
             accept    = 1
          else
             potential = pot0
             accept    = 0
             particle  = particleold
          endif

       end do
       print*,'done with particle:',i,potential
       
    enddo   
    call system_clock(T2,clock_rate,clock_max)


    PRINT*,'elapsed time placing particles',real(T2-T1)/real(clock_rate)
  end subroutine posit_LJ



  
  subroutine posit_mol
    implicit none
    double precision :: x1,y1,z1
    double precision :: r11,r12,r13,r21,r22,r23,r31,r32,r33
    double precision :: x,y,z
    double precision :: dx,dy,dz,dr2
    double precision :: pot0
    double precision :: delta,pref,rand1
    integer :: i,j,k,l
    integer :: particle,particleold
    integer :: accept
    integer :: angle_flag,angle_flag0
    integer :: T1,T2,clock_rate,clock_max
    
    call system_clock(T1,clock_rate,clock_max)
    

    potential = 0.0d0
    
    !---use to keep track of array index in anglelist
    angle_flag = 0

    particle = 1
    do i = 1,numMolType
       do j = 1, numMolArray(i)

             accept = 0
             angle_flag0 = angle_flag

             do while(accept.eq.0)
                particleold = particle
                pot0 = potential
                
         
                
                !---pick random displacement vector within box
                x1 = box*ran(); y1 = box*ran(); z1 = box*ran()
                
                
                !---put first particle of the molecule at displacement
                position(particle)%x = x1; position(particle)%y = y1; position(particle)%z = z1
                
                call build_spec_bond_list(i,particle)


                call build_rotation_matrix(r11,r12,r13,r21,r22,r23,r31,r32,r33)


                !---rotate other molecules
                do k = 2,numMolAtom(i)
                  
                  particle = particle +1
                   x = FC(i)%coords(k)%x; y = FC(i)%coords(k)%y; z = FC(i)%coords(k)%z
                   
                   position(particle)%x = x*r11 +y*r12 +z*r13 +x1
                   position(particle)%y = x*r21 +y*r22 +z*r23 +y1
                   position(particle)%z = x*r31 +y*r32 +z*r33 +z1


        
                enddo

                
                !---apply PBC to all new particles
                do l = particleold,particle
                   call pbc_particle(l)
                enddo

                !call check_geometry(particleold,particle)

                call bin(particle)
                call linked_cell_energy
                
                delta = potential - pot0
                
                pref = exp(-delta/(kboltz*temp))
                
                rand1 = ran()
                
                if(delta.lt.0.0d0)then
                   potential = potential + delta
                   accept    = 1
                   particle  = particle + 1
                   
                elseif(pref.gt.rand1)then
                   potential = potential + delta
                   accept    = 1
                   particle  = particle + 1
                else
                   potential = pot0
                   accept    = 0
                   particle  = particleold
                   angle_flag = angle_flag0
                endif
             end do
             !print*,'done with particle',i,particle,potential
          enddo
       enddo
       print*,'FINAL POTENTIAL FROM PLACING MOLECULES:',potential,rcut

       call bin(np)
       call linked_cell_energy
       print*,'FINAL POTENTIAL FROM PLACING MOLECULES:',potential


       call system_clock(T2,clock_rate,clock_max)

       print*,'ELAPSED TIME PLACING MOLECULES:',real(T2-T1)/real(clock_rate)
  end subroutine posit_mol


  

         
          
  subroutine build_rotation_matrix(r11,r12,r13,r21,r22,r23,r31,r32,r33)
    implicit none
    double precision :: r11,r12,r13,r21,r22,r23,r31,r32,r33
    double precision :: q0,q1,q2,q3,qnorm,qnormi

    !---pick random quaternion
    q0 = ran(); q1 = ran(); q2 = ran(); q3 = ran()
    qnorm = q0*q0+q1*q1+q2*q2+q3*q3
    qnormi = 1.0d0/sqrt(qnorm)


    
    q0 = q0*qnormi; q1 = q1*qnormi; q2 = q2*qnormi; q3 = q3*qnormi
            
    !---build rotation matrix
    r11 = q0*q0+q1*q1-q2*q2-q3*q3; r12 = 2.0d0*(q1*q2-q0*q3);     r13 = 2.0d0*(q1*q3+q0*q2)
    r21 = 2.0d0*(q1*q2+q0*q3);     r22 = q0*q0-q1*q1+q2*q2-q3*q3; r23 = 2.0d0*(q2*q3-q0*q1)
    r31 = 2.0d0*(q1*q3-q0*q2);     r32 = 2.0d0*(q2*q3+q0*q1);     r33 = q0*q0-q1*q1-q2*q2+q3*q3


  end subroutine build_rotation_matrix


  subroutine bin(particle)
    implicit none
    integer :: particle
    integer :: i,icell,step
    
    !    -- initialize head of chain
    HOC(:) = 0
   

    !    -- make linked list
    do i=1,particle
       
       !    -- determine cell number
       icell = xyz2bin(position(i)%x,position(i)%y,position(i)%z)
       
       !    -- update linked-list and head of chain
       LL(i) = HOC(icell)
       HOC(icell) = i
    enddo
    
  end subroutine bin

  subroutine linked_cell_energy
    implicit none
    double precision :: x1,y1,z1
    double precision :: x2,y2,z2
    double precision :: dx,dy,dz
    double precision :: dr2,dr2i,dr6i,dr12i
    double precision :: fudge_factor
    integer :: i,j,k,l
    integer :: particle1,particle2
    integer :: cell_neigh
    integer :: num1tmp,num2tmp,num3tmp


   

    potential = 0.0d0

    do i = 0,ncellT-1
       particle1 = HOC(i)
       
       do while(particle1.ne.0)
          
          x1 = position(particle1)%x
          y1 = position(particle1)%y
          z1 = position(particle1)%z

          num1tmp = num1bond(particle1)
          num2tmp = num2bond(particle1)
          num3tmp = num3bond(particle1)

          !---loop over neighboring cells
          do j = i*13,13*i+12
             cell_neigh = cnum(j)

             particle2 = HOC(cell_neigh)
             do while(particle2.ne.0)
                
                x2 = position(particle2)%x
                y2 = position(particle2)%y
                z2 = position(particle2)%z
                
                dx = x2-x1; dy = y2-y1; dz = z2-z1
                
                dx = dx-box*nint(ibox*dx)
                dy = dy-box*nint(ibox*dy)
                dz = dz-box*nint(ibox*dz)
                
                dr2 = dx*dx + dy*dy + dz*dz
                
                if(dr2.lt.rcut2)then
                   
                   
                   fudge_factor = 1.0d0
                   do l = 1,num3tmp
                      if(specbond(l,particle1).eq.particle2)then
                         fudge_factor = 0.0d0
                      endif
                   enddo
                   
                   
                   dr2i = 1.0d0/dr2
                   dr6i = dr2i*dr2i*dr2i
                   potential = potential + 4.0d0*fudge_factor*dr6i*(dr6i-1.0d0)
                endif
                
                particle2 = ll(particle2)
             enddo
          enddo

          particle2 = HOC(i)
          do while(particle2.ne.0)

             if(particle2.gt.particle1)then
                x2 = position(particle2)%x
                y2 = position(particle2)%y
                z2 = position(particle2)%z
                
                dx = x2-x1; dy = y2-y1; dz = z2-z1
                
                dx = dx-box*nint(ibox*dx)
                dy = dy-box*nint(ibox*dy)
                dz = dz-box*nint(ibox*dz)
                
                dr2 = dx*dx + dy*dy + dz*dz
                
                if(dr2.lt.rcut2)then
                   
                   
                   fudge_factor = 1.0d0
                   do l = 1,num3tmp
                      if(specbond(l,particle1).eq.particle2)then
                         fudge_factor = 0.0d0
                      endif
                   enddo
                   
                   
                   dr2i = 1.0d0/dr2
                   dr6i = dr2i*dr2i*dr2i
                   potential = potential + 4.0d0*fudge_factor*dr6i*(dr6i-1.0d0)
                endif
             endif
             particle2 = ll(particle2)
          enddo
          particle1 = ll(particle1)
       enddo
    enddo

  end subroutine linked_cell_energy
end module mod_posit_mol
