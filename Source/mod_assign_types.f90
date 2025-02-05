module mod_assign_types
  use global
  implicit none

  contains

    subroutine assign_atom_types
      implicit none
      integer :: tag
      integer :: i,j,k
      
      tag = 0
      do i = 1,numMolType
         do j = 1, numMolArray(i)
            do k = 1,numMolAtom(i)
               
               tag = tag+1
               position(tag)%type = FC(i)%coords(k)%type
               
            enddo
         enddo
      enddo
    end subroutine assign_atom_types
    
    subroutine assign_charge_types
      integer :: tag
      integer :: i,j,k
      double precision :: net_c
      
      net_c = 0.0d0
      do i =1 ,np
         tag = position(i)%type
         q(i) = qtype(tag)
         net_c = net_c + q(i)
      enddo
      print*,'WHAT IS NET CHARGE OF SYSTEM',net_c

      if(coul_flag.eq.2)then
         e_self = 0.0

         do i= 1,np
            e_self = e_self - coulpre*q(i)*q(i)*gewald/sqrt(acos(-1.0d0))
         enddo
         print*
         print*,'self energy check'
         print*,'gewald',gewald
         print*,'self',e_self
         print*
         
      elseif(coul_flag.eq.1)then
         e_self = 0.0
         
         print*
         print*,'self energy check'
         print*,'e_shift',e_shift
         print*,'alpha',alpha
         print*,'MY_PIS',MY_PIS
         PRINT*
         do i= 1,np
            e_self = e_self - (0.5*e_shift+alpha/MY_PIS)*&
                 q(i)*q(i)*coulpre
         enddo
      endif


    end subroutine assign_charge_types


    subroutine build_mol_list
      implicit none
      integer :: a,b,c
      integer :: num1,num2
      integer :: spot

      mol(:,:) = 0

      !---currently assigning all non hydrogen atoms in TTA as tracked
      !---need to think about the degeneracy of the oxygen atoms
      do a = 1,numMolArray(1)
         do b = 1,numatns
            spot = (a-1)*10 + b
            mol(1,spot) = 1
            mol(2,spot) = (a-1)*numatns + b
         enddo
         

      enddo

    end subroutine build_mol_list
    
  end module mod_assign_types
