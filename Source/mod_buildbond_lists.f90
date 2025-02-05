module mod_buildbond_lists

  use global
  contains


    subroutine build_bond_list_wrapper
      implicit none

      call build_spec_bond_list_once
      if(numbond.gt.0)then
         !---build 12 bondlist
         call build_12_bond_list
         
         if(numangle.gt.0)then
            !---build 13 anglelist
            call build_13_bond_list
            
            if(numdihed.gt.0)then
               call build_14_bond_list
            endif
            
            if(numimpro.gt.0)then
               call build_impro_bond_list
            endif
         endif
      endif

      if(nshake.gt.0)then
         call build_shake_list
      endif

      call build_integrate_flag

    end subroutine build_bond_list_wrapper

    
    !---------------------------------------------------------------------!
    ! build num1bond, num2bond, and num3bond
    ! build angle list
    ! format for angle list is that localatom ---> 3rd neighbor
    ! fill angle list such that localatom.gt.global_tag
    subroutine build_spec_bond_list(moltype,particle)
      implicit none
      integer :: moltype,particle
      integer :: i,j,k
      integer :: num1tmp,num2tmp,num3tmp
      integer :: global_tag,localatom,localatom2
      integer :: type
      
      do i = 1,numMolAtom(moltype)
         
         global_tag = particle + (i-1)
         
         num1tmp   = Mol2NumBonds(moltype)%numbonds(1,i) 
         num2tmp   = Mol2NumBonds(moltype)%numbonds(2,i)
         num3tmp   = Mol2NumBonds(moltype)%numbonds(3,i)
         
         
         num1bond(global_tag)   = num1tmp 
         num2bond(global_tag)   = num2tmp+num1tmp
         num3bond(global_tag)   = num2tmp+num1tmp+num3tmp
         
         
         do j =1,num3bond(global_tag)
            localatom =  Mol2SpecTotal(moltype)%specbond(j,i)%atom+(particle-1)
            
            specbond(j,global_tag) = localatom
         enddo
      enddo
      
    end subroutine build_spec_bond_list


    !---------------------------------------------------------------------!
    ! build num1bond, num2bond, and num3bond
    ! build angle list
    ! format for angle list is that localatom ---> 3rd neighbor
    ! fill angle list such that localatom.gt.global_tag
    subroutine build_spec_bond_list_once
      implicit none
      integer :: moltype,particle,particleold
      integer :: i,j,k,l
      integer :: num1tmp,num2tmp,num3tmp
      integer :: offset,localatom,localatom2
      integer :: type

         potential = 0.0d0
    
         !---use to keep track of array index in anglelist

         particle = 0
         do i = 1,numMolType
            do j = 1, numMolArray(i)

               !---rotate other molecules
               offset = particle
               do k = 1,numMolAtom(i)
                  
                  particle = particle +1       
                  num1tmp   = Mol2NumBonds(i)%numbonds(1,k) 
                  num2tmp   = Mol2NumBonds(i)%numbonds(2,k)
                  num3tmp   = Mol2NumBonds(i)%numbonds(3,k)
                  
                  
                  num1bond(particle)   = num1tmp 
                  num2bond(particle)   = num2tmp+num1tmp
                  num3bond(particle)   = num2tmp+num1tmp+num3tmp
                  
                  
                  do l =1,num3bond(particle)
                     localatom =  Mol2SpecTotal(i)%specbond(l,k)%atom+offset
                     
                     specbond(l,particle) = localatom
                  enddo
               enddo 
            enddo
         enddo

   

         open(unit = 1003, file='specbond.dat')
         do i =1,np
            write(1003,*)'checking bonds for atom:',i,position(i)%type
            do j = 1,num1bond(i)
               write(1003,*)specbond(j,i)
            end do
            write(1003,*)'------------------------------'
            write(1003,*)
         end do


       !  open(unit = 1004, file = 'molspec.dat')
       !  do i = 1,numMolAtom(1)

       !     num1tmp   = Mol2NumBonds(1)%numbonds(1,i) 
       !     num2tmp   = Mol2NumBonds(1)%numbonds(2,i) + num1tmp
       !     num3tmp   = Mol2NumBonds(1)%numbonds(3,i) + num2tmp


       !     write(1004,*)'------------------------------'
       !     write(1004,*)'atom',i

       !     do j = 1,num1tmp
       !        write(1004,*)Mol2SpecTotal(1)%specbond(j,i)%atom
       !     end do
       !     do j = num1tmp+1,num2tmp
       !        write(1004,*)Mol2SpecTotal(1)%specbond(j,i)%atom,Mol2SpecTotal(1)%specbond(j,i)%atom2
       !     end do
       !     do j = num2tmp+1,num3tmp
       !        write(1004,*)'dihed:',Mol2SpecTotal(1)%specbond(j,i)%atom
       !        write(1004,*)'dihed:',Mol2SpecTotal(1)%specbond(j,i)%atom2
       !        write(1004,*)'dihed:',Mol2SpecTotal(1)%specbond(j,i)%atom3
       !     end do
       !     write(1004,*)'------------------------------'
       !     write(1004,*)
       !  end do
       !  close(1004)


       !  open(unit = 1004, file = 'numbond_check.dat')
       !  do i = 1,np
       !     write(1004,*),i,num1bond(i),num2bond(i),num3bond(i)
       !  enddo
       !  close(1004)

       !  open(unit = 1004, file = 'molnumbond.dat')
       !  do i = 1,numMolAtom(1)
       !     num1tmp   = Mol2NumBonds(1)%numbonds(1,i) 

       !     write(1004,*)'------------------------------'
       !     write(1004,*)'atom',i
       !     write(1004,*)Mol2NumBonds(1)%numbonds(1,i)
       !     write(1004,*)Mol2NumBonds(1)%numbonds(2,i)
       !     write(1004,*)Mol2NumBonds(1)%numbonds(3,i) 

       !     write(1004,*)'------------------------------'
       !     write(1004,*)
       !  end do
       !  close(1004)


    end subroutine build_spec_bond_list_once



    subroutine build_impro_bond_list
      integer :: i,j,k,l
      integer :: num1tmp,num2tmp,num3tmp,num4tmp
      integer :: neigh,neigh2,neigh3
      integer :: total_bond
      integer :: global
      integer :: offset
      integer :: bondtype
      
      

      global = 0
      total_bond = 0
      do i = 1,numMolType
         do j = 1, numMolArray(i)
            
            offset = global
            do k = 1,numMolAtom(i)
               
               global = global + 1
               num3tmp = Mol2NumBonds(i)%numbonds(3,k)
               num4tmp = Mol2NumBonds(i)%numbonds(4,k)

               num3tmp = num3tmp+Mol2NumBonds(i)%numbonds(2,k)
               num3tmp = num3tmp+Mol2NumBonds(i)%numbonds(1,k)
               num4tmp = num3tmp+num4tmp


               do l = num3tmp+1,num4tmp
                  neigh    = Mol2SpecTotal(i)%specbond(l,k)%atom
                  neigh2   = Mol2SpecTotal(i)%specbond(l,k)%atom2
                  neigh3   = Mol2SpecTotal(i)%specbond(l,k)%atom3
                  bondtype = Mol2SpecTotal(i)%specbond(l,k)%type

                  
                  neigh   = neigh  + offset 
                  neigh2  = neigh2 + offset
                  neigh3  = neigh3 + offset
                  

                  if(global.gt.np)then
                     print*,'FOUND ERROR WITH GLOBAL IN MOD_BUILDBOND_LISTS'
                     stop
                  endif

                  if(neigh.gt.np)then
                     print*,'FOUND ERROR WITH neigh IN MOD_BUILDBOND_LISTS'
                     stop
                  endif
                  
                  if(neigh2.gt.np)then
                     print*,'FOUND ERROR WITH neigh2 IN MOD_BUILDBOND_LISTS'
                  endif
                  
                  if(neigh3.gt.np)then
                     print*,'FOUND ERROR WITH neigh3 IN MOD_BUILDBOND_LISTS',neigh3
                  endif
                  
                     
                  total_bond = total_bond + 1
                  
                  improlist(1,total_bond) = neigh
                  improlist(2,total_bond) = neigh2
                  improlist(3,total_bond) = global
                  improlist(4,total_bond) = neigh3
                  improlist(5,total_bond) = bondtype
         
                  
               enddo
            enddo
         enddo
      enddo
      
    end subroutine build_impro_bond_list
    
    subroutine build_14_bond_list
      implicit none
      
      integer :: i,j,k,l
      integer :: num1tmp,num2tmp,num3tmp
      integer :: neigh,neigh2,neigh3
      integer :: total_bond
      integer :: global
      integer :: offset
      integer :: bondtype
      
      
      open(unit = 1000, file = 'dihed.dat')
      global = 0
      total_bond = 0
      do i = 1,numMolType
         do j = 1, numMolArray(i)
            
            offset = global
            do k = 1,numMolAtom(i)
               
               global = global + 1
               num1tmp = num1bond(global); num2tmp = num2bond(global); num3tmp = num3bond(global)
               
             
               
               do l = num2tmp+1,num3tmp
                  neigh    = Mol2SpecTotal(i)%specbond(l,k)%atom
                  neigh2   = Mol2SpecTotal(i)%specbond(l,k)%atom2
                  neigh3   = Mol2SpecTotal(i)%specbond(l,k)%atom3
                  bondtype = Mol2SpecTotal(i)%specbond(l,k)%type
                  
                 ! if(neigh3.eq.neigh)then
                 !    print*,'WE ARE THROWING ERROR 1',total_bond,neigh,neigh3
                 ! endif

                  
                  neigh  = neigh  + offset 
                  neigh2 = neigh2 + offset
                  neigh3 = neigh3 + offset
                  

                  if(global.gt.np)then
                     print*,'FOUND ERROR WITH GLOBAL'
                     STOP
                  endif

                  if(neigh.gt.np)then
                     print*,'FOUND ERROR WITH neigh'
                     stop
                  endif
                  
                  if(neigh2.gt.np)then
                     print*,'FOUND ERROR WITH neigh2'
                     stop
                  endif
                  
                  if(neigh3.gt.np)then
                     print*,'FOUND ERROR WITH neigh3',neigh3
                     stop
                  endif
                  
                  if(neigh.gt.global)then
                     
                     total_bond = total_bond + 1
                     
                     if(global.lt.neigh)then
                        dihedlist(1,total_bond) = global
                        dihedlist(2,total_bond) = neigh2
                        dihedlist(3,total_bond) = neigh3
                        dihedlist(4,total_bond) = neigh
                        dihedlist(5,total_bond) = bondtype

                        write(1000,*)global,neigh2,neigh3,neigh,bondtype
                     endif

          
                     
                  endif
               enddo
            enddo
         enddo
      enddo
      close(1000)
      
    end subroutine build_14_bond_list
    
    
    
    subroutine build_13_bond_list
      implicit none
      
      integer :: i,j,k,l
      integer :: num1tmp,num2tmp,num3tmp
      integer :: neigh,neigh2
      integer :: total_bond
      integer :: global
      integer :: offset
      integer :: bondtype
      
      global = 0
      total_bond = 0
      do i = 1,numMolType
         do j = 1, numMolArray(i)
            
            offset = global
            do k = 1,numMolAtom(i)
               
               global = global + 1
               num1tmp = num1bond(global); num2tmp = num2bond(global); num3tmp = num3bond(global)
               
               
               
               do l = num1tmp+1,num2tmp
                  neigh    = Mol2SpecTotal(i)%specbond(l,k)%atom
                  neigh2   = Mol2SpecTotal(i)%specbond(l,k)%atom2
                  bondtype = Mol2SpecTotal(i)%specbond(l,k)%type

                  
                  neigh  = neigh  + offset 
                  neigh2 = neigh2 + offset
                  
                  if(neigh.gt.global)then
                     
                     total_bond = total_bond + 1
                     
                     anglelist(1,total_bond) = global
                     anglelist(2,total_bond) = neigh2
                     anglelist(3,total_bond) = neigh
                     anglelist(4,total_bond) = bondtype
                     
                     if(global.eq.0)then
                        print*,'FOUND ERROR IN BUILD_13_BONDLIST:',i,j,k,global,neigh2,neigh
                        stop
                     elseif(neigh2.eq.0)then
                        print*,'FOUND ERROR2 IN BUILD_13_BONDLIST:',i,j,k,global,neigh2,neigh
                        stop
                     elseif(neigh.eq.0)then
                        print*,'FOUND ERROR3 IN BUILD_13_BONDLIST:',i,j,k,global,neigh2,neigh
                        stop
                     endif

                                                
                  endif
               enddo
            enddo
         enddo
      enddo
      print*,'FINAL NUMBER OF BONDS',total_bond,numangle
      print*
      
    end subroutine build_13_bond_list
    
    
    
    subroutine build_12_bond_list
      implicit none
      integer :: i,j,k,l
      integer :: num1tmp,num2tmp,num3tmp
      integer :: neigh
      integer :: total_bond
      integer :: global
      integer :: offset
      integer :: bondtype
      
      
      
      global = 0
      total_bond = 0
      do i = 1,numMolType
         do j = 1, numMolArray(i)
            
            offset = global
            do k = 1,numMolAtom(i)
               
               global = global + 1
               num1tmp = num1bond(global); num2tmp = num2bond(global); num3tmp = num3bond(global)
               
               
               
               do l = 1,num1tmp
                  neigh = Mol2SpecTotal(i)%specbond(l,k)%atom
                  bondtype = Mol2SpecTotal(i)%specbond(l,k)%type
                  
                  neigh = neigh + offset 
                  
                  if(neigh.gt.global)then
                     
                     total_bond = total_bond + 1
                     
                     bondlist(1,total_bond) = global
                     bondlist(2,total_bond) = neigh
                     bondlist(3,total_bond) = bondtype
                  endif
               enddo
            enddo
         enddo
      enddo
      
    end subroutine build_12_bond_list


    subroutine build_shake_list
      implicit none
      double precision :: bondlength,bondlength2,bondlength3
      integer :: num
      integer :: counter
      integer :: i,j,k,l,m
      integer :: atom1,atom2,atom3,atom4
      integer :: neigh
      integer :: total_bond
      integer :: global
      integer :: offset
      integer :: bondtype
      
      
      open(unit = 9999, file = 'shake.dat')
      counter = 0
      global  = 0
      do i = 1,numMolType
         do j = 1, numMolArray(i)
            do k = 1,nummolshake(i)
               
               if(i.eq.1)then
                  offset =  k
               else
                  !offset =(i-1)*nummolshake(i-1) + k
                  offset = 0
                  do m = 1,i-1
                     offset = offset+nummolshake(m)
                  enddo
                  offset = offset + k
               endif

              
               
               num   = frac_shake(offset)%num

               global = 0
               if(i.ne.1)then
                  do m = 1,i-1
                     global = global + nummolatom(m)*nummolarray(m)
                  enddo
               endif
               global = global + nummolatom(i)*(j-1)

               if(num.eq.1)then
                  
                  atom1 = frac_shake(offset)%atom1
                  atom2 = frac_shake(offset)%atom2
                  bondlength = frac_shake(offset)%bond
                  
                  
                  !---what is number of this global molecule
          

                 ! if(i.eq.1)then
                 !    global = nummolatom(i)*(j-1)
                 ! else
                 !    global = nummolatom(i)*(j-1) +nummolatom(i-1)*nummolarray(i-1)
                 ! endif
                  atom1  = atom1 + global
                  atom2  = atom2 + global
                  
                  counter = counter + 1
                  
                  shakeatom(counter)%num   = num
                  shakeatom(counter)%atom1 = atom1
                  shakeatom(counter)%atom2 = atom2
                  shakeatom(counter)%bond  = bondlength
                  
                  write(9999,*)atom1,atom2,bondlength
                  write(9999,*)frac_shake(offset)%atom1,frac_shake(offset)%atom2,offset
                  write(9999,*)k,i
                  write(9999,*)

               elseif(num.eq.-3)then
                  
                  atom1 = frac_shake(offset)%atom1
                  atom2 = frac_shake(offset)%atom2
                  atom3 = frac_shake(offset)%atom3
                  bondlength  = frac_shake(offset)%bond
                  bondlength2 = frac_shake(offset)%bond2
                  bondlength3 = frac_shake(offset)%bond3
                  
                  !---what is number of this global molecule
                 ! if(i.eq.1)then
                 !    global = nummolatom(i)*(j-1)
                 ! else
                 !    global = nummolatom(i)*(j-1) +nummolatom(i-1)*nummolarray(i-1)
                 ! endif
                  atom1  = atom1 + global
                  atom2  = atom2 + global
                  atom3  = atom3 + global
                  
                  counter = counter + 1
                  
                  shakeatom(counter)%num   = num
                  shakeatom(counter)%atom1 = atom1
                  shakeatom(counter)%atom2 = atom2
                  shakeatom(counter)%atom3 = atom3
                  shakeatom(counter)%bond   = bondlength
                  shakeatom(counter)%bond2  = bondlength2
                  shakeatom(counter)%bond3  = bondlength3
                  
                  write(9999,*)atom1,atom2,atom3,bondlength,bondlength2,bondlength3
                  write(9999,*)frac_shake(offset)%atom1,frac_shake(offset)%atom2,frac_shake(offset)%atom3,offset
                  write(9999,*)k,i
                  write(9999,*)



               elseif(num.eq.3)then
                  
                  atom1 = frac_shake(offset)%atom1
                  atom2 = frac_shake(offset)%atom2
                  atom3 = frac_shake(offset)%atom3
                  bondlength  = frac_shake(offset)%bond
                  bondlength2 = frac_shake(offset)%bond2
                  
                  !---what is number of this global molecule
                !  if(i.eq.1)then
                !     global = nummolatom(i)*(j-1)
                !  else
                !     global = nummolatom(i)*(j-1) +nummolatom(i-1)*nummolarray(i-1)
                !  endif
                  atom1  = atom1 + global
                  atom2  = atom2 + global
                  atom3  = atom3 + global
                  
                  counter = counter + 1
                  
                  shakeatom(counter)%num   = num
                  shakeatom(counter)%atom1 = atom1
                  shakeatom(counter)%atom2 = atom2
                  shakeatom(counter)%atom3 = atom3
                  shakeatom(counter)%bond   = bondlength
                  shakeatom(counter)%bond2  = bondlength2
                  
                  write(9999,*)atom1,atom2,atom3,bondlength,bondlength2
                  write(9999,*)frac_shake(offset)%atom1,frac_shake(offset)%atom2,frac_shake(offset)%atom3,offset
                  write(9999,*)k,i
                  write(9999,*)



               elseif(num.eq.4)then
                  
                  atom1 = frac_shake(offset)%atom1
                  atom2 = frac_shake(offset)%atom2
                  atom3 = frac_shake(offset)%atom3
                  atom4 = frac_shake(offset)%atom4
                  bondlength  = frac_shake(offset)%bond
                  bondlength2 = frac_shake(offset)%bond2
                  bondlength3 = frac_shake(offset)%bond3
                  
                  !---what is number of this global molecule
                !  if(i.eq.1)then
                !     global = nummolatom(i)*(j-1)
                !  else
                !     global = nummolatom(i)*(j-1) +nummolatom(i-1)*nummolarray(i-1)
                !  endif
                  atom1  = atom1 + global
                  atom2  = atom2 + global
                  atom3  = atom3 + global
                  atom4  = atom4 + global
                  
                  counter = counter + 1
                  
                  shakeatom(counter)%num   = num
                  shakeatom(counter)%atom1 = atom1
                  shakeatom(counter)%atom2 = atom2
                  shakeatom(counter)%atom3 = atom3
                  shakeatom(counter)%atom4 = atom4
                  shakeatom(counter)%bond   = bondlength
                  shakeatom(counter)%bond2  = bondlength2
                  shakeatom(counter)%bond3  = bondlength3
                  
                  write(9999,*)atom1,atom2,atom3,bondlength,bondlength2,bondlength3
                  write(9999,*)frac_shake(offset)%atom1,frac_shake(offset)%atom2,frac_shake(offset)%atom3,&
                       frac_shake(offset)%atom4,offset
                  write(9999,*)k,i
                  write(9999,*)

               endif


            enddo
         enddo
      enddo
      
    end subroutine build_shake_list


    subroutine build_integrate_flag
      implicit none
      integer :: i

      integrate_flag(:)=1
      do i = 1,np

         if(i.le.3290)then
            integrate_flag(i)= 0
         endif
      enddo
    end subroutine build_integrate_flag

    
  end module mod_buildbond_lists
