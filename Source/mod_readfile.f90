module mod_readfile

  use global
  use mod_memory
  contains

    !-------------------------------------------------------------------!
    ! wrapped routine for reading all input files
    subroutine read_files_wrapper
      implicit none
      call read_input
      call read_gro_input
      call read_topol
      call check_gro_file
      call check_topol_file

    end subroutine read_files_wrapper
  
    !-----------------------------------------------!
    ! process input  file

    subroutine read_input
      implicit none
      print*,'we are in read input'
      open(unit=1,file='input_files/input')
      
      read(1,*)
      read(1,'(a)')ensemble
      read(1,*)
      read(1,*)

      read(1,'(a)')STA
      read(1,*)
      read(1,*)

      read(1,*)file_type_flag
      read(1,*)
      read(1,*)
      
      read(1,*)filename_mol
      read(1,*)
      read(1,*)
      filename_mol = adjustl(filename_mol)
      filename_mol = trim(filename_mol)

      read(1,*)filename_topol
      read(1,*)
      read(1,*)
      filename_topol=adjustl(filename_topol)
      filename_topol=trim(filename_topol)

      if(file_type_flag.eq.0)then
         read(1,*)filename_lammps
         filename_lammps = 'input_files/'//trim(filename_lammps)

         filename_lammps = adjustl(filename_lammps)
         filename_lammps = trim(filename_lammps)
      elseif(file_type_flag.eq.1)then
         print*,'we are reading in filename_gro'
         read(1,*)filename_gro
         filename_gro = 'input_files/'//trim(filename_gro)
         filename_gro = adjustl(filename_gro)
         filename_gro = trim(filename_gro)
      elseif(file_type_flag.eq.2)then
         print*,'we are reading in cluster restart files:file type flag is 2'
         read(1,*)filename_restart
         filename_restart = 'input_files/'//trim(filename_restart)
         filename_restart = adjustl(filename_restart)
         filename_restart = trim(filename_restart)
      endif
      read(1,*)
      read(1,*)
      write(1234,*)'=================================================='
      write(1234,*)' read in files'
      write(1234,*)trim(filename_topol)
      write(1234,*)trim(filename_mol)
      if(file_type_flag.eq.1)then
         write(1234,*)trim(filename_gro)
      endif

      read(1,*)vlist_flag
      read(1,*)
      read(1,*)

      read(1,*)intel_flag
      read(1,*)
      read(1,*)

      read(1,*)coul_flag
      read(1,*)
      read(1,*)

      read(1,*)ngrid_x,ngrid_y,ngrid_z
      read(1,*)
      read(1,*)
      read(1,*)fudge_12,fudge_13,fudge_14,fudge_14_coul
      read(1,*)
      read(1,*)
      !---store in higher precision for amber
      fudge_14_coul = 1.0d0/1.20d0
      
      read(1,*)newton3rd_flag
      read(1,*)
      read(1,*)

      read(1,*)dt
      read(1,*)
      read(1,*)

      read(1,*)rcut,cut_coul,neighbor_pad,mrcut,denscut,densmin
      read(1,*)
      read(1,*)

      read(1,*)cluster_type
      read(1,*)
      read(1,*)

      read(1,*)dens             ! thermodynamical conditions
      read(1,*)
      read(1,*)

      read(1,*)temp
      read(1,*)
      read(1,*)

      read(1,*)ptarget                ! read in presssure as atm
      read(1,*)
      read(1,*)

      read(1,*)pdamp
      read(1,*)
      read(1,*)

      read(1,*)alpha
      read(1,*)
      read(1,*)

      read(1,*)mdeq,mdsteps
      read(1,*)
      read(1,*)

      read(1,*)nsample
      read(1,*)
      read(1,*)

      read(1,*)nadjust
      read(1,*)
      read(1,*)

      read(1,*)numSweep
      read(1,*)
      read(1,*)

      read(1,*)numStep
      read(1,*)
      read(1,*)

      read(1,*)numSweepEquil
      read(1,*)
      read(1,*)

      read(1,*)kn
      read(1,*)
      read(1,*)

      read(1,*)anneal_temp
      close(1)

      if(rcut.gt.cut_coul)then
         print*,'vdw cutoff is less than coul cut'
      endif

      if(neighbor_pad.lt.1.0d0)then
         print*,'ERROR NEIGHBOR PAD IS LESS THAN ONE:',neighbor_pad
      endif

      print*,'succesfully read in input file',rank
      write(1234,*)'succesfully read in input file',rank
    end subroutine read_input




    !-----------------------------------------------!
    ! process input gro file
    ! processes fractional coordinates from input file
    ! assign atom types to each atom in position array
    subroutine read_gro_input
      implicit none
      integer :: i,j,k,tag,dum
      integer :: numatom
      


      open(unit = 1, file = 'input_files/'//trim(filename_mol))


      read(1,*)
      
      read(1,*)numMolType
     
      !---allocate moltype arrays
      call moltype_memory
      
      do i = 1,nummoltype
         
         read(1,*)
         read(1,*)numMolArray(i)
         read(1,*)numMolAtom(i)
         
         !---allocate fractional coordinate arrays
         !---allocate mol2numBonds array
         call moltype_datatype_memory(i)
         
         do j = 1,numMolAtom(i)
            read(1,*)FC(i)%coords(j)%type,dum,FC(i)%coords(j)%x,FC(i)%coords(j)%y,FC(i)%coords(j)%z
         enddo
         
      enddo
      
      
      close(1)
      

    end subroutine read_gro_input
    
    subroutine read_topol
      implicit none
      double precision :: local_eps
      double precision :: k1,theta1,k2,r1,theta0
      double precision :: KK1,KK2,KK3,KK4
      double precision :: bondlength,bondlength2,bondlength3
      integer :: i,j,k,D,tag
      integer :: track1,track2
      integer :: moltype,bondtype,angletype,dihedtype
      integer :: atom1,atom2,atom3,atom4
      integer :: offset
      integer :: off1,off2
      integer :: loopflag
      integer :: loop,loop2
      integer :: tag1,tag2,tag3
      integer :: numbond
      integer :: numangle
      integer :: intid
      integer :: num1tmp,total
      integer :: numshake1,numshake2,numincons
      double precision :: local_k
      integer :: sh_one,mult_one
      character(len=150) :: id,rank_char,filename

      print*,'LETS READ TOPOL'
      !--- numbonds(1,*) ---> 1-2 bonds
      !--- numbonds(2,*) ---> 1-3 bonds
      !--- numbonds(3,*) ---> 1-4 bonds
      !--- numbonds(4,*) ---> improper
      !---zero numbonds for each molecule
      do i = 1,numMolType
         do j = 1,numMolAtom(i)
            do k = 1,4
               Mol2NumBonds(i)%numbonds(k,j) = 0
            enddo
         enddo
      end do

      open(unit = 1, file = 'input_files/'//trim(filename_topol))
   
      !---read in nonbonded interaction parameters
      read(1,*)
      read(1,*)numAtomType
      
      
      call atomtype_memory
      
      read(1,*)


      
      do i = 1,numAtomType
  

         read(1,*)tag,ltype2ftype(i),type2string(i) ,sig(i),eps(i),qtype(i),mass(i) 

         type2string(i) = adjustl(type2string(i))
         type2string(i) = trim(type2string(i))
         
      enddo
    

      !---read in bonded interaction parameters
      read(1,*)
      read(1,*)
      read(1,*)numbondtype,loopflag
      read(1,*)
      
      if(numbondtype.gt.0)then
         call bondtype_memory
      endif

      do i = 1,loopflag
         read(1,*)moltype,bondtype,atom1,atom2,k1,r1

         
         bondcoeff(1,bondtype) = k1
         bondcoeff(2,bondtype) = r1

         !---UPDATE NUMBER OF BONDS FOR EACH ATOM
         track1 = Mol2NumBonds(moltype)%numbonds(1,atom1) + 1
         track2 = Mol2NumBonds(moltype)%numbonds(1,atom2) + 1
         Mol2NumBonds(moltype)%numbonds(1,atom1) = track1
         Mol2NumBonds(moltype)%numbonds(1,atom2) = track2

         !---TRACK WHICH ATOMS ARE BONDED
         Mol2Spectotal(moltype)%specbond(track1,atom1)%atom = atom2
         Mol2SpecTotal(moltype)%specbond(track1,atom1)%type = bondtype

         Mol2Spectotal(moltype)%specbond(track2,atom2)%atom = atom1
         Mol2SpecTotal(moltype)%specbond(track2,atom2)%type = bondtype


      end do
   
      !--------------------------------------------------------------------------
      !--------------------------------------------------------------------------
      !---read in angle bending parameters
      read(1,*)
      read(1,*)
      read(1,*)numangletype,loopflag
      read(1,*)

      if(numangletype.gt.0)then
         call angletype_memory
      endif


      do i = 1, loopflag
         read(1,*)moltype,angletype,atom1,atom2,atom3,k1,theta1

         !---UPDATE NUMBER OF BONDS FOR EACH ATOM
         track1 = Mol2NumBonds(moltype)%numbonds(2,atom1) + 1
         track2 = Mol2NumBonds(moltype)%numbonds(2,atom3) + 1
         
         Mol2NumBonds(moltype)%numbonds(2,atom1) = track1
         Mol2NumBonds(moltype)%numbonds(2,atom3) = track2

         track1 = Mol2NumBonds(moltype)%numbonds(1,atom1) + track1
         track2 = Mol2NumBonds(moltype)%numbonds(1,atom3) + track2

         !---TRACK WHICH ATOMS ARE BONDED
         Mol2Spectotal(moltype)%specbond(track1,atom1)%atom  = atom3
         Mol2Spectotal(moltype)%specbond(track1,atom1)%atom2 = atom2
         Mol2Spectotal(moltype)%specbond(track1,atom1)%type  = angletype

         Mol2Spectotal(moltype)%specbond(track2,atom3)%atom  = atom1
         Mol2Spectotal(moltype)%specbond(track2,atom3)%atom2 = atom2
         Mol2Spectotal(moltype)%specbond(track2,atom3)%type  = angletype


         !---store parameters
         anglecoeff(1,angletype)  = k1
         anglecoeff(2,angletype)  = theta1*(2.0d0*acos(-1.0d0)/(360.0d0))

      end do


      !--------------------------------------------------------------
      !-------------------------------------------------------------
      !---allocate 1-4 bond type arrays/ read in 1-4 bond parameters
      read(1,*)
      read(1,*)
      read(1,*)numdihedtype,loopflag
      read(1,*)

      if(numangletype.gt.0)then
         call dihedtype_memory
      endif

      do i = 1, loopflag
         read(1,*)moltype,dihedtype,atom1,atom2,atom3,atom4,local_k,mult_one,sh_one

         !---UPDATE NUMBER OF BONDS FOR EACH ATOM
         track1 = Mol2NumBonds(moltype)%numbonds(3,atom1) + 1
         track2 = Mol2NumBonds(moltype)%numbonds(3,atom4) + 1
         
         Mol2NumBonds(moltype)%numbonds(3,atom1) = track1
         Mol2NumBonds(moltype)%numbonds(3,atom4) = track2

         track1 = track1 +  Mol2NumBonds(moltype)%numbonds(1,atom1) + Mol2NumBonds(moltype)%numbonds(2,atom1) 
         track2 = track2 +  Mol2NumBonds(moltype)%numbonds(1,atom4) + Mol2NumBonds(moltype)%numbonds(2,atom4) 


         !---TRACK WHICH ATOMS ARE BONDED
         Mol2Spectotal(moltype)%specbond(track1,atom1)%atom  = atom4
         Mol2Spectotal(moltype)%specbond(track1,atom1)%atom2 = atom2
         Mol2Spectotal(moltype)%specbond(track1,atom1)%atom3 = atom3
         Mol2Spectotal(moltype)%specbond(track1,atom1)%type  = dihedtype

         Mol2Spectotal(moltype)%specbond(track2,atom4)%atom  = atom1
         Mol2Spectotal(moltype)%specbond(track2,atom4)%atom2 = atom2
         Mol2Spectotal(moltype)%specbond(track2,atom4)%atom3 = atom3
         Mol2Spectotal(moltype)%specbond(track2,atom4)%type  = dihedtype


         !---store interaction parameters
      
         !---need to multiply k by 1/2
         k_dihed(dihedtype) = local_k
         shift_dihed(dihedtype) = sh_one
         cos_shift(dihedtype) = cos(acos(-1.0d0)*sh_one/180.0)
         sin_shift(dihedtype) = sin(acos(-1.0d0)*sh_one/180.0)
         multiplicity(dihedtype) = mult_one


      end do
      

      !--------------------------------------------------------------
      !-------------------------------------------------------------
      !---allocate improper 1-4 bond type arrays/ read in 1-4 bond parameters
      read(1,*)
      read(1,*)
      read(1,*)numimprotype,loopflag
      read(1,*)
      
      if(numimprotype.gt.0)then
         call improdihedtype_memory
      endif

      do i = 1,loopflag
         read(1,*)moltype,bondtype,atom1,atom2,atom3,atom4,k1,D,mult_one

         !---CONVERT ANGLE TO RADIANS
         improcoeff(1,bondtype) = k1
         improcoeff(2,bondtype) = theta0*(acos(-1.0d0)/180.0d0)

         !---UPDATE NUMBER OF BONDS FOR EACH ATOM
         track1 = Mol2NumBonds(moltype)%numbonds(4,atom3) + 1
         !track2 = Mol2NumBonds(moltype)%numbonds(4,atom2) + 1

         Mol2NumBonds(moltype)%numbonds(4,atom3) = track1
         !Mol2NumBonds(moltype)%numbonds(4,atom2) = track2


         track1 = track1 + Mol2NumBonds(moltype)%numbonds(3,atom3) 
         track1 = track1 + Mol2NumBonds(moltype)%numbonds(2,atom3) 
         track1 = track1 + Mol2NumBonds(moltype)%numbonds(1,atom3) 


         !---TRACK WHICH ATOMS ARE BONDED
         Mol2Spectotal(moltype)%specbond(track1,atom3)%atom  = atom1
         Mol2Spectotal(moltype)%specbond(track1,atom3)%atom2 = atom2
         Mol2Spectotal(moltype)%specbond(track1,atom3)%atom3 = atom4
         Mol2SpecTotal(moltype)%specbond(track1,atom3)%type  = bondtype
         
         k_impro(bondtype) = k1
         sign_impro(bondtype) = D
         multiplicity_impro(bondtype) = mult_one


      end do
      



      !--------------------------------------------------------------
      !-------------------------------------------------------------
      !---shake constraints
      read(1,*)
      read(1,*)
      read(1,*)
      total = 0
      nshake = 0
      do i = 1,nummoltype
         read(1,*)nummolshake(i)
         total = total + nummolshake(i)
         nshake = nummolarray(i)*nummolshake(i)
      enddo
      read(1,*)
      read(1,*)
      allocate(frac_shake(total))

      
      ncons = 0
      do i = 1,total
         read(1,*)numincons
         
         if(numincons.eq.1)then
            frac_shake(i)%num = numincons
            read(1,*)moltype,atom1,atom2,bondlength
            frac_shake(i)%atom1 = atom1
            frac_shake(i)%atom2 = atom2
            frac_shake(i)%bond = bondlength
            
            ncons = ncons + nummolarray(moltype)*1

          !  if(i.le.nummolshake(1))then
          !     ncons = ncons + nummolarray(1)*1
          !  else
          !     ncons = ncons + nummolarray(2)*1
          !  endif

            !--CH2 group
         elseif(numincons.eq.3)then
            frac_shake(i)%num = numincons
            read(1,*)moltype,atom1,atom2,atom3,bondlength,bondlength2
            frac_shake(i)%atom1 = atom1
            frac_shake(i)%atom2 = atom2
            frac_shake(i)%atom3 = atom3
            frac_shake(i)%bond = bondlength
            frac_shake(i)%bond2 = bondlength2

            ncons = ncons + nummolarray(moltype)*2

          !  if(i.le.nummolshake(1))then
          !     ncons = ncons + nummolarray(1)*2
          !  else
          !     ncons = ncons + nummolarray(2)*2
          !  endif

            !---CH3 group
         elseif(numincons.eq.4)then
            frac_shake(i)%num = numincons
            read(1,*)moltype,atom1,atom2,atom3,atom4,bondlength,bondlength2,bondlength3
            frac_shake(i)%atom1 = atom1
            frac_shake(i)%atom2 = atom2
            frac_shake(i)%atom3 = atom3
            frac_shake(i)%atom4 = atom4
            frac_shake(i)%bond = bondlength
            frac_shake(i)%bond2 = bondlength2
            frac_shake(i)%bond3 = bondlength3

            ncons = ncons + nummolarray(moltype)*3
            !if(i.le.nummolshake(1))then
            !   ncons = ncons + nummolarray(1)*3
            !else
            !   ncons = ncons + nummolarray(2)*3
            !endif

         elseif(numincons.eq.-3)then

            frac_shake(i)%num = numincons
            read(1,*)moltype,atom1,atom2,atom3,bondlength,bondlength2,bondlength3
            frac_shake(i)%atom1 = atom1
            frac_shake(i)%atom2 = atom2
            frac_shake(i)%atom3 = atom3
            frac_shake(i)%bond  = bondlength
            frac_shake(i)%bond2 = bondlength2
            frac_shake(i)%bond3 = bondlength3


           ncons = ncons + nummolarray(moltype)*3
           !if(i.le.nummolshake(1))then
           !   ncons = ncons + nummolarray(1)*3
           !else
           !   ncons = ncons + nummolarray(2)*3
           !endif
           
         endif
      enddo
      
  

    end subroutine read_topol




    subroutine check_gro_file
      implicit none
      double precision :: dx,dy,dz
      integer :: i,j
      
      print*,'done with reading in gro file'
      
      print*,'LETS CHECK THE GRO FILE'
      
      do i =1,nummoltype
         do j = 1,numMolAtom(i)
         enddo
      enddo


    end subroutine check_gro_file

    subroutine check_topol_file
      implicit none
      double precision :: local_eps,local_sig
      integer :: loop1,loop2
      integer :: total,loop,offset
      
      open(unit = 501, file = 'check_topol.dat')
      PRINT*,'FORCE FIELD PARAMETERS FOR SIMULATION,sig,eps,q'
      do loop = 1 , numAtomType
         do loop2 = 1,numAtomType
            
            offset = numAtomType*(loop-1)+loop2
            
            local_eps = sqrt(eps(loop)*eps(loop2))

            local_sig = (sig(loop)+sig(loop2))
            lj1(offset) = 12.0d0*local_eps*(local_sig)**12
            lj2(offset) = 12.0d0*local_eps*(local_sig)**6
            lj3(offset) = local_eps*(local_sig)**12
            lj4(offset) = 2.0d0*local_eps*(local_sig)**6


            !lj1(offset) = 48.0d0*local_eps*(local_sig)**12
            !lj2(offset) = 24.0d0*local_eps*(local_sig)**6
            !lj3(offset) = 4.0d0*local_eps*(local_sig)**12
            !lj4(offset) = 4.0d0*local_eps*(local_sig)**6

            !--- 1-4 interaction components
            lj14_1(loop,loop2)   = 0.50d0*12.0d0*local_eps*(local_sig)**12
            lj14_2(loop,loop2)   = 0.50d0*12.0d0*local_eps*(local_sig)**6
            lj14_3(loop,loop2)   = 0.50d0*local_eps*(local_sig)**12
            lj14_4(loop,loop2)   = local_eps*(local_sig)**6


            !lj14_1(loop,loop2)   = 0.50d0*48.0d0*local_eps*(local_sig)**12
            !lj14_2(loop,loop2)   = 0.50d0*24.0d0*local_eps*(local_sig)**6
            !lj14_3(loop,loop2) = 0.50d0*4.0d0*local_eps*(local_sig)**12
            !lj14_4(loop,loop2) = 0.50d0*4.0d0*local_eps*(local_sig)**6

            write(501,*)'pair',loop,loop2
            write(501,*)lj1(offset),lj2(offset),lj3(offset),lj4(offset)
            write(501,*)lj14_1(loop,loop2),lj14_2(loop,loop2),&
                 lj14_3(loop,loop2),lj14_4(loop,loop2)
            write(501,*)

            
         enddo
      enddo
           
  
      close(501)
    end subroutine check_topol_file
 
    

  end module mod_readfile
