module mod_restart_datadump
  use get_started
  use mod_memory
  use mod_readfile
  use mod_init_posit
  use mod_force
  use mod_file
  use mod_LRC
  use mod_neighbor
  use mod_pbc
  use mod_prune
  use mod_GroFile
  use mod_buildbond_lists
  use mod_assign_types
  use mod_nondimensionalize
  use mod_sample
  use mod_cluster
  use mod_neighbor_cluster
  use mod_build_ptta
  use mod_datadump
  use mod_shake
  use mod_ensemble
  use mod_langevin
  use mod_makefolders
  use ifport
  implicit none

  contains

    subroutine restart_datadump
      implicit none
      interface
         subroutine init_pme(pme_cutoff, lx, ly, lz, nx, ny, nz, n) bind (c)
           use iso_c_binding
           integer ( c_int ), VALUE :: n,nx,ny,nz
           real ( c_float ), VALUE :: pme_cutoff,lx,ly,lz
         end subroutine init_pme
      end interface
      double precision :: sum,calctemp
      integer :: i,j,k,num,greatest
      integer :: T1,T2,clock_rate,clock_max,status
      character(len=300) :: exec,file,newton,auto,rms,filename,intchar
      
      !---calculate interaction parameters
      print*,'HELLO FROM MPI AMBER RUNS FINAL DIFF FLUSH'
      !-----------------------------------
      ! determine path variable
      status = getcwd( path_string )
      path_string = adjustl(path_string)
      path_string = trim(path_string)
      path_string = trim(path_string)//'/'
      if(status.ne.0)then
         print*,'ERROR DETERMINING STATUS'
      endif
      print*,'what is current path',trim(path_string)
      

      !===================================================================
      ! intiialize restart mult
      restart_mult = 0

      !=====================================================================
      ! read in restart mult
      open(unit = 99999,file = 'restart_file.restart')
      read(99999,*)restart_mult
      close(99999)

      print*,'we are now in restart datadump'
      call make_MPI_files
      call system_clock(T1,clock_rate,clock_max)
      call read_files_wrapper
   


      timeclock = 0.0d0
 
      !=====================================================================
      ! read in restart mult
      open(unit = 99999,file = 'restart_file.restart')
      read(99999,*)restart_mult

      write(filename,'(i8)')restart_mult*mdeq
      
      filename = adjustl(filename)
      filename = trim(filename)
      filename = 'restart_'//trim(filename)
      

      if(restart_mult.eq.0.and.file_type_flag.eq.1)then

         print*,'reading in grofile box'
         print*,'filename',trim(filename_gro)
         call read_in_grofile_box( trim(filename_gro) )
         print*,'read in box length from gro file',box
         timeclock = 0.0d0

      elseif(restart_mult.eq.0.and.file_type_flag.eq.2)then
         print*,'READING IN RESTART FILE FILE TYPE FLAG IS 2'
  

         open(unit =1000,file=trim(filename_restart))
         read(1000,*)box
         close(1000)
         print*,'what is box length',box



      elseif(restart_mult.eq.0.and.file_type_flag.eq.0)then
         print*,'READING IN RESTART FILE FILE TYPE FLAG IS 2'
         box = 77.6077
         print*,'what is box length',box

      else
         !======================================================================
         ! FILENAME SPECIFICATION FOR FILE RESTART PRODUCTION RUNS

         print*,'READING IN RESTART FILE'
         print*,'NUMBER OF PARTICLES:;',np,nummolarray(1)
         write(intchar,'(i8)')rank
         intchar = adjustl(intchar)
         filename = 'MPIrestart_final'//'/'//trim(intchar)//'/'//trim(filename)

         open(unit =1000,file=trim(filename))
         read(1000,*)box
         close(1000)
         print*,'what is box length',box
         
         !=======================================================================
      endif


      print*,'heres our filename',trim(filename)
      
      !---calculate number of particles/bonds/angles/dihedrals


      

      call init_variables_wrapper

      print*,'box to grid ratio',box/ngrid_x
      if(box/ngrid_x.gt.1.0)then
         print*,'box to grid ratio is too big'
         stop
      endif

      call init_pme(cut_coul, box, box, box, ngrid_x, ngrid_y, ngrid_z, np)
     
  
      
      !---initialize mncellT and mncellD
      call init_list_cluster
      call init_list
      
      !---allocate arrays
      call memory
      
      !---calculate neighboring cells for clustering vlist
      call mol_neigh
      
      
      !---print values to screen/file
      call print_input
      
      call assign_atom_types
      call assign_charge_types
      call build_mol_list
      call build_bond_list_wrapper
      
      print*,'what is restart_mult',restart_mult
      if(restart_mult.eq.0.and.file_type_flag.eq.1)then

         print*,'=================================='
         print*,'reading initil grofile coordinates'
         print*
         call read_in_grofile(trim(filename_gro))
         call init_velocity
         timeclock = 0.0d0
      elseif(restart_mult.eq.0.and.file_type_flag.eq.2)then
         print*,'=================================='
         print*,'reading in filename restart'
         print*
         call init_velocity
         call read_dump_global_id_v(trim(filename_restart))
         timeclock = 0.0d0
      elseif(restart_mult.eq.0.and.file_type_flag.eq.0)then
         print*,'=================================='
         print*,'reading in filename restart'
         print*
         call init_velocity
         call read_in_lammps(trim(filename_lammps))
         timeclock = 0.0d0
      else
         ncons = ncons + 3
         call read_dump_global_id_v(trim(filename))
      endif
    
      call pbc(0)

      
      call mol_GroFile('restart_datadump.gro')
      

      !---initialize global id
      do i = 1,np
         global_id(i)=i
      enddo
         
      num = get_greatestCluster(0,1,greatest)
      
    
      call make_grofile_cluster_step(0,greatest,num)
          
      !---initialize cell neighbors
      if(vlist_flag.ne.0.and.vlist_flag.ne.4)then
         if(newton3rd_flag.eq.1)then
            call neighcell_newton
         elseif(newton3rd_flag.eq.0)then
            call neighcell_nonewton
         endif


         call neighbor_wrapper(0)
      endif
      
      !---initialize random number generator
      call init_rand_stream
      call long_range_corr_init
      call calc_long_range_correction
      call force_driver(0,vlist_flag)
      if(nshake.gt.0)then
         shake_time = 0.0d0
         call shake(0)
      endif
      call calculate_ke
      call calculate_pressure

      write(1234,*)'potential',potential
      write(1234,*)'14 term',pot_14
      write(1234,*)'total vdw',(potential+pot_14)
      write(1234,*)'electric 14:',e_coul_14
      write(1234,*)'electric non bond',e_coul
      write(1234,*)'e_nonbond+e_coul_14',e_coul+e_coul_14
      write(1234,*)'self energy include in long',e_self
      write(1234,*)'long range energy',e_coul_long
      if(coul_flag.eq.2)then
         write(1234,*)'total electric',(e_coul+e_coul_14+e_coul_long)
      else
         write(1234,*)'total electric',(e_coul+e_coul_14+e_self)
      endif
      write(1234,*)'e_bond',e_bond
      write(1234,*)'e_angle',e_angle
      write(1234,*)'e_dihedral',e_dihedral
      write(1234,*)'e_imp',e_improper
      write(1234,*)'virialx',virialx
      write(1234,*)'virialy',virialy
      write(1234,*)'virialz',virialz
      write(1234,*)'KE',ke
      write(1234,*)'long range vdw',vdw_corr
      write(1234,*)'long range pres',pres_corr
      write(1234,*)'pressure',pressure
      tot_e = potential + pot_14 + e_coul + e_coul_14 + &
           e_bond + e_angle + e_dihedral + e_improper + ke
       

      print*,'=============================================================='
      print*,'PRINTING ENERGY INFORMATION FOR RANK',rank

      print*,'---------------------------------------------------'
      print*,'Welcome to Sofia: a molecular software code for nucleation'
      print*,' -------------------------------------------------'
      
      print*,'monte carlo ensemble:',ensemble
      print*,' Lennard-Jones parameters:'
      print*,'   np=',np
      print*,'   rcut= ',rcut
      print*,'   rcut2=',rcut2
      print*,'   rcut= ',cut_coul
      print*,'   rcut2=',cut_coulsq
      print*,'   neighbor par',rvrc
      
      print*,'   total number cons',ncons
      
      
      
      print*,'   mrcut= ',mrcut
      print*,'   mrcut2=',mrcut2
      print*,'   rv:',rv
      print*,'   rvrc2',rvrc2
      print*,'   Simulation box length: ',box
      print*,'   Volume               : ',vol
      print*,'   box                   :' ,box
      print*,'   hbox                  :' ,hbox
      print*,'   Temperature          : ',temp      
      print*,'   # MD equil steps     : ',mdeq
      print*,'   Sampling freq  : ',nsample
     
      print*
      print*
      print*,'restart filename',trim(filename)
      print*,'WHAT IS CURRENT BOX SIZE',box
      print*,'what is cluster size num',num
      print*,'what is timeclock',timeclock
      print*,'this is initial KE',sum*amu2e
      print*,'kboltz',kboltz
      print*,'ncons',ncons
      print*,'this is set point temperature',temp
      print*,'potential',potential
      print*,'14 term',pot_14
      print*,'total vdw',(potential+pot_14)
      print*,'with corr',(potential+pot_14+vdw_corr)
      print*,'electric 14:',e_coul_14
      print*,'electric non bond',e_coul
      print*,'e_nonbond+e_coul_14',e_coul+e_coul_14
      print*,'self energy include in long',e_self
      print*,'long range energy',e_coul_long
      if(coul_flag.eq.2)then
         print*,'total electric',(e_coul+e_coul_14+e_coul_long)
      else
         print*,'total electric',(e_coul+e_coul_14+e_self)
      endif
      print*,'e_bond',e_bond
      print*,'e_angle',e_angle
      print*,'e_dihedral',e_dihedral
      print*,'e_imp',e_improper
      print*,'virialx',virialx
      print*,'virialy',virialy
      print*,'virialz',virialz
      print*,'KE',ke
      print*,'long range vdw',vdw_corr
      print*,'long range pres',pres_corr
      print*,'pressure',pressure

      print*,'=========================================================='
      print*
      call system_clock(T2,clock_rate,clock_max)
      time_initialize = real(T2-T1)/real(clock_rate)
      

      tot_e = potential + pot_14 + e_coul + e_coul_14 + &
           e_bond + e_angle + e_dihedral + e_improper + ke
       




      !---open sample files
      call open_sample_files
     


      call system_clock(T2,clock_rate,clock_max)
      time_initialize = real(T2-T1)/real(clock_rate)

      print*,'what is our time to initialize',time_initialize
     
      
    end subroutine restart_datadump







    
    subroutine restart_production
      implicit none
      interface
         subroutine init_pme(pme_cutoff, lx, ly, lz, nx, ny, nz, n) bind (c)
           use iso_c_binding
           integer ( c_int ), VALUE :: n,nx,ny,nz
           real ( c_float ), VALUE :: pme_cutoff,lx,ly,lz
         end subroutine init_pme
      end interface
      logical :: e,result
      double precision :: sum,calctemp
      integer :: i,k,num,greatest,ngrid_x,ngrid_y,ngrid_z
      integer :: j,jflag
      integer :: T1,T2,clock_rate,clock_max
      character(len=50) :: exec,file,newton,auto,rms,filename
      !---calculate interaction parameters


      print*,'we are restarting production'
      if(rank.eq.0)then
         temp = 270
         filename = '270_CLUS_AMBER.restart'
      endif
      if(rank.eq.1)then
         temp = 300
         filename = '300_CLUS_AMBER.restart'
      endif
      if(rank.eq.2)then
         temp = 330
         filename = '300_CLUS_AMBER.restart'
      endif
      if(rank.eq.3)then
         filename = '360_CLUS_AMBER.restart'
         temp = 360
      endif
      print*,'what is filename',filename


      open(unit =1000,file=filename)
      read(1000,*)box
      close(1000)
      print*,'what is box length',box     
      vol = box*box*box
      hbox = 0.50*box
      ibox = 1.0/box

    


      call assign_atom_types
      call assign_charge_types
      call build_mol_list
      call build_bond_list_wrapper
      call read_dump_global_id_v(filename)
      !call read_in_grofile(filename)
      call build_ptta_wrap
      call pbc(0)
      call init_velocity
     
      
      !---initialize cell neighbors
      if(vlist_flag.ne.0.and.vlist_flag.ne.4)then
         if(newton3rd_flag.eq.1)then
            call neighcell_newton
         elseif(newton3rd_flag.eq.0)then
            call neighcell_nonewton
         endif
         !---initialize global id
         do i = 1,np
            global_id(i)=i
         enddo
         call neighbor_wrapper(0)
      endif
      
     
      
      call force_driver(0,vlist_flag)
      call long_range_corr_init
      call calc_long_range_correction
      if(nshake.gt.0)then
         shake_time = 0.0d0
       !  call shake(0)
      endif
       call calculate_ke
      call calculate_pressure

      write(1234,*)'restarting simulation'
      write(1234,*)'potential',potential
      write(1234,*)'14 term',pot_14
      write(1234,*)'total vdw',(potential+pot_14)
      write(1234,*)'electric 14:',e_coul_14
      write(1234,*)'electric non bond',e_coul
      write(1234,*)'e_nonbond+e_coul_14',e_coul+e_coul_14
      write(1234,*)'self energy include in long',e_self
      write(1234,*)'long range energy',e_coul_long
      if(coul_flag.eq.2)then
         write(1234,*)'total electric',(e_coul+e_coul_14+e_coul_long)
      else
         write(1234,*)'total electric',(e_coul+e_coul_14+e_self)
      endif
      write(1234,*)'e_bond',e_bond
      write(1234,*)'e_angle',e_angle
      write(1234,*)'e_dihedral',e_dihedral
      write(1234,*)'e_imp',e_improper
      write(1234,*)'virialx',virialx
      write(1234,*)'virialy',virialy
      write(1234,*)'virialz',virialz
      write(1234,*)'KE',ke
      write(1234,*)'long range vdw',vdw_corr
      write(1234,*)'long range pres',pres_corr
      write(1234,*)'pressure',pressure

      print*,'potential',potential
      print*,'14 term',pot_14
      print*,'total vdw',(potential+pot_14)
      print*,'with corr',(potential+pot_14+vdw_corr)
      print*,'electric 14:',e_coul_14
      print*,'electric non bond',e_coul
      print*,'e_nonbond+e_coul_14',e_coul+e_coul_14
      print*,'self energy include in long',e_self
      print*,'long range energy',e_coul_long
      if(coul_flag.eq.2)then
         print*,'total electric',(e_coul+e_coul_14+e_coul_long)
      else
         print*,'total electric',(e_coul+e_coul_14+e_self)
      endif
      print*,'e_bond',e_bond
      print*,'e_angle',e_angle
      print*,'e_dihedral',e_dihedral
      print*,'e_imp',e_improper
      print*,'virialx',virialx
      print*,'virialy',virialy
      print*,'virialz',virialz
      print*,'KE',ke
      print*,'long range vdw',vdw_corr
      print*,'long range pres',pres_corr
      print*,'pressure',pressure
    




      tot_e = potential + pot_14 + e_coul + e_coul_14 + &
           e_bond + e_angle + e_dihedral + e_improper + ke
       
      

      call system_clock(T2,clock_rate,clock_max)
      time_initialize = real(T2-T1)/real(clock_rate)
      

      tot_e = potential + pot_14 + e_coul + e_coul_14 + &
           e_bond + e_angle + e_dihedral + e_improper + ke
       

      call system_clock(T2,clock_rate,clock_max)
      time_initialize = real(T2-T1)/real(clock_rate)

      print*,'what is our time to initialize',time_initialize
      call sample
      
    end subroutine restart_production










  end module mod_restart_datadump
