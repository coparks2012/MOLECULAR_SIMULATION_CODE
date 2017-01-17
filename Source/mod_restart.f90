module mod_restart
  use get_started
  use mod_memory
  use mod_readfile
  use mod_init_posit
  use mod_force
  use mod_file
  use mod_LRC
  use mod_neighbor
  use mod_pbc
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
  use mod_LRC
  use mod_makefolders
  use mod_langevin
  use ifport
  implicit none

  contains

    subroutine restart
      implicit none
      interface
         subroutine init_pme(pme_cutoff, lx, ly, lz, nx, ny, nz, n) bind (c)
           use iso_c_binding
           integer ( c_int ), VALUE :: n,nx,ny,nz
           real ( c_float ), VALUE :: pme_cutoff,lx,ly,lz
         end subroutine init_pme
      end interface
      double precision :: sum
      integer :: id,jflag
      integer :: i,j,k,num,greatest
      integer :: T1,T2,clock_rate,clock_max
      character(len=50) :: exec,file,newton,auto,rms,filename
      !---calculate interaction parameters
      
      print*,'WELCOME TO RESTART AND RANK',rank


      !---make folders
      call make_MPI_files

      print*,'hello from restart'
      call system_clock(T1,clock_rate,clock_max)
      call read_files_wrapper
      
      print*,'NUMBER OF MOLECULES',numMolArray(1)
      
      
      filename = 'restart_files/beta_1.85nm_solv_em.gro'
      filename = 'restart_files/beta_1.5nm_solv.gro'
      box = 97.99706

      timeclock = 0.0d0
      if(rank.eq.0)then
         temp = 275
      elseif(rank.eq.1)then
         temp = 280
      elseif(rank.eq.2)then
         temp = 285
      elseif(rank.eq.3)then
         temp = 290
      elseif(rank.eq.4)then
         temp = 295
      elseif(rank.eq.5)then
         temp = 300
      endif

      call init_variables_wrapper
      !ngrid_x,y,z should be multipler of 2,3 or 5
      ngrid_x = 100
      ngrid_y = 100
      ngrid_z = 100


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
      

      !---assign types,charges
      call assign_atom_types
      call assign_charge_types

      !---build mol list and bond list
      call build_mol_list
      call build_bond_list_wrapper

      !---read in positions
      call read_in_grofile(filename)
      !call read_in_lammps(filename)
      !call read_in_lammps_beta(filename)
      !call read_in_lammps_velocity(filename)
      !call center_coords


      !---build initial PTTA. MUST BE BEFORE PBC
      !call build_ptta

      
      call pbc(0)
      !call pbc_triclinic

      !---iitialize cluster neighbor lists
      !call prune_template
      !call center_array_template(1)
      !call Grofile_template
      !call ptta_gro('ptta-start.gro')
       call mol_GroFile('restart_one.gro')
      
      !---initialize global id
      do i = 1,np
         global_id(i)=i
      enddo

      num = get_greatestCluster(0,1,greatest)
      print*,'what the hell is num',num
      call make_grofile_cluster_step(0,greatest,num)

          
      
      
  
      !---assign initial random velocities
 
      call init_velocity


      !---initialize cluster neighbor lists
       call center_array_template(1)
       call Grofile_template
      
     
      !---initialize cell neighbors
      if(vlist_flag.ne.0.and.vlist_flag.ne.4)then
         if(newton3rd_flag.eq.1)then
            call neighcell_newton
         elseif(newton3rd_flag.eq.0)then
            call neighcell_nonewton
         endif
         call neighbor_wrapper(0)
      endif
      
    
      !--- initialize random number generator
      call init_rand_stream
      call calculate_ke
      call long_range_corr_init
      call calc_long_range_correction

      time_force=0.0d0
      call force_driver(0,vlist_flag)
      if(nshake.gt.0)then
         shake_time = 0.0d0
        call shake(0)
      endif
      
      call calculate_pressure

      

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
      print*,'temperature set point',temp

      if(coul_flag.eq.2)then
         tot_e = potential + e_bond + e_angle + e_dihedral + e_improper + pot_14 + &
              e_coul + e_coul_14 +e_coul_long + vdw_corr+ke
      elseif(coul_flag.eq.1)then
         tot_e = potential + e_bond + e_angle + e_dihedral + e_improper + pot_14 + &
              e_coul + e_coul_14 +e_self + vdw_corr+ke
      endif
    
      print*,'initial total PE',tot_e-ke
      print*,'intitial total energy',tot_e
      call system_clock(T2,clock_rate,clock_max)
      time_initialize = real(T2-T1)/real(clock_rate)
      print*,'this is our time initializing',time_initialize
      print*
      print*
      print*
      
     

      num = get_greatestCluster(0,1,greatest)
      print*,'what the hell is num',num
      !call dump_global_velocity_id_mpi_final(0)
      stop
    end subroutine restart

  end module mod_restart
