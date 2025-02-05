module mod_integrate
  
  use global
  use mod_pbc
  use mod_ensemble
  use mod_force
  use mod_sample
  use mod_velScale
  use mod_print_results
  use mod_file
  use mod_neighbor
  use mod_GroFile
  use mod_sample
  use mod_cluster
  use mod_datadump
  use mod_shake
  use mod_Gaussian
  use mod_LRC
  use mod_restart_datadump
  use mpi
  implicit none
   
contains
  
  subroutine integrate
    implicit none
    double precision ::  dx,dy,dz
    double precision ::  sum, calctemp
    integer :: i,nflag,j,loc,w1
    integer:: n0,k
    integer :: T1,T2,clock_rate,clock_max
    integer :: greatest,num,numt
    
    
    !---initialize value of force before loop
    call force_driver(0,vlist_flag)
    if(nshake.gt.0)then
       call shake(1)
    endif
    call sample_eq


   dx = ranmars(10)
   timeclock = 0.0d0
   time_bond = 0.0d0
   time_nonbond = 0.0d0
   time_angle = 0.0d0
   time_dihed = 0.0d0
   time_impro =0.0d0
   time_data = 0.0d0
   time_neigh = 0.0d0
   time_cache =  0.0d0
   time_sort = 0.0d0
   time_lan = 0.0d0
   time_ver = 0.0d0
   time_pme = 0.0d0
   time_force = 0.0d0
   nn_builds = 0
   print*,'running prouction for rank',rank
   call system_clock(T1,clock_rate,clock_max)
   do i= 1,mdeq
      
      !---first half on time integration step
      if(ensemble.eq.'NVE')then
         call nve_v(dthalf)
         call nve_x
         
         call pbc(i)  
         
         if(vlist_flag.ne.0.and.vlist_flag.ne.4)then
            call check_neighbor(nflag,i)
            
            if(nflag.eq.1)then
               call neighbor_wrapper(i)
            endif
            
         endif
         
         call force_driver(i,vlist_flag)
         if(nshake.gt.0)then
            call shake(1)
         endif
         
         call nve_v(dthalf)
         
         !call calc_long_range_correction
         !call calculate_pressure
         !call berendsen
      endif
      
      timeclock = timeclock + dt
      
      if(mod(i,nsample).eq.0)then 
         
         call sample_eq
         call timeseries_ptta_GroFile
      endif
      
      if(mod(i,50000).eq.0)then
         call dump_global_velocity_id_mpi(i,num)
      endif


   enddo
   
   print*,'TOTAL TIME BOND',time_bond
   print*,'TOTAL TIME PME',time_pme
   print*,'TOTAL TIME NONBOND',time_nonbond
   print*,'TOTAL TIME IN FORCE',time_force
   call system_clock(T2,clock_rate,clock_max)
   print*,'ELAPSED TIME INTEGRATING',real(T2-T1)/real(clock_rate)
   time_integrate = real(T2-T1)/real(clock_rate)
   
   call print_equil_results(i)
   

  
 end subroutine integrate

  subroutine integrate_nvt
    implicit none
    double precision ::  dx,dy,dz
    double precision ::  sum, calctemp
    integer :: i,nflag,j,loc,w1
    integer:: n0,k
    integer :: T1,T2,clock_rate,clock_max
    integer :: greatest,num,numt
    
    
    !---initialize value of force before loop
    call force_driver(0,vlist_flag)
    if(nshake.gt.0)then
       call shake(1)
    endif
    call sample_eq
    
    
    dx = ranmars(10)
    timeclock = 0.0d0
    time_bond = 0.0d0
    time_nonbond = 0.0d0
    time_angle = 0.0d0
    time_dihed = 0.0d0
    time_impro =0.0d0
    time_data = 0.0d0
    time_neigh = 0.0d0
    time_cache =  0.0d0
    time_sort = 0.0d0
    time_lan = 0.0d0
    time_ver = 0.0d0
    time_pme = 0.0d0
    time_force = 0.0d0
    nn_builds = 0
    print*,'running production for rank',rank
    call system_clock(T1,clock_rate,clock_max)
    do i= 1,mdeq
       
       !---first half on time integration step
       if(ensemble.eq.'NVE')then
          call nve_v(dthalf)
          call nve_x
          
          call pbc(i)  
          
          if(vlist_flag.ne.0.and.vlist_flag.ne.4)then
             call check_neighbor(nflag,i)
             
             if(nflag.eq.1)then
                call neighbor_wrapper(i)
             endif
             
          endif
          
          call force_driver(i,vlist_flag)
          if(nshake.gt.0)then
             call shake(1)
          endif
          
          call nve_v(dthalf)
          
          
       endif
       
       timeclock = timeclock + dt
       
       if(mod(i,nsample).eq.0)then 
          call calc_long_range_correction
          call calculate_pressure
          call sample_eq
       

          num = get_greatestCluster(0,1,greatest)
          call record_num_eq(num)
       endif

       !if(mod(i,10000).eq.0)then
       !   call timeseries_ptta_grofile
       !endif

       if(mod(i,25000).eq.0)then
          call make_grofile_cluster_step(i,greatest,num)
       endif
       if(mod(i,25000).eq.0)then
          call dump_global_velocity_id_mpi(i,num)
       endif
             
    enddo
    


    print*,'TOTAL TIME BOND',time_bond
    print*,'TOTAL TIME PME',time_pme
    print*,'TOTAL TIME NONBOND',time_nonbond
    print*,'TOTAL TIME IN FORCE',time_force
    call system_clock(T2,clock_rate,clock_max)
    print*,'ELAPSED TIME INTEGRATING',real(T2-T1)/real(clock_rate)
    time_integrate = real(T2-T1)/real(clock_rate)

    call dump_global_velocity_id_mpi_final((restart_mult)*mdeq+mdeq,num)

  end subroutine integrate_nvt
  

  subroutine integrate_npt
    implicit none
    double precision ::  dx,dy,dz
    double precision ::  sum, calctemp
    integer :: i,nflag,j,loc,w1
    integer:: n0,k
    integer :: T1,T2,clock_rate,clock_max
    integer :: greatest,num,numt
    
    print*,'================================================'
    print*,'beginning integration for ensemble NPT'
    
    !---initialize value of force before loop
    call force_driver(0,vlist_flag)
    if(nshake.gt.0)then
       call shake(1)
    endif


    call sample_eq
    num = get_greatestCluster(0,1,greatest)
    call record_num_eq(num)
    
    
    dx = ranmars(10)
    time_clus = 0.0d0
    time_bond = 0.0d0
    time_bond_total = 0.0d0
    time_nonbond = 0.0d0
    time_angle = 0.0d0
    time_dihed = 0.0d0
    time_impro =0.0d0
    time_data = 0.0d0
    time_neigh = 0.0d0
    time_cache =  0.0d0
    time_sort = 0.0d0
    time_lan = 0.0d0
    time_ver = 0.0d0
    time_pme = 0.0d0
    time_force = 0.0d0
    time_12corr = 0.0d0
    time_13corr = 0.0d0
    time_14corr =0.0d0
    nn_builds = 0
    print*,'running production for rank',rank
    call system_clock(T1,clock_rate,clock_max)
    do i= 1,mdeq
       
       !---first half on time integration step
       if(ensemble.eq.'NVE')then
          call nve_v(dthalf)
          call nve_x
          
          call pbc(i)  
          
          if(vlist_flag.ne.0.and.vlist_flag.ne.4)then
             call check_neighbor(nflag,i)
             
             if(nflag.eq.1)then
                call neighbor_wrapper(i)
             endif
             
          endif
          
          call force_driver(i,vlist_flag)
          if(nshake.gt.0)then
             call shake(1)
          endif
          
          call nve_v(dthalf)
          
          call calc_long_range_correction
          call calculate_pressure
          call berendsen
       endif

       
       timeclock = timeclock + dt
       
       if(mod(i,nsample).eq.0)then 
          call sample_eq
  
       endif

    enddo
    call system_clock(T2,clock_rate,clock_max)
    time_integrate = real(T2-T1)/real(clock_rate)
    print*,'ELAPSED TIME INTEGRATING FOR RANK',rank,real(T2-T1)/real(clock_rate)

    print*,'TOTAL TIME bond',time_bond
    print*,'total time bond total',time_bond_total
    print*,'TOTAL TIME angle',time_angle
    print*,'TOTAL TIME dihed',time_dihed
    print*,'TOTAL TIME impro',time_impro
    print*,'TOTAL TIME 12corr',time_12corr
    print*,'TOTAL TIME 13cirr',time_13corr
    print*,'TOTAL TIME 14corr',time_14corr
    print*,'TOTAL TIME PME',time_pme
    print*,'TOTAL TIME NONBOND',time_nonbond
    print*,'TOTAL TIME IN FORCE',time_force
    print*,'TOTAL NEIGH TIME',time_neigh
    print*,'total shake time',shake_time
    print*,'total time lang',time_lan
 
    

    write(1234,*)'TOTAL TIME bond',time_bond
    write(1234,*)'total time bond total',time_bond_total
    write(1234,*)'TOTAL TIME angle',time_angle
    write(1234,*)'TOTAL TIME dihed',time_dihed
    write(1234,*)'TOTAL TIME impro',time_impro
    write(1234,*)'TOTAL TIME 12corr',time_12corr
    write(1234,*)'TOTAL TIME 13cirr',time_13corr
    write(1234,*)'TOTAL TIME 14corr',time_14corr
    write(1234,*)'TOTAL TIME PME',time_pme
    write(1234,*)'TOTAL TIME NONBOND',time_nonbond
    write(1234,*)'TOTAL TIME IN FORCE',time_force
    write(1234,*)'TOTAL NEIGH TIME',time_neigh
    write(1234,*)'total shake time',shake_time
    write(1234,*)'total time lang',time_lan
 

    
    call print_equil_results(i)
    
    
    
  end subroutine integrate_npt

  
  


  subroutine integrate_f
    implicit none
    double precision ::  dx,dy,dz
    double precision ::  sum, calctemp
    integer :: i,nflag,j,loc,w1
    integer:: n0,k,jflag
    integer :: T1,T2,clock_rate,clock_max
    integer :: greatest,num,numt
    
    
    !---initialize vlue of force before loop
    call force_driver(0,vlist_flag)
    call langevin(dthalf)
    if(nshake.gt.0)then
       call shake(0)
    endif


    print*
    print*
    print*,'we are beginning integration with integration flag'

    dx = ranmars(10)
    timeclock = 0.0d0
    time_bond = 0.0d0
    time_nonbond = 0.0d9
    time_angle = 0.0d0
    time_dihed = 0.0d0
    time_impro =0.0d0
    time_data = 0.0d0
    time_neigh = 0.0d0
    time_cache =  0.0d0
    time_sort = 0.0d0
    time_lan = 0.0d0
    time_ver = 0.0d0
    time_single_update = 0.0d0
    time_video = 0.0d0
    nn_builds = 0

    do i= 1,mdeq
       call system_clock(T1,clock_rate,clock_max)

       !---first half on time integration step
       if(ensemble.eq.'NVE')then

          
          call nve_v_integrate_flag(dthalf)
          call nve_x_integrate_flag
          call pbc(i)  


          if(vlist_flag.ne.0.and.vlist_flag.ne.4)then
             call check_neighbor(nflag,i)

             if(nflag.eq.1)then
                call neighbor_wrapper(i)
             endif
           endif
       
          call force_driver(i,vlist_flag)
          call langevin(dthalf)
          if(nshake.gt.0)then
             call shake(1)
          endif
      

     
          
          call nve_v_integrate_flag(dthalf)
          call calc_long_range_correction
          call calculate_pressure
          call berendsen
          
       endif
       
       timeclock = timeclock + dt
 
       if(mod(i,nsample).eq.0)then 
         
          call sample_eq
          num = get_greatestCluster(0,1,greatest)
          call record_num_eq(num)

       endif
       if(mod(i,1000).eq.0)then
          call make_grofile_cluster_step(i,greatest,num)
       endif

       if(mod(i,10000).eq.0)then
          call dump_global_velocity_id_mpi(i,num)
       endif

       if(mod(i,1000).eq.0)then
          call timeseries_ptta_grofile
       endif

       call system_clock(T2,clock_rate,clock_max)
       time_integrate = time_integrate + real(T2-T1)/real(clock_rate)
    enddo

 
   ! call print_equil_results(i)
  end subroutine integrate_f
  
  
 



  subroutine print_equil_results(i)
    implicit none
    double precision :: sum,calctemp
    integer :: i,num,greatest
    integer :: T1,T2,clock_rate,clock_max,IERR
    character(len=1500) :: temp_string,rf
    
    write(1234,*)'----------------------------------------------------------'
    write(1234,*)'----------------------------------------------------------'
    write(1234,*)'----------------------------------------------------------'
    write(1234,*)'FINAL RUN ANALSYSIS FOR RANK',rank

    sum = 0.0d0
    do i = 1, np
       sum = sum + 0.50d0*mass(position(i)%type)*(v(i)%x*v(i)%x+v(i)%y*v(i)%y+v(i)%z*v(i)%z)
    enddo
    calctemp = 2.0d0*sum*amu2e/(kboltz*(3.0d0*np-ncons))
    
    !=========================================================
    num = get_greatestCluster(0,1,greatest)
    write(1234,*)'FINAL CLUSTER SIZE',num
  
    !========================================================
    ! lets write restart files
    ! only rank 0 should do this

    print*,'rank has made it to MPI barrier',rank
    write(1234,*)'rank has made it to barrier',rank
    call mpi_barrier(MPI_COMM_WORLD,ierr)

    write(1234,*)'all ranks have made it to barrier',rank
    print*,'at mpi barrier'
    if(rank.eq.0)then
       open(unit = 999990,status='replace',file = 'restart_file.restart')
       write(999990,*)restart_mult+1
       close(999990)
    endif


    write(1234,*)'lets print data summary for simulation',rank
    write(1234,*)'WHAT IS FINAL BOX LENGTH FOR RESTART',box
    write(1234,*)'this is final temp after equil',calctemp
    write(1234,*)'lets check integrate',time_integrate
    write(1234,*)'nonbond time',time_nonbond
    write(1234,*)'neighbor time',time_neigh
    write(1234,*)'time data',time_data
    write(1234,*)'heres out total shake time',shake_time
    write(1234,*)'total time langevin',time_lan
    
    if(time_integrate.ne.0.0)then
       write(1234,*)'% of time spent in shake',100*(shake_time)/time_integrate
       write(1234,*)'% of time spent in force total',100*(time_force)/time_integrate

       write(1234,*)'% of time spent in nonbond',100*(time_nonbond)/time_integrate
       write(1234,*)'% of tume spent in PME',100*(time_pme)/time_integrate
       write(1234,*)'% of time spent in bond',100*(time_bond)/time_integrate
       write(1234,*)'% of time spent in angle',100*(time_angle)/time_integrate
       write(1234,*)'% of time spent in dihed',100*(time_dihed)/time_integrate
       write(1234,*)'% of time spent in impro',100*(time_impro)/time_integrate
       write(1234,*)'% of time spent in 12corr',100*(time_12corr)/time_integrate
       write(1234,*)'% of time spent in 13corr',100*(time_13corr)/time_integrate
       write(1234,*)'% of time spent in 14corr',100*(time_14corr)/time_integrate
       write(1234,*)'% of time spend in single',100*(time_single_update)/time_integrate
       write(1234,*)'$ of time spent in cache',100*(time_cache)/time_integrate
       write(1234,*)'% of time spent in data',100*(time_data)/time_integrate
       write(1234,*)'% of time spent in neigh',100*(time_neigh)/time_integrate
       write(1234,*)'% of time spent in sort',100*(time_sort)/time_integrate
       write(1234,*)'% of time spent in langevin',100*(time_lan)/time_integrate
       write(1234,*)'% of time spent in VER:',100*(time_ver)/time_integrate
       write(1234,*)'% of time spentin video',100*(time_video)/time_integrate
       write(1234,*)'% of time spentin clustering',100*(time_clus)/time_integrate

    endif
    
    if(nn_builds.ne.0)then
       write(1234,*)'total number of builds',nn_builds
       write(1234,*)'average time between builds',real(mdeq/nn_builds)
    endif
    write(1234,*)'current neigh alloc',neigh_alloc
    write(1234,*)'current nlistsize',nlistsize
    write(1234,*)'current cache size',cache_size
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
    


  end subroutine print_equil_results

    












 
    
   
end module mod_integrate
