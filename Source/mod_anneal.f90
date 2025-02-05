module mod_anneal
  use global
  use mod_file
  use mod_pbc
  use mod_ensemble
  use mod_force
  use mod_neighbor
  use mod_cluster
  use mod_sample
  use mod_shake
  use mod_cluster
  use mod_Gaussian
  use mod_datadump
  contains
    subroutine anneal
      implicit none
      double precision :: rr,calctemp,calctemp2
      integer :: num,greatest
      integer :: i,nflag,count
      integer :: T1,T2,clock_rate,clock_max
      integer :: flag1, flag2, flag3, flag4
      print*
      print*
      print*
      print*,'hello from anneal'
      call system_clock(T1,clock_rate,clock_max)



      PRINT*,'what is anneal_temp',anneal_temp*rtemp,temp*rtemp
      call open_anneal_files


      rr = ranmars(10)
      timeclock = 0.0d0
      time_bond = 0.0d0
      time_nonbond = 0.0d0
      time_angle = 0.0d0
      time_dihed = 0.0d0
      time_impro =0.0d0
      time_data = 0.0d0
      time_neigh = 0.0d0
      time_cache =  0.0d0
      time_clus = 0.0d0
      nn_builds = 0
     
      !call NVT_equil
      print*,'WE ARE NOW DONE WITH THE FIRST EQUIL'
      !call dump_global_velocity_id('global_id_344.xyz')
      !call dump_ptta('dump_ptta_344.xyz')

      open(unit = 5673, file = 'num.dat')

   
      count = 0 

      flag1 = 0; flag2 = 0; flag3 = 0; flag4 = 0

      !---initialize calctemp
      call get_temp(calctemp)

      do while(temp.gt.anneal_temp)

         do while(calctemp.gt.temp)
            do i = 1, 100
               count = count + 1
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
               
               if(mod(i,nsample).eq.0)then
                  call sample_anneal
                  num = get_greatestCluster(0,1,greatest)
                  write(5673,*)num
               endif
            enddo

            call get_temp(calctemp)
         enddo
           
         !--- update thermostat variables
         temp = 0.99d0*temp
         STDtemp = sqrt(temp)
         gconst =  sqrt(12.0*kboltz*temp*fric/dthalf)
         
         if(calctemp*rtemp.le.294.9.and.flag1.eq.0)then
            flag1 = 1
            print*,'count for 294'
            call dump_global_velocity_id('global_id_294test.xyz')
            call dump_ptta('294ptta.xyz')
            
            call get_temp(calctemp2)
            call force_driver(0,vlist_flag)
            
            print*,'-------------------------------------------'
            print*,'294 anneal'
            print*,'calctemp',calctemp2*rtemp
            print*,'potential',potential*reps
            print*,'14 term',pot_14*reps
            print*,'electric brute :',(e_coul_14+e_coul)*reps
            print*,'e_bond',e_bond*reps
            print*,'e_angle',e_angle*reps
            print*,'e_dihedral',e_dihedral*reps
            print*,'e_imp',e_improper*reps
            print*,'-------------------------------------------'
   
            
            !call dump_global_velocity_id('global_id_294anneal.xyz')
            !call dump_ptta('dump_ptta_294anneal.xyz')
            
        
            
         endif
         
         if(calctemp*rtemp.le.274.9.and.flag2.eq.0)then
            flag2 = 1
            print*,'count for 274'
            !call dump_global_velocity_id('global_id_274anneal.xyz')
            !call dump_ptta('dump_ptta_274anneal.xyz')
            
            call dump_global_velocity_id('global_id_274test.xyz')
            call dump_ptta('global_id_274ptta.xyz')


            call get_temp(calctemp2)
            call force_driver(0,vlist_flag)
            print*,'-------------------------------------------'
            print*,'274 anneal'
            print*,'calctemp',calctemp2*rtemp
            print*,'potential',potential*reps
            print*,'14 term',pot_14*reps
            print*,'electric brute :',(e_coul_14+e_coul)*reps
            print*,'e_bond',e_bond*reps
            print*,'e_angle',e_angle*reps
            print*,'e_dihedral',e_dihedral*reps
            print*,'e_imp',e_improper*reps
            print*,'-------------------------------------------'
         endif
         
         if(calctemp*rtemp.le.264.0.and.flag3.eq.0)then
            flag3 = 1
            print*,'count for 260'
            !call dump_global_velocity_id('global_id_260anneal.xyz')
            !call dump_ptta('dump_ptta_260anneal.xyz')

            call dump_global_velocity_id('global_id_260test.xyz')
            call dump_ptta('260ptta.xyz')

            call get_temp(calctemp2)
            call force_driver(0,vlist_flag)
            print*,'-------------------------------------------'
            print*,'264 anneal'
            print*,'calctemp',calctemp2*rtemp
            print*,'potential',potential*reps
            print*,'14 term',pot_14*reps
            print*,'electric brute :',(e_coul_14+e_coul)*reps
            print*,'e_bond',e_bond*reps
            print*,'e_angle',e_angle*reps
            print*,'e_dihedral',e_dihedral*reps
            print*,'e_imp',e_improper*reps
            print*,'-------------------------------------------'

         endif
      end do

      stop


      !call dump_global_velocity_id('global_id_244anneal.xyz')
      !call dump_ptta('dump_ptta_244anneal.xyz')
      call NVT_equil


      call system_clock(T2,clock_rate,clock_max)
      print*,'what the hell is T2 T1 now',T2,T1

      num = get_greatestCluster(0,1,greatest)
      print*,'what the hell is the final num',num

      call make_grofile_cluster_mol(greatest,num*numatns,'cluster.gro')



      print*,'ELAPSED TIME INTEGRATING',real(T2-T1)/real(clock_rate)
      time_integrate = real(T2-T1)/real(clock_rate)
      
      
      print*,'lets check integrate',time_integrate
      print*,'nonbond time',time_nonbond
      print*,'neighbor time',time_neigh
      print*,'time data',time_data
      print*,'heres out total shake time',shake_time
      print*,'% of time spent in shake',100*(shake_time)/time_integrate
      print*,'% of time spent in nonbond',100*(time_nonbond)/time_integrate
      print*,'% of time spent in bond',100*(time_bond)/time_integrate
      print*,'% of time spent in angle',100*(time_angle)/time_integrate
      print*,'% of time spent in dihed',100*(time_dihed)/time_integrate
      print*,'% of time spent in impro',100*(time_impro)/time_integrate
      print*,'$ of time spent in cache',100*(time_cache)/time_integrate
      print*,'% of time spent in data',100*(time_data)/time_integrate
      print*,'% of time spent in neigh',100*(time_neigh)/time_integrate
      print*,'% of time spent in clus ',100*(time_clus)/time_integrate
      print*,'total number of builds',nn_builds
      print*,'average time between builds',real((count+2*mdeq)/nn_builds)
      
      
      
      print*,'current neigh alloc',neigh_alloc
      print*,'current nlistsize',nlistsize
      print*,'current cache size',cache_size
      print*,'potential',potential*reps
      print*,'14 term',pot_14*reps
      print*,'electric brute :',(e_coul_14+e_coul)*reps
      print*,'e_bond',e_bond*reps
      print*,'e_angle',e_angle*reps
      print*,'e_dihedral',e_dihedral*reps
      print*,'e_imp',e_improper*reps
      
      
      
      !call system_clock(T1,clock_rate,clock_max)
      !call dump_global_velocity_id('global_id_244final.xyz')
      !call dump_ptta('dump_ptta_244final.xyz')
      !call system_clock(T2,clock_rate,clock_max)
      !call mol_GroFile('anneal.gro')




    end subroutine anneal


    subroutine NVT_equil
      implicit none
      double precision :: sum
      integer :: i,nflag
      integer :: greatest,num

       timeclock = 0.0d0
       do i= 1,mdeq
          
          !---first half on time integration step
          
          call nve_v(dthalf)
          call nve_x
          call pbc(i)  
          call check_neighbor(nflag,i)
          
          if(vlist_flag.ne.0.and.nflag.eq.1)then
             call neighbor_wrapper(i)
          endif
          
          call force_driver(i,vlist_flag)
          if(nshake.gt.0)then
               call shake(1)
            endif
          call nve_v(dthalf)          
          
          
          timeclock = timeclock + dt
          
          if(mod(i,nsample).eq.0)then
             call sample_anneal

             densmin = 6
             num = get_greatestCluster(0,1,greatest)
             write(78900,*)num
             
             densmin = 7
             num = get_greatestCluster(0,1,greatest)
             write(78901,*)num
             
             densmin = 8
             num = get_greatestCluster(0,1,greatest)
             write(78902,*)num
             
             densmin = 9
             num = get_greatestCluster(0,1,greatest)
             write(78903,*)num

          endif


          if(mod(i,500).eq.0)then
          !   call timeseries_ptta_grofile
          endif

       enddo
       
  end subroutine NVT_equil


  subroutine get_temp(calctemp)
    implicit none
    double precision :: sum,calctemp
    integer :: i
    sum = 0.0d0

      do i = 1, np
         sum = sum + 0.50d0*mass(position(i)%type)*(v(i)%x*v(i)%x+v(i)%y*v(i)%y+v(i)%z*v(i)%z)

      enddo
      calctemp = 2.0d0*sum/(kboltz*(3.0d0*np-ncons))
    end subroutine get_temp

end module mod_anneal
