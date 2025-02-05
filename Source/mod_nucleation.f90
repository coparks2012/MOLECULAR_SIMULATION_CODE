module mod_nucleation

  use global
  use mod_neighbor
  use mod_cluster
  use mod_force
  use mod_pbc
  use mod_ensemble
  use mod_shake
  use mod_force
  use mod_sample
  use mod_adjust
  use mod_datadump

  use mod_assign_types
  use mod_buildbond_lists
  
  !use mod_file_manager
  !use mpi
  !use mod_MPI
  implicit none
contains

  
  !---------------------------------------------------------------------
  !*** obtains cluster of size ngoal
  !*** 
  !***
  subroutine grow_cluster(ngoal,criterium)
    integer :: ngoal,criterium
    integer :: n0
    integer :: T1,T2,clock_rate,clock_max
    integer :: w0,w1,greatest

    call system_clock(T1,clock_rate,clock_max)

    naccept = 0
    nreject = 0
    cc      = 0
    print*
    print*
    print*,'hello from gro_cluster!'
    !call print_cluster_size(693,w1,greatest,0,0,0.0d0)
    call print_cluster_size(w1,greatest)
    print*,ngoal
    print*,'initial cluster size:',w1
    print*

    
    call obtain_cluster_size(ngoal,criterium)

    print*,'final cluster size from gro cluster:'
    !call print_cluster_size(693,w0,greatest,0,0,0.0d0)
    call print_cluster_size(w1,greatest)
    call make_grofile_cluster(greatest,w1,'obtain_clus.gro')
    !call WriteDataToFile('Gro_Data/final_growcluster.dat')


    call close_GroFiles



    call system_clock(T2,clock_rate,clock_max)
    print*,'ELAPSED TIME IN GRO CLUSTER:',real(T2-T1)/real(clock_rate)

  end subroutine grow_cluster




  !---------------------------------------------------------------------
  !*** obtains a cluster of size n0
  !***
  subroutine obtain_cluster_size(n0,criterium)
    implicit none
    double precision :: pref
    integer :: n0,w1,w0,criterium,greatest
    integer :: count
    character(len=40)  :: datadump
  
   
 

    
    !call print_cluster_size(692,w0,greatest,0,0,0.0d0)
    call print_cluster_size(w0,greatest)


    print*,'hello from obtain cluster size!:'
    print*,'initial size,goal,criterium:',w0,n0,criterium,kn,numStep,numSweep
    print*,'hello from obtain cluster size!:'
    print*,'initial size,goal,criterium:',w0,n0,criterium
    print*,'number of MC moves per bias move:', numStep
    print*,'spring constant:',kn
    print*
    print*
    
    count = 0
    do while(w0.ne.n0)

       call nuc_move(w0,n0,criterium)
       count = count + 1

       if(mod(count,20).eq.0)then
       !  call adjust
       endif

       if(mod(count,100).eq.0)then
          !print*,w0,rank
          print*,count,w0
       endif
    enddo

    !---print simulation information to file
    !if(n0.ne.rank)then
    !   write(datadump,'(i8)')n0
    !   datadump = ADJUSTL(datadump)
    !   datadump = trim(datadump)//'.dat'
    !   call WriteDataToFile(datadump)
    !endif
    
    print*
    !print*,'FINISHING OBTAIN CLUSTER SIZE:rank,w0',rank,w0
    print*,'FINISHING OBTAIN CLUSTER SIZE:',w0
    print*
    
 
    
  end subroutine obtain_cluster_size


  !---------------------------------------------------------------------
  !*** perform umbrella sampling. Parallel tempering called here
  !***
  !***
  !subroutine equil_umbrella_sample(n0,criterium,potential)
  !  implicit none
  !  double precision :: potential
  !  double precision :: E0,E1
  !  double precision :: delta
  !  double precision :: rand1
  !  integer :: criterium
  !  integer :: n0
  !  integer :: disp
  !  integer :: loc
  !  integer :: i,loop
  !  integer :: w0,greatest
  !  integer :: msgtag
  !  integer :: n0temp,n1temp
  !  integer :: p0temp,p1temp
  !  integer :: w0temp,w1temp
  !  integer :: ierr
  !  integer :: T1,T2,clock_rate,clock_max



  !  print*,'WHAT UP FROM RANK',rank,'WERE IN equil UMBRELLA SAMPLE'
  !  call mpi_barrier(mpi_comm_world,ierr)
    
  !  write(695,*)'BEGINNING EQUIL umbrella SAMPLE:',rank
  !  write(695,*)

    
  !  call system_clock(T1,clock_rate,clock_max)

  !  call aos2xyz()
  !  w0 = get_greatestCluster(criterium,greatest,0)


  !  !---change numSweep and numStep
  !  call system_clock(T1,clock_rate,clock_max)


  !  write(695,*)'numSweepEquil,numStep,kn,n0:',numSweepEquil,numStep,kn,n0

  !  print*,'numSweepEquil,numStep,kn,n0:',numSweepEquil,numStep,kn,n0

  !  nuchisto(:) = 0
  !  do i=1,numSweepEquil
       
  !     call nuc_move(w0,n0,criterium,potential)
  !     
       
  !     write(689,*)potential+ecorr
  !     write(700,*)potential+ecorr
  !     write(690,*)vol
       
       
  !     if(mod(i,100).eq.0)then
  !        call parallel_temper(w0,n0)
  !     endif
       
       
  !  enddo
    
  !  call system_clock(T2,clock_rate,clock_max)
  !  write(695,*),'elapsed time umbrella sample:',real(T2-T1)/real(clock_rate),rank
    
    

  !  call dump_umbrella(0,2)

  !end subroutine equil_umbrella_sample







  !---------------------------------------------------------------------
  !*** perform umbrella sampling. Parallel tempering called here
  !***
  !***
  !subroutine umbrella_sample(n0,criterium,potential)
  !  implicit none
  !  double precision :: potential
  !  double precision :: E0,E1
  !  double precision :: delta
  !  double precision :: rand1
  !  integer :: criterium
  !  integer :: n0
  !  integer :: disp
  !  integer :: loc
  !  integer :: i,loop
  !  integer :: w0,greatest
  !  integer :: msgtag
  !  integer :: n0temp,n1temp
  !  integer :: p0temp,p1temp
  !  integer :: w0temp,w1temp
  !  integer :: ierr
  !  integer :: T1,T2,clock_rate,clock_max



  !  print*,'WHAT UP FROM RANK',rank,'WERE IN UMBRELLA SAMPLE'
  !  call mpi_barrier(mpi_comm_world,ierr)
    
  !  write(695,*)'BEGINNING UMBRELLA SAMPLE:',rank
  !  write(695,*)

    
  !  call system_clock(T1,clock_rate,clock_max)

  !  call aos2xyz()
  !  w0 = get_greatestCluster(criterium,greatest,0)

  !  nuchisto(:) = 0

  !  !---change numSweep and numStep

  !  call system_clock(T1,clock_rate,clock_max)


  !   write(695,*),numSweep,numStep,kn,n0

  !   do i=1,numSweep

  !     call nuc_move(w0,n0,criterium,potential)
       
  !     if(w0.le.nucMax)then
  !        nuchisto(w0) = nuchisto(w0) + 1
  !     else
  !        print*,'error: nucMax lt w0: nucMax,w0,rank:',nucMax,w0,rank
  !     endif


  !     write(689,*)potential
  !     write(690,*)vol
       
       
  !     if(mod(i,100).eq.0)then
  !        call parallel_temper(w0,n0)
  !     endif
       
  !     if(mod(i,5000).eq.0)then
  !        call dump_umbrella(i,0)
  !     endif
       
  !  enddo
    
  !  call system_clock(T2,clock_rate,clock_max)
  !  write(695,*),'elapsed time umbrella sample:',real(T2-T1)/real(clock_rate),rank
    
    

  !  call dump_umbrella(0,1)
  !  call create_time_series_file(n0)
  !  call close_US

  !end subroutine umbrella_sample




  !---------------------------------------------------------------------
  !*** perform monte carlo. Accept/reject based off change in biasing potential
  !***
  !***
  subroutine nuc_move(w0,n0,criterium)
    implicit none
    double precision :: ran1
    double precision :: deltaE,delta0,E0,E1
    double precision :: sum,calctemp
    integer :: nflag
    integer :: criterium
    integer :: ar,num
    integer :: ar_vol
    integer :: T1, T2, clock_rate, clock_max
    integer :: j,i,k,rep
    integer :: w0,w1,n0,greatest,ierr,flag
    integer :: ini,end
    character(len=50) :: file1,file2,rank_char



    !--- store initial coordinates
    do i =1,np
       p0(i)%x=position(i)%x; p0(i)%y = position(i)%y; p0(i)%z=position(i)%z
       p0(i)%type = position(i)%type
    enddo

    do i = 1,np
       v0(i)%x = v(i)%x; v0(i)%y = v(i)%y; v0(i)%z = v(i)%z
    enddo

    do i = 1,nshake
       num = shakeatom(i)%num
       shakeatom0(i)%num = num
       
       if(num.eq.1)then
          shakeatom0(i)%bond  = shakeatom(i)%bond
          shakeatom0(i)%atom1 = shakeatom(i)%atom1
          shakeatom0(i)%atom2 = shakeatom(i)%atom2
          
       elseif(num.eq.-3)then
          
          shakeatom0(i)%bond  = shakeatom(i)%bond
          shakeatom0(i)%bond2 = shakeatom(i)%bond2
          shakeatom0(i)%bond3 = shakeatom(i)%bond3
          shakeatom0(i)%atom1 = shakeatom(i)%atom1   
          shakeatom0(i)%atom2 = shakeatom(i)%atom2   
          shakeatom0(i)%atom3 = shakeatom(i)%atom3
          
          
       elseif(num.eq.3)then
          
          shakeatom0(i)%bond  = shakeatom(i)%bond
          shakeatom0(i)%bond2 = shakeatom(i)%bond2
          shakeatom0(i)%atom1 = shakeatom(i)%atom1   
          shakeatom0(i)%atom2 = shakeatom(i)%atom2   
          shakeatom0(i)%atom3 = shakeatom(i)%atom3


          
       elseif(num.eq.4)then
   
          shakeatom0(i)%bond  = shakeatom(i)%bond
          shakeatom0(i)%bond2 = shakeatom(i)%bond2
          shakeatom0(i)%bond3 = shakeatom(i)%bond3
          shakeatom0(i)%atom1 = shakeatom(i)%atom1   
          shakeatom0(i)%atom2 = shakeatom(i)%atom2   
          shakeatom0(i)%atom3 = shakeatom(i)%atom3
          shakeatom0(i)%atom4 = shakeatom(i)%atom4
       endif
    enddo
    
   

    do i = 1,np
       if(mol(1,i).eq.1)then
          ptta0(mol(2,i))%x = ptta(mol(2,i))%x
          ptta0(mol(2,i))%y = ptta(mol(2,i))%y
          ptta0(mol(2,i))%z = ptta(mol(2,i))%z
       endif
    enddo

    
    global_id0(:) = global_id(:)
    mol0(:,:) = mol(:,:)
    q0(:) = q(:)
    bondlist0(:,:)  = bondlist(:,:)
    anglelist0(:,:) = anglelist(:,:)
    dihedlist0(:,:) = dihedlist(:,:)
    improlist0(:,:) = improlist(:,:)
    num1bond0(:)   = num1bond(:)
    num2bond0(:)   = num2bond(:)
    num3bond0(:)   = num3bond(:)
    specbond0(:,:) = specbond(:,:)
    
    ini = w0
    !--- store initial data before sweep    
    E0 = 0.5d0*kn*real((w0-n0)**2)
    write(2098,*)'this is the initial energy:',E0

    do i=1,numStep
       

       cc = cc + 1
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
          flag = 0
          call sample_grocluster(flag)
       endif
       
       if(mod(cc,10).eq.0)then
          call timeseries_ptta_grofile
       endif


    enddo


    
    w1 = get_greatestCluster(i,1,greatest)
    E1 = 0.5d0*kn*real((w1-n0)**2)
    end = w1
    
    
    deltaE = E1-E0
    delta0 = deltaE
    deltaE = exp(-deltaE/(kboltz*temp))
    
    ran1 = rand()

    
    if(deltaE.lt.ran1)then
       nreject = nreject + 1
    
 
       do i =1,np
          position(i)%x=p0(i)%x; position(i)%y = p0(i)%y; position(i)%z=p0(i)%z
          position(i)%type = p0(i)%type
       enddo
       

       do i = 1,np
          v(i)%x = v0(i)%x; v(i)%y = v0(i)%y; v(i)%z = v0(i)%z
       enddo
      
       do i = 1,np
          if(mol0(1,i).eq.1)then
             ptta(mol0(2,i))%x = ptta0(mol0(2,i))%x
             ptta(mol0(2,i))%y = ptta0(mol0(2,i))%y
             ptta(mol0(2,i))%z = ptta0(mol0(2,i))%z
          endif
       enddo


       do i = 1,nshake
          num = shakeatom0(i)%num
          shakeatom(i)%num = num
         
          if(num.eq.1)then
             shakeatom(i)%bond  = shakeatom0(i)%bond
             shakeatom(i)%atom1 = shakeatom0(i)%atom1
             shakeatom(i)%atom2 = shakeatom0(i)%atom2
            
          elseif(num.eq.-3)then
            
             shakeatom(i)%bond  = shakeatom0(i)%bond
             shakeatom(i)%bond2 = shakeatom0(i)%bond2
             shakeatom(i)%bond3 = shakeatom0(i)%bond3
             shakeatom(i)%atom1 = shakeatom0(i)%atom1   
             shakeatom(i)%atom2 = shakeatom0(i)%atom2   
             shakeatom(i)%atom3 = shakeatom0(i)%atom3
             
            
          elseif(num.eq.3)then
            
             shakeatom(i)%bond  = shakeatom0(i)%bond
             shakeatom(i)%bond2 = shakeatom0(i)%bond2
             shakeatom(i)%atom1 = shakeatom0(i)%atom1   
             shakeatom(i)%atom2 = shakeatom0(i)%atom2   
             shakeatom(i)%atom3 = shakeatom0(i)%atom3

         elseif(num.eq.4)then
  
             shakeatom(i)%bond  = shakeatom0(i)%bond
             shakeatom(i)%bond2 = shakeatom0(i)%bond2
             shakeatom(i)%bond3 = shakeatom0(i)%bond3
             shakeatom(i)%atom1 = shakeatom0(i)%atom1   
             shakeatom(i)%atom2 = shakeatom0(i)%atom2   
             shakeatom(i)%atom3 = shakeatom0(i)%atom3
             shakeatom(i)%atom4 = shakeatom0(i)%atom4
          endif
       enddo
     
       mol(:,:) = mol0(:,:)
       q(:) = q0(:)
       bondlist(:,:)  = bondlist0(:,:)
       anglelist(:,:) = anglelist0(:,:)
       dihedlist(:,:) = dihedlist0(:,:)
       improlist(:,:) = improlist0(:,:)
       num1bond(:)    = num1bond0(:)
       num2bond(:)    = num2bond0(:)
       num3bond(:)    = num3bond0(:)
       specbond(:,:)  = specbond0(:,:)  
       global_id(:)   = global_id0(:)
       call thermalize


       call neighbor_wrapper(i)
      
       !---dont do this for production runs!
       call force_driver(0,vlist_flag)
       if(nshake.gt.0)then
          call shake(1)
       endif
       call langevin(dthalf)

       
       write(2098,*)'move rejected:',w1
       write(2098,*)'E0 and E1',E0,E1,deltaE
       write(2098,*)'accept/reject',naccept,nreject
       

       
    else
       naccept = naccept + 1
       w0 = w1
       sum = 0.0d0
       do i = 1, np
          sum = sum + 0.50d0*mass(position(i)%type)*(v(i)%x*v(i)%x+v(i)%y*v(i)%y+v(i)%z*v(i)%z)
       enddo
       calctemp = 2.0d0*sum/(kboltz*(3.0d0*np-ncons))
       
       write(2098,*)'in accepted KE',sum
       write(2098,*)'move accept temp',calctemp
       write(2098,*)'move accepted:',w0
       write(2098,*)'E0 and E1',E0,E1,deltaE
       write(2098,*)'acceot/reject',naccept,nreject
    endif
    write(2098,*)
    write(2098,*)'-----------------------------------------'

    !call ptta_gro('ptta-nucmove.gro')


    print*,'-----------------------------------------'
    print*,kn,naccept,nreject
    print*,ini,end
    print*,E1,E1
    print*,deltaE
    print*,ran1
    print*,'-------------------------------------------'
    print*
  end subroutine nuc_move





    
  !---------------------------------------------------------------------
  !*** open times series file for free energy calculation
  !*** 
  !***
  subroutine create_time_series_file(n0)
    implicit none
    integer, parameter :: rk = selected_real_kind(p=15), ik=4
    integer :: i,n0
    integer :: Oh,j
    integer :: num,disp,loc
    integer :: count,W
 
   ! write(1000,*)'final bias potential of processor:',n0
   ! count = 0
   ! do i = 1,nucMax
   !    num = nuchisto(i)
   !    
   !    do j = 1,num
   !       count = count + 1
   !       write(1000,*)count,i
   !    enddo
   ! enddo

       
       
  end subroutine create_time_series_file





  !---------------------------------------------------------------------
  !*** open files for nucleation simulation
  !*** all files go in folder rank
  !***
  subroutine open_nucleation_files()
    implicit none
    logical :: e
    character(len=150) :: rank_char,rank_char2
    character(len=40)  :: file1,file2,file3,file4
    character(len=40)  :: file5,file6,file7,file8,file9
    character(len=40)  :: file10,file11,file12,file13

    

    file13 = 'total_energy.dat'
    file1 = 'cluster_size.dat'
    file2 = 'obtain_en.dat'
    file3 = 'obtain_vol.dat'
    file4 = 'umbrella_en.dat'
    file5 = 'umbrella_vol.dat'
    file6 = 'time_series.dat'
    file12 ='umbrella_info.dat'
    file7 = 'obtain_Nuc.dat'
    file8 = 'obtain_Info.dat'
    file9 = 'GroInfo.dat'
    file10 = 'Gro_Data'
    file11 = 'US_Data'

    file10 = ADJUSTL(file10)
    file10 = trim(file10)//'/'
    file11 = ADJUSTL(file11)
    file11 = trim(file11)//'/'

  !  write(rank_char,'(i8)')rank
  !  rank_char = ADJUSTL(rank_char)
  !  rank_char = trim(rank_char)//'/'

  !  open(unit = 700, file = trim(rank_char)//trim(file13), status = 'replace')

  !  rank_char = trim(rank_char)//trim(file10)//'/'

  !  INQUIRE(file=trim(rank_char),EXIST=e)
  !  if(.not.e)then
  !     call system('mkdir '//trim(rank_char))
  !  endif
  !  open(unit = 687, file = trim(rank_char)//trim(file2), status='replace')
  !  open(unit = 688, file = trim(rank_char)//trim(file3), status='replace')
  !  open(unit = 691, file = trim(rank_char)//trim(file7), status='replace')
  !  open(unit = 692, file = trim(rank_char)//trim(file8), status='replace')
  !  open(unit = 693, file = trim(rank_char)//trim(file9), status='replace')


  !  write(rank_char2,'(i8)')rank
  !  rank_char2 = ADJUSTL(rank_char2)
  !  rank_char2 = trim(rank_char2)//'/'
  !  rank_char2 = trim(rank_char2)//trim(file11)//'/'

  !  INQUIRE(file=trim(rank_char2),EXIST=e)
  !  if(.not.e)then
  !     call system('mkdir '//trim(rank_char2))
  !  endif
  !  open(unit = 689, file = trim(rank_char2)//trim(file4), status='replace')
  !  open(unit = 690, file = trim(rank_char2)//trim(file5), status='replace')
  !  open(unit = 1000,file = trim(rank_char2)//trim(file6), status='replace')
  !  open(unit = 695 ,file = trim(rank_char2)//trim(file12), status='replace')

  end subroutine open_nucleation_files


  subroutine close_Grofiles
    implicit none
    

    close(687)
    close(688)
    close(691)
    close(692)
    close(693)

  end subroutine close_Grofiles



  subroutine close_US
    implicit none
    close(689)
    close(690)
    close(1000)
    close(695)
  end subroutine close_US


end module mod_nucleation
