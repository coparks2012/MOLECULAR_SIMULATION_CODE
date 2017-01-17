
module mod_sample
  use global
  contains


    subroutine record_num_eq(num)
      implicit none
      integer :: num
      logical :: exist

     
      write(clus_f_id,*)timeclock,num

   
    end subroutine record_num_eq


    
    subroutine sample_eq
      implicit none
      integer :: T1,T2,clock_rate,clock_max
      call system_clock(T1,clock_rate,clock_max)
      
      call get_KE_eq
      call get_PE_eq
      call get_totE_eq
      call get_box
      call record_pressure
      call record_volume
      call system_clock(T2,clock_rate,clock_max)
      time_data = time_data + real(T2-T1)/real(clock_rate)
    end subroutine sample_eq



    subroutine get_box
      implicit none
      logical :: exist

 
      write(box_f_id,*)timeclock,box

    end subroutine get_box

    subroutine record_pressure
      implicit none
      logical :: exist

      write(pressure_f_id,*)timeclock,pressure


    end subroutine record_pressure
    
    subroutine record_volume
      implicit none
      logical :: exist

 
      write(vol_f_id,*)timeclock,vol      

    end subroutine record_volume


    subroutine get_KE_eq
      implicit none
      logical :: exist
      double precision :: sum
      real*4 :: calcTemp
      integer :: i

      sum = 0.0d0
      do i = 1, np
         sum = sum + 0.50d0*mass(position(i)%type)*(v(i)%x*v(i)%x+v(i)%y*v(i)%y+v(i)%z*v(i)%z)
      enddo

      calctemp = 2.0d0*sum*amu2e/(kboltz*(3.0d0*np-ncons))
      
   
      write(temp_f_id,*)timeclock,calctemp
      write(KE_f_id,*)timeclock,sum

      

    end subroutine get_KE_eq

    subroutine get_PE_eq
      implicit none
      logical :: exist
      double precision :: modifyPot
      
      modifypot = potential + e_bond + e_angle + e_dihedral + e_improper + pot_14 + e_coul + e_coul_14 + e_coul_long+vdw_corr

      write(PE_f_id,*)timeclock,modifyPot
     


    end subroutine get_PE_eq


    subroutine get_totE_eq
      implicit none
      double precision :: sum
      double precision :: modifyPot
      logical :: exist

      integer :: i
      sum = 0.0d0
      do i = 1, np
         sum = sum + 0.50d0*mass(position(i)%type)*(v(i)%x*v(i)%x+v(i)%y*v(i)%y+v(i)%z*v(i)%z)
      enddo
      sum = sum*amu2e
    

      if(coul_flag.eq.2)then
         modifypot = potential + e_bond + e_angle + e_dihedral + e_improper + pot_14 + &
              e_coul + e_coul_14 +e_coul_long + vdw_corr
      elseif(coul_flag.eq.1)then
         modifypot = potential + e_bond + e_angle + e_dihedral + e_improper + pot_14 + &
              e_coul + e_coul_14 +e_self + vdw_corr

      endif

      write(total_energy_f_id,*)timeclock,modifyPot+sum

    end subroutine get_totE_eq



    subroutine record_msd
      implicit none
      character(len=*),parameter :: FMT3 = "(3f12.6)"
      integer :: i,j,off
      logical :: exist

      inquire(file=trim(path_string)//trim(data_files_path)//'/'//'msd2.dat',exist=exist)
      if(exist)then
         open(unit = 9999, file = trim(path_string)//trim(data_files_path)//'/'//'msd2.dat',status='old',position='append',action='write')
      else
         open(unit = 99999, file = trim(path_string)//trim(data_files_path)//'/'//'msd2.dat',status='new',action='write')
      endif


      msd_record = msd_record + 1
      write(9999,*)
      do i = 1,nummolarray(1)
         do j = 1,numatns

            off = (i-1)*numatns + j
            write(9999,FMT3)ptta(off)%x,ptta(off)%y,ptta(off)%z
         enddo
      enddo
      close(9999)
    end subroutine record_msd







    subroutine sample_seeded
      implicit none
      integer :: T1,T2,clock_rate,clock_max
      call system_clock(T1,clock_rate,clock_max)
      
      call get_KE_seeded
      call get_PE_seeded
      call get_totE_seeded
      call get_box_seeded
      call record_pressure_seeded
      call record_volume_seeded
      call system_clock(T2,clock_rate,clock_max)
      time_data = time_data + real(T2-T1)/real(clock_rate)
    end subroutine sample_seeded


    subroutine record_num_seeded(num)
      implicit none
      integer :: num
      logical :: exist

      inquire(file=trim(path_string)//trim(data_files_path)//'/'//trim(sim_string)//'/'//'cluster_num.dat',exist=exist)
      if(exist)then
         open(unit = 99999, file = trim(path_string)//trim(data_files_path)//'/'//trim(sim_string)//'/'//'cluster_num.dat',status='old',position='append',action='write')
      else
         open(unit = 99999, file = trim(path_string)//trim(data_files_path)//'/'//trim(sim_string)//'/'//'cluster_num.dat',status='new',action='write')
      endif
      write(99999,*)timeclock,num
      close(99999)


    end subroutine record_num_seeded



    subroutine get_box_seeded
      implicit none
      logical :: exist

      !write(8899,*)timeclock,box


      inquire(file=trim(path_string)//trim(data_files_path)//'/'//trim(sim_string)//'/'//'box_equil2.dat',exist=exist)


      if(exist)then
         open(unit = 99999, file = trim(path_string)//trim(data_files_path)//'/'//trim(sim_string)//'/'//'box_equil2.dat',status='old',position='append',action='write')
      else
         open(unit = 99999, file = trim(path_string)//trim(data_files_path)//'/'//trim(sim_string)//'/'//'box_equil2.dat',status='new',action='write')
      endif
      write(99999,*)timeclock,box
      close(99999)

    end subroutine get_box_seeded

    subroutine record_pressure_seeded
      implicit none
      logical :: exist

      !write(8887,*)timeclock,pressure


      inquire(file=trim(path_string)//trim(data_files_path)//'/'//trim(sim_string)//'/'//'pressure_equil2.dat',exist=exist)
      if(exist)then
         open(unit = 99999, file = trim(path_string)//trim(data_files_path)//'/'//trim(sim_string)//'/'//'pressure_equil2.dat',status='old',position='append',action='write')
      else
         open(unit = 99999, file = trim(path_string)//trim(data_files_path)//'/'//trim(sim_string)//'/'//'pressure_equil2.dat',status='new',action='write')
      endif
      write(99999,*)timeclock,pressure
      close(99999)


    end subroutine record_pressure_seeded
    
    subroutine record_volume_seeded
      implicit none
      logical :: exist

      !write(8896,*)timeclock,vol

      inquire(file=trim(path_string)//trim(data_files_path)//'/'//trim(sim_string)//'/'//'volume_equil2.dat',exist=exist)
      if(exist)then
         open(unit = 99999, file = trim(path_string)//trim(data_files_path)//'/'//trim(sim_string)//'/'//'volume_equil2.dat',status='old',position='append',action='write')
      else
         open(unit = 99999, file = trim(path_string)//trim(data_files_path)//'/'//trim(sim_string)//'/'//'volume_equil2.dat',status='new',action='write')
      endif
      write(99999,*)timeclock,vol
      close(99999)

      

    end subroutine record_volume_seeded


    subroutine get_KE_seeded
      implicit none
      logical :: exist
      double precision :: sum,calcTemp
      integer :: i

      sum = 0.0d0
      do i = 1, np
         sum = sum + 0.50d0*mass(position(i)%type)*(v(i)%x*v(i)%x+v(i)%y*v(i)%y+v(i)%z*v(i)%z)
      enddo

      calctemp = 2.0d0*sum*amu2e/(kboltz*(3.0d0*np-ncons))
      
      inquire(file=trim(path_string)//trim(data_files_path)//'/'//trim(sim_string)//'/'//'temp_equil2.dat',exist=exist)
      if(exist)then
         open(unit = 99999, file = trim(path_string)//trim(data_files_path)//'/'//trim(sim_string)//'/'//'temp_equil2.dat',status='old',position='append',action='write')
      else
         open(unit = 99999, file = trim(path_string)//trim(data_files_path)//'/'//trim(sim_string)//'/'//'temp_equil2.dat',status='new',action='write')
      endif
      write(99999,*)timeclock,calctemp
      close(99999)
      
      inquire(file=trim(path_string)//trim(data_files_path)//'/'//trim(sim_string)//'/'//'KE_equil2.dat',exist=exist)
      if(exist)then
         open(unit = 99999, file = trim(path_string)//trim(data_files_path)//'/'//trim(sim_string)//'/'//'KE_equil2.dat',status='old',position='append',action='write')
      else
         open(unit = 99999, file = trim(path_string)//trim(data_files_path)//'/'//trim(sim_string)//'/'//'KE_equil2.dat',status='new',action='write')
      endif
      sum = sum*amu2e
      write(99999,*)timeclock,sum
      close(99999)

      
      !print*,calctemp
      !write(10,*)timeclock,calcTemp
      !write(80,*)timeclock,calctemp
      !write(20,*)timeclock,sum

    end subroutine get_KE_seeded

    subroutine get_PE_seeded
      implicit none
      logical :: exist
      double precision :: modifyPot
      
      modifypot = potential + e_bond + e_angle + e_dihedral + e_improper + pot_14 + e_coul + e_coul_14 + e_coul_long+vdw_corr


      
      !write(30,*)timeclock,modifyPot


      inquire(file=trim(path_string)//trim(data_files_path)//'/'//trim(sim_string)//'/'//'PE_equil2.dat',exist=exist)
      if(exist)then
         open(unit = 99999, file = trim(path_string)//trim(data_files_path)//'/'//trim(sim_string)//'/'//'PE_equil2.dat',status='old',position='append',action='write')
      else
         open(unit = 99999, file = trim(path_string)//trim(data_files_path)//'/'//trim(sim_string)//'/'//'PE_equil2.dat',status='new',action='write')
      endif
      write(99999,*)timeclock,modifyPot
      close(99999)
     


    end subroutine get_PE_seeded


    subroutine get_totE_seeded
      implicit none
      double precision :: sum
      double precision :: modifyPot
      logical :: exist

      integer :: i
      sum = 0.0d0
      do i = 1, np
         sum = sum + 0.50d0*mass(position(i)%type)*(v(i)%x*v(i)%x+v(i)%y*v(i)%y+v(i)%z*v(i)%z)
      enddo
      sum = sum*amu2e
    

      if(coul_flag.eq.2)then
         modifypot = potential + e_bond + e_angle + e_dihedral + e_improper + pot_14 + &
              e_coul + e_coul_14 +e_coul_long + vdw_corr
      elseif(coul_flag.eq.1)then
         modifypot = potential + e_bond + e_angle + e_dihedral + e_improper + pot_14 + &
              e_coul + e_coul_14 +e_self + vdw_corr

      endif


      inquire(file=trim(path_string)//trim(data_files_path)//'/'//trim(sim_string)//'/'//'energy_drift2.dat',exist=exist)
      if(exist)then
         open(unit = 99999, file = trim(path_string)//trim(data_files_path)//'/'//trim(sim_string)//'/'//'energy_drift2.dat',status='old',position='append',action='write')
      else
         open(unit = 99999, file = trim(path_string)//trim(data_files_path)//'/'//trim(sim_string)//'/'//'energy_drift2.dat',status='new',action='write')
      endif
      write(99999,*)timeclock,(modifypot+sum-tot_e)/tot_e
      close(99999)


      inquire(file=trim(path_string)//trim(data_files_path)//'/'//trim(sim_string)//'/'//'total_energy2.dat',exist=exist)
      if(exist)then
         open(unit = 99999, file = trim(path_string)//trim(data_files_path)//'/'//trim(sim_string)//'/'//'total_energy2.dat',status='old',position='append',action='write')
      else
         open(unit = 99999, file = trim(path_string)//trim(data_files_path)//'/'//trim(sim_string)//'/'//'total_energy2.dat',status='new',action='write')
      endif
      write(99999,*)timeclock,modifyPot+sum
      close(99999)


    end subroutine get_totE_seeded
















    subroutine sample
      implicit none
      
      call get_KE
      call get_PE
      call get_pressure
      call get_totE

      
    end subroutine sample

    subroutine get_KE
      implicit none
      double precision :: sum,calcTemp
      integer :: i

      sum = 0.0d0
      do i = 1, np
         sum = sum + 0.50d0*mass(position(i)%type)*(v(i)%x*v(i)%x+v(i)%y*v(i)%y+v(i)%z*v(i)%z)
      enddo

      calctemp = 2.0d0*sum*amu2e/(kboltz*(3.0d0*np-ncons))

      write(50,*)timeclock,calcTemp
      write(80,*)timeclock,calctemp
      write(60,*)timeclock,sum

    end subroutine get_KE

    subroutine get_pressure
      implicit none
      
      write(100,*)pressure
    end subroutine get_pressure


    subroutine get_PE
      implicit none
      double precision :: modifyPot
      
      modifypot  = (potential+ e_bond + e_angle + e_dihedral + e_improper + pot_14)
      write(70,*)timeclock,modifyPot
    end subroutine get_PE


    subroutine get_totE
      implicit none
      double precision :: sum
      double precision :: modifyPot
      integer :: i
      sum = 0.0d0
      do i = 1, np
         sum = sum + 0.50d0*mass(position(i)%type)*(v(i)%x*v(i)%x+v(i)%y*v(i)%y+v(i)%z*v(i)%z)
      enddo
      sum = sum*amu2e
      modifyPot = potential + e_bond + e_angle + e_dihedral + e_improper + pot_14

      write(40,*)timeclock,(modifyPot+sum)
    end subroutine get_totE


    subroutine sample_anneal
      integer :: T1,T2,clock_rate,clock_max
      
      call system_clock(T1,clock_rate,clock_max)
      call get_temp_anneal
      call get_potential_anneal
      call system_clock(T2,clock_rate,clock_max)
      
      time_data = time_data + real(T2-T1)/real(clock_rate)

    end subroutine sample_anneal

    subroutine get_temp_anneal
      double precision :: sum,calcTemp
      integer :: i

      sum = 0.0d0
      do i = 1, np
         sum = sum + 0.50d0*mass(position(i)%type)*(v(i)%x*v(i)%x+v(i)%y*v(i)%y+v(i)%z*v(i)%z)
      enddo
      sum = sum*amu2e
      
      calcTemp = 2.0d0*sum/(kboltz*(3.0d0*np-ncons))
      
      write(8975,*)calcTemp
    end subroutine get_temp_anneal


    subroutine get_potential_anneal
      implicit none
      double precision :: modifyPot
      
      modifypot = potential + e_bond + e_angle + e_dihedral + e_improper + pot_14 + e_coul + e_coul_14


      
      write(8001,*)modifyPot
    end subroutine get_potential_anneal



    subroutine sample_grocluster(flag)
      integer :: T1,T2,clock_rate,clock_max
      integer :: flag
      call system_clock(T1,clock_rate,clock_max)
      call get_temp_grocluster(flag)
      call get_potential_grocluster
      call get_box_gro
      call record_pressure_gro
      call record_volume_gro
      call system_clock(T2,clock_rate,clock_max)
      

    end subroutine sample_grocluster

    subroutine get_box_gro
      implicit none
      logical :: exist
      !write(88990,*)timeclock,box

      inquire(file=trim(path_string)//trim(data_files_path)//'/'//'box_groclus2.dat',exist=exist)
      if(exist)then
         open(unit = 99999, file = trim(path_string)//trim(data_files_path)//'/'//'box_groclus2.dat',status='old',position='append',action='write')
      else
         open(unit = 99999, file = trim(path_string)//trim(data_files_path)//'/'//'box_groclus2.dat',status='new',action='write')
      endif
      write(99999,*)timeclock,box
      close(99999)
         
    end subroutine get_box_gro

    subroutine record_pressure_gro
      implicit none
      logical :: exist
      !write(88870,*)timeclock,pressure

      inquire(file=trim(path_string)//trim(data_files_path)//'/'//'pressure_groclus2.dat',exist=exist)
      if(exist)then
         open(unit = 99999, file = trim(path_string)//trim(data_files_path)//'/'//'pressure_groclus2.dat',status='old',position='append',action='write')
      else
         open(unit = 99999, file = trim(path_string)//trim(data_files_path)//'/'//'pressure_groclus2.dat',status='new',action='write')
      endif
      write(99999,*)timeclock,pressure
      close(99999)


    end subroutine record_pressure_gro
    
    subroutine record_volume_gro
      implicit none
      logical :: exist
      !write(88960,*)timeclock,vol

      inquire(file=trim(path_string)//trim(data_files_path)//'/'//'vol_groclus2.dat',exist=exist)
      if(exist)then
         open(unit = 99999, file = trim(path_string)//trim(data_files_path)//'/'//'vol_groclus2.dat',status='old',position='append',action='write')
      else
         open(unit = 99999, file = trim(path_string)//trim(data_files_path)//'/'//'vol_groclus2.dat',status='new',action='write')
      endif
      write(99999,*)timeclock,vol
      close(99999)


    end subroutine record_volume_gro


    subroutine get_temp_grocluster(flag)
      double precision :: sum,calcTemp
      integer :: i,flag
      logical :: exist

      flag = 0
      sum = 0.0d0
      do i = 1, np
         sum = sum + 0.50d0*mass(position(i)%type)*(v(i)%x*v(i)%x+v(i)%y*v(i)%y+v(i)%z*v(i)%z)
      enddo
      sum = sum*amu2e

      calcTemp = 2.0d0*sum/(kboltz*(3.0d0*np-ncons))
      


      !write(9999,*)calcTemp

      inquire(file=trim(path_string)//trim(data_files_path)//'/'//'temp_groclus2.dat',exist=exist)
      if(exist)then
         open(unit = 99999, file = trim(path_string)//trim(data_files_path)//'/'//'temp_groclus2.dat',status='old',position='append',action='write')
      else
         open(unit = 99999, file = trim(path_string)//trim(data_files_path)//'/'//'temp_groclus2.dat',status='new',action='write')
      endif
      write(99999,*)timeclock,calcTemp
      close(99999)
      
    end subroutine get_temp_grocluster


    subroutine get_potential_grocluster
      implicit none
      logical :: exist
      double precision :: modifyPot
      
      modifypot = potential + e_bond + e_angle + e_dihedral + e_improper + pot_14 + e_coul + e_coul_14

      !write(9991,*)modifyPot
      inquire(file=trim(path_string)//trim(data_files_path)//'/'//'potential_groclus2.dat',exist=exist)
      if(exist)then
         open(unit = 99999, file = trim(path_string)//trim(data_files_path)//'/'//'potential_groclus2.dat',status='old',position='append',action='write')
      else
         open(unit = 99999, file = trim(path_string)//trim(data_files_path)//'/'//'potential_groclus2.dat',status='new',action='write')
      endif
      write(99999,*)timeclock,modifypot
      close(99999)


    end subroutine get_potential_grocluster




  end module mod_sample

