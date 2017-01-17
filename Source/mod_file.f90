module mod_file
  use global
  
contains
  
  subroutine open_sample_files
    implicit none
    logical :: exist

    temp_f_id = 99990
    KE_f_id   = 99991
    PE_f_id   = 99992
    total_energy_f_id = 99993
    clus_f_id = 99994
    box_f_id = 99995
    vol_f_id = 99996
    pressure_f_id= 99997
    
    
    inquire(file=trim(path_string)//trim(data_files_path)//'/'//'temp_equil2.dat',exist=exist)
    if(exist)then
       open(unit = temp_f_id, file = trim(path_string)//trim(data_files_path)//'/'//'temp_equil2.dat',status='old',position='append',action='write')
    else
       open(unit = temp_f_id, file = trim(path_string)//trim(data_files_path)//'/'//'temp_equil2.dat',status='new',action='write')
    endif

    inquire(file=trim(path_string)//trim(data_files_path)//'/'//'KE_equil2.dat',exist=exist)
    if(exist)then
       open(unit = KE_f_id, file = trim(path_string)//trim(data_files_path)//'/'//'KE_equil2.dat',status='old',position='append',action='write')
    else
       open(unit = KE_f_id, file = trim(path_string)//trim(data_files_path)//'/'//'KE_equil2.dat',status='new',action='write')
    endif

    inquire(file=trim(path_string)//trim(data_files_path)//'/'//'potential_equil2.dat',exist=exist)
    if(exist)then
       open(unit = PE_f_id, file = trim(path_string)//trim(data_files_path)//'/'//'potential_equil2.dat',status='old',position='append',action='write')
    else
       open(unit = PE_f_id, file = trim(path_string)//trim(data_files_path)//'/'//'potential_equil2.dat',status='new',action='write')
    endif


    inquire(file=trim(path_string)//trim(data_files_path)//'/'//'total_energy.dat',exist=exist)
    if(exist)then
       open(unit = total_energy_f_id, file = trim(path_string)//trim(data_files_path)//'/'//'total_energy.dat',status='old',position='append',action='write')
    else
       open(unit = total_energy_f_id, file = trim(path_string)//trim(data_files_path)//'/'//'total_energy.dat',status='new',action='write')
    endif



    inquire(file=trim(path_string)//trim(data_files_path)//'/'//'total_energy.dat',exist=exist)
    if(exist)then
       open(unit = clus_f_id, file = trim(path_string)//trim(data_files_path)//'/'//'total_energy.dat',status='old',position='append',action='write')
    else
       open(unit = clus_f_id, file = trim(path_string)//trim(data_files_path)//'/'//'total_energy.dat',status='new',action='write')
    endif



    inquire(file=trim(path_string)//trim(data_files_path)//'/'//'cluster_num.dat',exist=exist)
    if(exist)then
       open(unit = clus_f_id, file = trim(path_string)//trim(data_files_path)//'/'//'cluster_num.dat',status='old',position='append',action='write')
    else
       open(unit = clus_f_id, file = trim(path_string)//trim(data_files_path)//'/'//'cluster_num.dat',status='new',action='write')
    endif


    inquire(file=trim(path_string)//trim(data_files_path)//'/'//'volume_equil2.dat',exist=exist)
    if(exist)then
       open(unit = vol_f_id, file = trim(path_string)//trim(data_files_path)//'/'//'volume_equil2.dat',status='old',position='append',action='write')
    else
       open(unit = vol_f_id, file = trim(path_string)//trim(data_files_path)//'/'//'volume_equil2.dat',status='new',action='write')
    endif



    inquire(file=trim(path_string)//trim(data_files_path)//'/'//'box_equil2.dat',exist=exist)
    if(exist)then
       open(unit = box_f_id, file = trim(path_string)//trim(data_files_path)//'/'//'box_equil2.dat',status='old',position='append',action='write')
    else
       open(unit = box_f_id, file = trim(path_string)//trim(data_files_path)//'/'//'box_equil2.dat',status='new',action='write')
    endif


    inquire(file=trim(path_string)//trim(data_files_path)//'/'//'pressure_equil2.dat',exist=exist)
    if(exist)then
       open(unit = pressure_f_id, file = trim(path_string)//trim(data_files_path)//'/'//'pressure_equil2.dat',status='old',position='append',action='write')
    else
       open(unit = pressure_f_id, file = trim(path_string)//trim(data_files_path)//'/'//'pressure_equil2.dat',status='new',action='write')
    endif


  end subroutine open_sample_files

  subroutine open_anneal_files
    implicit none
    
    open(unit = 8975, file = 'temp_anneal.dat')
    open(unit = 8001, file = 'potential_anneal.dat')
  end subroutine open_anneal_files



  subroutine close_files
    implicit none
    print*,'closing files:'

    
    close(temp_f_id)
    close(KE_f_id)
    close(PE_f_id)
    close(total_energy_f_id)
    close(clus_f_id)
    close(box_f_id)
    close(vol_f_id)
    close(pressure_f_id)



  end subroutine close_files
  
end module mod_file
