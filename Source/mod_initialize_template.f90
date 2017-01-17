module mod_initialize_template
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
 

contains


  subroutine initialize_template
    implicit none
    double precision :: rcut0,rcut20
    integer :: i,j,k,res
    character(len=50) :: exec,file,newton,auto,rms,filename


   
  end subroutine initialize_template


end module mod_initialize_template
