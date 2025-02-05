!include 'mkl_vsl.f90'
program MD
  !USE MKL_VSL_TYPE
  use global 
  use mod_file
  use mod_initialize
  use mod_initialize_template
  use mod_restart
  use mod_integrate
  use mod_anneal
  use mod_nucleation
  use mod_restart_datadump
  use mpi
  implicit none
  integer :: i
  !integer :: T1,T2,clock_rate,clock_max
  integer :: ierr,msgtag,n0temp
  integer :: status(MPI_STATUS_SIZE)
  integer :: lnodename, node_comm, n, color,color2, key, rank_in_group
  character (len=1) :: separator = '.'
  character (len=50) :: nodename



  !******initialize MPI and rank_char
  call mpi_init(ierr)
  call mpi_comm_size(MPI_COMM_WORLD,size,ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)
  
  
 
  
  i = rank
  ourmic = mod(rank,2)
  
  print*,'hello from rank',rank,size

  write(rank_char,'(i8)')rank
  rank_char = ADJUSTL(rank_char)
  rank_char = trim(rank_char)
  
  
  
  call srand(rank)
  open(unit = 1000, file = 'input')
  read(1000,*)
  read(1000,*)
  read(1000,*)
  read(1000,*)
  read(1000,*)sta
  close(1000)
  if(sta.eq.'TEM')then
     call initialize_template
  elseif(sta.eq.'XYZ')then
     call restart
  elseif(sta.eq.'DUM')then
     call restart_datadump
  else
     call initialize
  endif

  !call system_clock(T1,clock_rate,clock_max)


  if(ensemble.eq.'NVE')then
     call integrate_npt
  endif

  if(ensemble.eq.'NUC')then
     call anneal
  endif

  if(ensemble.eq.'NRS')then
     call grow_cluster(100,1)     
  endif

  call close_files

  call mpi_finalize(ierr)


  !call system_clock(T2,clock_rate,clock_max)

  !print*,'TOTAL ELAPSED TIME IN MD PROGRAM:',rank,real(T2-T1)/real(clock_rate)
  
end program MD
