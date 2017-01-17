module mod_makefolders
  use global
  contains
    
    subroutine make_MPI_files
      
      implicit none
      integer, parameter :: rk = selected_real_kind(p=15), ik=4
      logical :: e,result
      integer :: i
      integer :: Oh,j
      integer :: disp,loc
      integer :: count
      character(len=150) :: folder,folder2
      character(len=40)  :: file
      character(len=200) :: filename,num

      print*,'in make mpi files'
      INQUIRE(file = 'DATA'//trim(rank_char),EXIST=e)
      if(.not.e)then
         call system('mkdir '//'DATA')
      endif

      INQUIRE(file = 'DATA'//trim(rank_char),EXIST=e)
      if(.not.e)then
         call system('mkdir '//'DATA'//'/'//trim(rank_char))
      endif
      
      grofile_path = 'DATA'//'/'//trim(rank_char)//'/'//'grofiles'
      grofile_path = ADJUSTL(grofile_path)
      grofile_path = trim(grofile_path)
      grofile_path = adjustl(grofile_path)
      
      INQUIRE(file = trim(grofile_path),EXIST=e)
      if(.not.e)then
         call system('mkdir '//trim(grofile_path))
         !result = MakeDirQQ(trim(folder))
      endif
      
      data_files_path = 'data_files'
      data_files_path = ADJUSTL(data_files_path)
      data_files_path = 'DATA'//'/'//trim(rank_char)//'/'//data_files_path
      data_files_path = trim(data_files_path)
      data_files_path = adjustl(data_files_path)
      !folder = trim(path)//'/'//folder
      
      INQUIRE(file = trim(data_files_path),EXIST=e)
      if(.not.e)then
         call system('mkdir '//trim(data_files_path))
         !result = MakeDirQQ(trim(folder))
      endif

      INQUIRE(file = 'MPIrestart',EXIST=e)
      if(.not.e)then
         call system('mkdir '//'MPIrestart')
         !result = MakeDirQQ(trim(folder))
      endif

      restart_files_path = trim('MPIrestart')//'/'//trim(rank_char)
      INQUIRE(file = trim(restart_files_path),EXIST=e)
      if(.not.e)then
         call system('mkdir '//trim(restart_files_path))
         !result = MakeDirQQ(trim(folder))
      endif


      INQUIRE(file = 'MPIrestart_final',EXIST=e)
      if(.not.e)then
         call system('mkdir '//'MPIrestart_final')
         !result = MakeDirQQ(trim(folder))
      endif

      restart_files_path_final = trim('MPIrestart_final')//'/'//trim(rank_char)
      INQUIRE(file = trim(restart_files_path),EXIST=e)
      if(.not.e)then
         call system('mkdir '//trim(restart_files_path_final))
         !result = MakeDirQQ(trim(folder))
      endif
      
      !---open record file
      write(num,'(i8)')restart_mult
      num = adjustl(num)
      num = trim(num)
      filename = 'log.sofia'//'_'//trim(num)
      open(unit = 1234, file = 'DATA'//'/'//trim(rank_char)//'/'//trim(filename))
      print*,'leaving make mpi files'
    end subroutine make_MPI_files

  END module mod_makefolders
