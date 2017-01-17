module mod_datadump
  use global
  use mod_memory
  contains
    subroutine datadump_file(filename)
      double precision :: x,y,z,c
      integer :: i,typ,J
      integer :: a1,a2,a3
      integer :: counter
      character(len=3) :: char
      character(len=*),parameter :: FM1 = "(i6,a5,f12.6,2f12.6,5i6)"
      character(len=*),parameter :: FM2 = "(i6)"
      character(len=*) :: filename


      open(unit = 500,file = filename )

      write(500,*)np
      write(500,*)numMolType
      do i=1,numMolType
         write(500,*)numMolArray(i)
      enddo
      do i = 1,np
         typ  = position(i)%type
         x = position(i)%x
         y = position(i)%y
         z = position(i)%z
         c = q(i)
         
         char = type2string(typ)
         char = adjustl(char)
         char = trim(char)

        ! typ  = ltype2ftype(typ)
         
         write(500,FM1)i,char,x,y,z,typ,num1bond(i),num2bond(i),num3bond(i)
         write(500,*)c
         counter = 0
         do while(counter.lt.num3bond(i))
            counter = counter+1
            write(500,*)specbond(counter,i)
         enddo         
      end do
      print*,'WE HAVE DONE A DATA DUMP'
      close(500)
    end subroutine datadump_file


    subroutine partial_read_in_datadump(filename)
      double precision :: x,y,z
      integer :: i,typ,J,junk
      integer :: a1,a2,a3
      integer :: counter
      character(len=3) :: char
      character(len=*),parameter :: FM1 = "(i6,a5,f12.6,2f12.6,5i6)"
      character(len=*),parameter :: FM2 = "(i6)"
      character(len=*) :: filename


      open(unit = 500,file = filename )

      read(500,*)np
      read(500,*)numMolType
      do i =1,numMolType
         read(500,*)numMolArray(i)
      enddo

      close(500)
    end subroutine partial_read_in_datadump



    subroutine dump_global_id(filename)
      implicit none
      integer :: i,j
      double precision :: x,y,z
      character(len=*) :: filename

      open(unit = 5001,file=filename)
      write(5001,*)np

      do i = 1,np
         j = global_id(i)
         x = position(i)%x; y = position(i)%y; z=position(i)%z
         write(5001,*)j,x,y,z
      enddo
      close (5001)
    end subroutine dump_global_id


    subroutine dump_global_velocity_id(filename)
      implicit none
      character(len=*) :: filename
      integer :: i,j
     ! real*4:: x,y,z,vx,vy,vz
      double precision :: x,y,z,vx,vy,vz
      character(len=*),parameter :: FMT3 = "(i5,6f12.6)"


      open(unit = 5001, file = filename)
      write(5001,*)box
      do i = 1,np
         j = global_id(i)
         x = position(i)%x; y = position(i)%y; z= position(i)%z
         vx = v(i)%x; vy = v(i)%y; vz = v(i)%z
         
         write(5001,FMT3)j,x,y,z,vx,vy,vz
      enddo
      close(5001)
    end subroutine dump_global_velocity_id


    subroutine dump_global_velocity_id_mpi(step,num)
      implicit none
      character(len=50) :: filename,ints
      integer :: i,j,step,num
      real*4:: x,y,z,vx,vy,vz
      character(len=*),parameter :: FMT3 = "(i5,6f12.6)"


     
      call flush(temp_f_id)
      call flush(KE_f_id)
      call flush(PE_f_id)
      call flush(total_energy_f_id)
      call flush(clus_f_id)
      call flush(box_f_id)
      call flush(vol_f_id)
      call flush(pressure_f_id)

      

      write(ints,'(i8)')step
      ints = adjustl(ints)
      ints = trim(ints)
      filename= trim(restart_files_path)//'/'//'restart_'
      filename = adjustl(filename)
      filename = trim(filename)
      filename = trim(filename)//trim(ints)
      open(unit = 5001,status='replace',file = filename)
      write(5001,*)box
      write(5001,*)temp
      write(5001,*)num
      write(5001,*)timeclock
      do i = 1,np
         j = global_id(i)
         x = position(i)%x; y = position(i)%y; z= position(i)%z
         vx = v(i)%x; vy = v(i)%y; vz = v(i)%z
         
         write(5001,FMT3)j,x,y,z,vx,vy,vz
      enddo
      close(5001)
    end subroutine dump_global_velocity_id_mpi


    subroutine dump_global_velocity_id_mpi_final(step,num)
      implicit none
      character(len=50) :: filename,ints
      integer :: i,j,step,num
      real*4:: x,y,z,vx,vy,vz
      character(len=*),parameter :: FMT3 = "(i5,6f12.6)"


      call flush(temp_f_id)
      call flush(KE_f_id)
      call flush(PE_f_id)
      call flush(total_energy_f_id)
      call flush(clus_f_id)
      call flush(box_f_id)
      call flush(vol_f_id)
      call flush(pressure_f_id)




      write(ints,'(i8)')step
      ints = adjustl(ints)
      ints = trim(ints)
      filename= trim(restart_files_path_final)//'/'//'restart_'
      filename = adjustl(filename)
      filename = trim(filename)
      filename = trim(filename)//trim(ints)
      open(unit = 5001, file = filename)
      write(5001,*)box
      write(5001,*)temp
      write(5001,*)num
      write(5001,*)timeclock
      do i = 1,np
         j = global_id(i)
         x = position(i)%x; y = position(i)%y; z= position(i)%z
         vx = v(i)%x; vy = v(i)%y; vz = v(i)%z
         
         write(5001,FMT3)j,x,y,z,vx,vy,vz
      enddo
      close(5001)
    end subroutine dump_global_velocity_id_mpi_final

    subroutine read_dump_global_id(filename)
      implicit none
      integer :: i,j
      double precision :: x,y,z
      character(len=*) :: filename

      open(unit = 5001,file=filename)
      read(5001,*)np

      do i = 1,np
         read(5001,*)j,x,y
         read(5001,*)z

         position(j)%x=x; position(j)%y=y; position(j)%z=z
      enddo
    end subroutine read_dump_global_id


    subroutine read_dump_global_id_v(filename)
      implicit none
      integer :: i,j,read_num
      real*4  :: x,y,z,vx,vy,vz,junk
      character(len=*) :: filename
      
      print*,'TIME TO OPEN RESTART FILE'
      open(unit = 5001,file=filename)
      print*,'OPENED RESTART FILE'
      read(5001,*)box
      read(5001,*)junk
      read(5001,*)read_num
      read(5001,*)timeclock
      hbox = box/2.0
      ibox = 1.0/box
      vol = box*box*box

      print*,'========================================='
      print*,'READIN IN COORDINATES:',np
      PRINT*,'we have read in our box size',box
      print*,'filename is',trim(filename)
      print*,'cluster size:',read_num
      print*
      print*
      print*,'=========================================='
      do i = 1,np
         read(5001,*)j,x,y,z,vx,vy,vz
         !read(5001,*)vy,vz

         position(j)%x=x; position(j)%y=y; position(j)%z=z
         v(j)%x = vx; v(j)%y = vy; v(j)%z = vz
         
         !v(i)%x = vx; v(i)%y = vy; v(i)%z = 
      enddo
      close(5001)
      
    end subroutine read_dump_global_id_v


    subroutine read_in_datadump(filename)
      double precision :: x,y,z,c
      integer :: i,typ,J,junk
      integer :: a1,a2,a3
      integer :: counter
      character(len=3) :: char
      character(len=*),parameter :: FM1 = "(i6,a5,f12.6,2f12.6,5i6)"
      character(len=*),parameter :: FM2 = "(i6)"
      character(len=*) :: filename


     open(unit = 500,file = filename )

      read(500,*)np
      read(500,*)numMolType
      do i =1,numMolType
         read(500,*)numMolArray(i)
      enddo

      do i =1,np
         read(500,FM1)junk,char,x,y,z,typ,num1bond(i),num2bond(i),num3bond(i)
         read(500,*)c
         position(i)%x = x
         position(i)%y = y
         position(i)%z = z
         position(i)%type = typ
         q(i) = c

        

         do counter = 1,num3bond(i)
            read(500,*)specbond(counter,i)
         enddo
      enddo
         
       
      print*,'WE HAVE READ IN A DATA DUMP'
      close(500)
    end subroutine read_in_datadump

    
    subroutine datadump_bondlist(filename)
      implicit none
      integer :: step
      integer :: i
      integer :: m,i1,i2,itype,ifactor
      double precision:: delx,dely,delz,rsq,r,dr,rk,force
      character(len=*) :: filename
      
      open(unit = 500, file = filename)
      write(500,*)numbond
      do m = 1,numbond
         
         i1 = bondlist(1,m)
         i2 = bondlist(2,m)
         itype = bondlist(3,m)
         write(500,*)i1,i2,itype
      end do
      close(500)
      
    end subroutine datadump_bondlist


    subroutine read_in_bondlist(filename)
      implicit none
      integer :: step
      integer :: i
      integer :: m,i1,i2,itype,ifactor
      double precision:: delx,dely,delz,rsq,r,dr,rk,force
      character(len=*) :: filename
      
      print*,'reading in bondlist'
      open(unit = 500, file = filename)
      read(500,*)numbond
      call totalbond_memory
      do m = 1,numbond
         
         i1 = bondlist(1,m)
         i2 = bondlist(2,m)
         itype = bondlist(3,m)
         read(500,*)bondlist(1,m),bondlist(2,m),bondlist(3,m)
      end do
      close(500)
      
    end subroutine read_in_bondlist
    
    
    subroutine datadump_anglelist(filename)
      implicit none
      integer :: i1,i2,i3,itype
      integer :: m
      character(len=*) :: filename
      

      open(unit = 500, file = filename)
      write(500,*)numangle
      do m = 1,numangle
         i1 = anglelist(1,m)
         i2 = anglelist(2,m)
         i3 = anglelist(3,m)
         itype = anglelist(4,m)
         write(500,*)i1,i2,i3,itype
      enddo
      close(500)
    end subroutine datadump_anglelist


    subroutine read_in_anglelist(filename)
      implicit none
      integer :: i1,i2,i3,itype
      integer :: m
      character(len=*) :: filename
      
      
      print*,'reading in angle list'
      open(unit = 500, file = filename)
      read(500,*)numangle
      call totalangle_memory
      do m = 1,numangle
         read(500,*)anglelist(1,m),anglelist(2,m),anglelist(3,m),anglelist(4,m)
      enddo
      close(500)
    end subroutine read_in_anglelist

       
    subroutine datadump_dihedlist(filename)
      implicit none
      integer :: i1,i2,i3,i4,itype
      integer :: m
      character(len=*) :: filename
      

      open(unit = 500, file = filename)
      write(500,*)numdihed
      do m = 1,numangle
         i1 = dihedlist(1,m)
         i2 = dihedlist(2,m)
         i3 = dihedlist(3,m)
         i4 = dihedlist(4,m)
         itype = dihedlist(5,m)
         write(500,*)i1,i2,i3,i4,itype
      enddo
      close(500)
    end subroutine datadump_dihedlist

    subroutine read_in_dihedlist(filename)
      implicit none
      integer :: i1,i2,i3,i4,itype
      integer :: m
      character(len=*) :: filename
      
      print*,'reading in dihedlist'

      open(unit = 500, file = filename)
      read(500,*)numdihed
      call totaldihed_memory
      do m = 1,numdihed
         read(500,*)i1,i2,i3,i4,itype
         dihedlist(1,m) = i1
         dihedlist(2,m) = i2
         dihedlist(3,m) = i3
         dihedlist(4,m) = i4
         dihedlist(5,m) = itype
      enddo
      close(500)
    end subroutine read_in_dihedlist
           
    subroutine datadump_improlist(filename)
      implicit none
      integer :: i1,i2,i3,i4,itype
      integer :: m
      character(len=*) :: filename
      

      open(unit = 500, file = filename)
      write(500,*)numimpro
      do m = 1,numimpro
         i1 = improlist(1,m)
         i2 = improlist(2,m)
         i3 = improlist(3,m)
         i4 = improlist(4,m)
         itype = improlist(5,m)
         write(500,*)i1,i2,i3,i4,itype
      enddo
      close(500)
    end subroutine datadump_improlist
    
    subroutine read_in_improlist(filename)
      implicit none
      integer :: i1,i2,i3,i4,itype
      integer :: m
      character(len=*) :: filename
      
      print*,'reading in impro'
      open(unit = 500, file = filename)
      read(500,*)numimpro
      call improdihed_memory
      do m = 1,numimpro
         read(500,*)i1,i2,i3,i4,itype
         improlist(1,m) = i1
         improlist(2,m) = i2
         improlist(3,m) = i3
         improlist(4,m) = i4
         improlist(5,m) = itype
      enddo
      close(500)
    end subroutine read_in_improlist




    subroutine dump_mollist(filename)
      implicit none
      character(len=*) :: filename
      integer :: i,j,k

     open(unit = 500, file = filename)
     do i =1,np
        write(500,*)mol(1,i),mol(2,i)
     enddo
     close(500)

   end subroutine dump_mollist
    

    subroutine read_in_mollist(filename)
      implicit none
      character(len=*) :: filename
      integer :: i,j,k

     open(unit = 500, file = filename)
     do i =1,np
        read(500,*)mol(1,i),mol(2,i)
     enddo
     close(500)

   end subroutine read_in_mollist



   subroutine dump_ptta(filename)
     implicit none
     character(len=*) :: filename
     character(len=*),parameter :: FMT3 = "(3f12.6)"
     integer :: i,j
     
     open(unit = 500, file = filename)
     do i = 1,numMolArray(1)*numatns
        write(500,FMT3)ptta(i)%x,ptta(i)%y,ptta(i)%z
     enddo
     close(500)
   end subroutine dump_ptta


    subroutine read_in_ptta(filename)
     implicit none
     character(len=*) :: filename
     integer :: i,j
     
     print*,'reading in ptta'
     print*,'check',numMolArray(1),numatns
     open(unit = 500, file = filename)
     do i = 1,numMolArray(1)*numatns
        read(500,*)ptta(i)%x,ptta(i)%y,ptta(i)%z
     enddo
     close(500)
   end subroutine read_in_ptta



   subroutine read_in_velocity(filename)
     implicit none
     character(len=*) :: filename
     integer :: i,j,k
     real*4 :: vx,vy,vz

     open(unit = 500, file = filename)

     do i = 1,np
        read(500,*)vx,vy,vz
        v(i)%x = vx
        v(i)%y = vy
        v(i)%z = vz
     enddo
     close(500)
        
   end subroutine read_in_velocity

    subroutine datadump_velocity(filename)
     implicit none
     character(len=*) :: filename
     integer :: i,j,k
     real*4 :: vx,vy,dz

     open(unit = 500, file = filename)

     do i = 1,np
        write(500,*)v(i)%x,v(i)%y,v(i)%z
     enddo
     close(500)
        
   end subroutine datadump_velocity


  end module mod_datadump
