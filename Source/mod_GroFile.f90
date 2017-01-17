module mod_GroFile
  use global
  
  implicit none
  contains

    subroutine template_generator(part,filename,filename2)
      implicit none
      integer :: part,neigh,i
      double precision :: x,y,z
      integer :: length_clus,type,cluster,numAt
      character(len=*),parameter :: FM1 = "(i5)"
      character(len=*),parameter :: FMT2 = "( i5 , 2a5 , i5 , 3f8.3)"
      character(len=*),parameter :: FMT3 = "(f10.5,f10.5,f10.5)"
      character(len=*) :: filename,filename2
      integer :: loop,typ,id,j,jflag
      character(len=3) :: char
      integer :: T1,T2,clock_rate,clock_max,count


      open(unit = 1000, file = filename)
      write(1000,*)'Lennard Jones'
      write(1000,FM1)denstrack(part)*10+10

      open(unit = 1001, file =filename2)

      count = 0
      do i =(part-1)*10+1,(part-1)*10+10
            
         typ  = position(i)%type
         char = type2string(typ)
         char = adjustl(char)
         char = trim(char)
         
         typ  = ltype2ftype(typ)
         
         x = position(i)%x/10.0d0; y = position(i)%y/10.0d0; z = position(i)%z/10.0d0

         write(1001,*)position(i)%x,position(i)%y,position(i)%z
         
         count = count + 1
         write(1000,FMT2)count,'TTA',char,count,x,y,z
      enddo

      do i = 1,denstrack(part)
         neigh = densvl( (part-1)*densvlalloc + i )
         
         do j = (neigh-1)*10+1,(neigh-1)*10+10
            
            typ  = position(j)%type
            char = type2string(typ)
            char = adjustl(char)
            char = trim(char)
            
            typ  = ltype2ftype(typ)
            
            x = position(j)%x/10.0d0; y = position(j)%y/10.0d0; z = position(j)%z/10.0d0
            
            write(1001,*)position(j)%x,position(j)%y,position(j)%z
            count = count + 1
            write(1000,FMT2)count,'TTA',char,count,x,y,z
         enddo
      enddo

      write(1000,FMT3)box/10.0d0,box/10.0d0,box/10.0d0
      close(1000)
    end subroutine template_generator

    subroutine read_in_test(filename)
      implicit none
      integer :: i,j
      double precision :: x,y,z
      character(len=*) :: filename
      
      open(unit = 500, file = filename)
   
      do i = 1,np
         read(500,*)j,x,y,z
         !read(500,*)x,y,z
         position(j)%x = x
         position(j)%y = y
         position(j)%z = z

         
      end do
      close(500)
    end subroutine read_in_test
      
    subroutine read_in_Grofile_box(filename)
      implicit none
      character(len=30) :: jum
      character(len=3) :: char
      character(len=3) :: char2
      character(len=*) :: filename
      character(len=*),parameter :: FMT3 = "(f10.5,f10.5,f10.5)"
      real*4 :: x,y,z,dummy
      integer :: i,j,k,flag

      print*,'WE ARE READING IN GROFILE'
      open(unit = 5000, file = filename)
      read(5000,*)
      read(5000,*)j
      print*,'what is j',j
      do i = 1,j
         
         read(5000,*)
      enddo
      read(5000,FMT3)box
      box = box*10.0d0

      close(5000)
   

         
    end subroutine read_in_Grofile_box
    
     subroutine read_in_Grofile(filename)
      implicit none
      character(len=30) :: jum
      character(len=3) :: char
      character(len=3) :: char2
      character(len=*) :: filename
      real*4 :: x,y,z,dummy
      integer :: i,j,k,flag

      print*,'WE ARE READING IN GROFILE'
      open(unit = 5000, file = filename)
      read(5000,*)
      read(5000,*)j
        do i = 1,np
         if(i.le.9999.or.i.gt.99999.and.i.le.109999)then
            read(5000,*)char,char2,k,position(i)%x,position(i)%y,position(i)%z
            
            position(i)%x = position(i)%x*10.0d0
            position(i)%y = position(i)%y*10.0d0
            position(i)%z = position(i)%z*10.0d0
            
           ! if(i.eq.1000000)then
           ! print*,'read in',pcharosition(i)%x,position(i)%y,position(i)%z
           
          elseif(i.gt.9999.or.i.gt.109999)then
             read(5000,*)char,char2,position(i)%x,position(i)%y,position(i)%z
           
             if(i.eq.9999)then
                
                print*,'trouble line dont forget!'
                print*,position(i)%x,position(i)%y,position(i)%z
                stop
             endif
            position(i)%x = position(i)%x*10.0d0
            position(i)%y = position(i)%y*10.0d0
            position(i)%z = position(i)%z*10.0d0
         endif
      enddo
      close(5000)
   

         
    end subroutine read_in_Grofile







    subroutine Grofile_prune(filename)
      implicit none
      double precision :: x,y,z,box_units
      character(len=*) :: filename
      character(len=*),parameter :: FM1 = "(i5)"
      character(len=*),parameter :: FMT2 = "( i5 , 2a5 , i5 , 3f8.3)"
      character(len=*),parameter :: FMT3 = "(f10.5,f10.5,f10.5)"
      integer :: i,flag,spot,j
      integer :: typ,type
      character(len=3) :: char

      open(unit = 7000, file = trim(grofile_path)//'/'//filename)
      write(7000,*)'xyz'

      write(7000,FM1)numMolArray(1)*numatompertemp
  
      do i = 1,numMolArray(1)
         do j = 1,numatompertemp
            
            spot = (i-1)*10+j
            typ  = position(spot)%type
            char = type2string(typ)
            char = adjustl(char)
            char = trim(char)          
            
            spot = (i-1)*numatompertemp+j
            x = pruned(spot)%x/10.0d0
            y = pruned(spot)%y/10.0d0
            z = pruned(spot)%z/10.0d0
            
            
            write(7000,FMT2)i,'TTA',char,i, x,y,z
         enddo
      end do

      box_units = box/(10.0d0)
      
      write(7000,FMT3)box_units,box_units,box_units
      close(7000)

    end subroutine Grofile_prune

    subroutine debug_gro(xx,filename)
      double precision :: box_units
      double precision :: x,y,z
      double precision :: xx(3,numatompertemp)
      character(len=*) :: filename
      character(len=*),parameter :: FM1 = "(i5)"
      character(len=*),parameter :: FMT2 = "( i5 , 2a5 , i5 , 3f8.3)"
      character(len=*),parameter :: FMT3 = "(f10.5,f10.5,f10.5)"
      integer :: i,flag
      integer :: typ,type
      character(len=3) :: char

      open(unit = 7000, file = trim(grofile_path)//'/'//filename)
      write(7000,*)'xyz'

      write(7000,FM1)numatompertemp
  
      do i = 1,numatompertemp
      
         if(i.eq.3)char='N'
         if(i.eq.2)char='CA'
         if(i.eq.1)char='C'
         if(i.eq.4)char='O'
         x  = xx(1,i)/10.0d0
         y  = xx(2,i)/10.0d0
         z  = xx(3,i)/10.0d0
         
         write(7000,FMT2)i,'TTA',char,i, x,y,z
      enddo



      box_units = box/(10.0d0)
      
      write(7000,FMT3)box_units,box_units,box_units
      close(7000)
    end subroutine debug_gro

    subroutine debug_template(filename)
      double precision :: box_units
      double precision :: x,y,z
      double precision :: xx(3,numatompertemp)
      character(len=*) :: filename
      character(len=*),parameter :: FM1 = "(i5)"
      character(len=*),parameter :: FMT2 = "( i5 , 2a5 , i5 , 3f8.3)"
      character(len=*),parameter :: FMT3 = "(f10.5,f10.5,f10.5)"
      integer :: i,flag
      integer :: typ,type
      character(len=3) :: char

      open(unit = 7000, file = trim(grofile_path)//'/'//filename)
      write(7000,*)'xyz'

      write(7000,FM1)numatompertemp
  
      do i = 1,numatompertemp
      
         if(i.eq.3)char='N'
         if(i.eq.2)char='CA'
         if(i.eq.1)char='C'
         if(i.eq.4)char='O'
         x  = BOPtemp_prune(i)%x/10.0d0
         y  = BOPtemp_prune(i)%y/10.0d0
         z  = BOPtemp_prune(i)%z/10.0d0
         
         write(7000,FMT2)i,'TTA',char,i, x,y,z
      enddo



      box_units = box/(10.0d0)
      
      write(7000,FMT3)box_units,box_units,box_units
      close(7000)
    end subroutine debug_template

    subroutine xyz_template_Gro(local_template,filename,flag)
      implicit none
      double precision :: local_template(3,150)
      double precision :: box_units
      double precision :: x,y,z
      character(len=*) :: filename
      character(len=*),parameter :: FM1 = "(i5)"
      character(len=*),parameter :: FMT2 = "( i5 , 2a5 , i5 , 3f8.3)"
      character(len=*),parameter :: FMT3 = "(f10.5,f10.5,f10.5)"
      integer :: i,flag
      integer :: typ,type
      character(len=3) :: char

      open(unit = 7000, file = trim(grofile_path)//'/'//filename)
      write(7000,*)'xyz'

      write(7000,FM1)2*numatompertemp
    ! if(flag.eq.1)then
    !     write(7000,FM1)xyz_number 
    !  else
    !     write(7000,FM1)xyz_number +150
    !  endif
      
   !   do i = 1,xyz_number
      do i = 1,numatompertemp
         typ  = position(i)%type
         char = type2string(typ)
         char = adjustl(char)
         char = trim(char)          
         
         
         x = xyz(1,i)/10.0d0
         y = xyz(2,i)/10.0d0
         z = xyz(3,i)/10.0d0
         
         write(7000,FMT2)i,'TTA',char,i, x,y,z
      enddo

      if(flag.eq.0)then
     !    do i = 1,150
         do i = 1,numatompertemp
            typ  = position(i)%type
            char = type2string(typ)
            char = adjustl(char)
            char = trim(char)          
            
            
            x = local_template(1,i)/10.0d0
            y = local_template(2,i)/10.0d0
            z = local_template(3,i)/10.0d0
            
            write(7000,FMT2)i,'TTA',char,i, x,y,z
         enddo
      endif


      box_units = box/(10.0d0)
      
      write(7000,FMT3)box_units,box_units,box_units
      close(7000)
    end subroutine xyz_template_Gro


    subroutine ptta_gro(filename)
      implicit none
      double precision :: local_template(3,150)
      double precision :: box_units
      double precision :: x,y,z
      character(len=*),parameter :: FM1 = "(i5)"
      character(len=*),parameter :: FMT2 = "( i5 , 2a5 , i5 , 3f8.3)"
      character(len=*),parameter :: FMT3 = "(f10.5,f10.5,f10.5)"
      character(len=*) :: filename
      integer :: i,j,spot
      integer :: typ
      character(len=3) :: char

      open(unit = 7000, file = trim(grofile_path)//'/'//filename)
      write(7000,*)'ptta'

      write(7000,FM1)numMolArray(1)*numatns

      do i = 1,numMolArray(1)
         do j = 1, numatns
                      
            if(j.eq.1)char='C'
            if(j.eq.2)char='CA'
            if(j.eq.3)char='N'
            if(j.eq.4)char='O'
            if(j.eq.5)char='O'
     
            
            spot = (i-1)*numatns+j
            x = ptta(spot)%x/10.0d0
            y = ptta(spot)%y/10.0d0
            z = ptta(spot)%z/10.0d0
            
            write(7000,FMT2)i,'TTA',char,i, x,y,z
         enddo
      enddo
      box_units = box/(10.0d0)
      
      write(7000,FMT3)box_units,box_units,box_units
      close(7000)
    end subroutine ptta_gro


    subroutine XYZ_2_GRO(filename1,filename2,length)
      implicit none
      double precision :: local_template(3,150)
      double precision :: box_units
      double precision :: x,y,z
      character(len=*),parameter :: FM1 = "(i5)"
      character(len=*),parameter :: FMT2 = "( i5 , 2a5 , i5 , 3f8.3)"
      character(len=*),parameter :: FMT3 = "(f10.5,f10.5,f10.5)"
      character(len=*) :: filename1, filename2
      integer :: i,j,spot,length
      integer :: typ
      character(len=3) :: char

      open(unit = 7000, file = trim(grofile_path)//'/'//filename1)
      open(unit = 7001, file = trim(grofile_path)//'/'//filename2)

      write(7001,*)'bulk crystal structure'

      write(7001,FM1)length

      do i = 1,length
            
         read(7000,*)x,y,z
         x = x/10.0d0; y = y/10.0d0; z = z/10.0d0

         typ  = position(i)%type
         char = type2string(typ)
         char = adjustl(char)
         char = trim(char)          
         
        
         
         write(7001,FMT2)i,'TTA',char,i, x,y,z
      enddo
      box_units = box/(10.0d0)
      
      write(7001,FMT3)box_units,box_units,box_units
      close(7000)
      close(7001)
    end subroutine XYZ_2_GRO




     subroutine Grofile_template
      implicit none
      double precision :: x,y,z
      double precision :: box_units
      integer :: length_clus,type,cluster,numAt
      character(len = 50) :: filename,numchar,char
      character(len=*),parameter :: FM1 = "(i5)"
      character(len=*),parameter :: FMT2 = "( i5 , 2a5 , i5 , 3f8.3)"
      character(len=*),parameter :: FMT3 = "(f10.5,f10.5,f10.5)"
      integer :: loop1,loop2
      integer :: typ,j,i
      integer :: spot,count
            
      do j = 1,numTemPlate

      
         if(j.eq.1)then
            open(unit = 1003, file = trim(grofile_path)//'/'//'alpha0_np.gro')
            open(unit = 1004, file = trim(grofile_path)//'/'//'alpha0.gro')

         elseif(j.eq.2)then
            open(unit = 1003, file = trim(grofile_path)//'/'//'alpha1_np.gro')
            open(unit = 1004, file = trim(grofile_path)//'/'//'alpha1.gro')

         elseif(j.eq.3)then
            open(unit = 1003, file = trim(grofile_path)//'/'//'alpha2_np.gro')
            open(unit = 1004, file = trim(grofile_path)//'/'//'alpha2.gro')

         elseif(j.eq.4)then
            open(unit = 1003, file = trim(grofile_path)//'/'//'alpha3_np.gro')
            open(unit = 1004, file = trim(grofile_path)//'/'//'alpha3.gro')

         elseif(j.eq.5)then
            open(unit = 1003, file = trim(grofile_path)//'/'//'beta0_np.gro')
            open(unit = 1004, file = trim(grofile_path)//'/'//'beta0.gro')

         elseif(j.eq.6)then
            open(unit = 1003, file = trim(grofile_path)//'/'//'beta1_np.gro')
            open(unit = 1004, file = trim(grofile_path)//'/'//'beta1.gro')

         end if




         write(1003,*)'Lennard Jones'
         write(1003,FM1)numatns*nummoltemp

         count = 0
         do loop1 = 1,nummoltemp            
            do loop2 = 1,numatns
               
               count = count + 1
               spot = (loop1-1)*10+loop2
          

               x = tempa_noprune(j)%temp(1,count)/10.0d0
               y = tempa_noprune(j)%temp(2,count)/10.0d0
               z = tempa_noprune(j)%temp(3,count)/10.0d0
               
               if(loop2.eq.3)char='N'
               if(loop2.eq.2)char='CA'
               if(loop2.eq.1)char='C'
               if(loop2.eq.4)char='O'
               if(loop2.eq.5)char='O'
      
              
               write(1003,FMT2)spot,'TT',char,spot, x,y,z
               
            enddo
         enddo

         write(1004,*)'Lennard Jones'
         write(1004,FM1)numatomtemp

         count = 0
         do loop1 = 1,nummoltemp            
            do loop2 = 1,numatompertemp
               
               count = count + 1
               spot = (loop1-1)*10+loop2
          

               x = tempa(j)%temp(1,count)/10.0d0
               y = tempa(j)%temp(2,count)/10.0d0
               z = tempa(j)%temp(3,count)/10.0d0
               
               if(loop2.eq.3)char='N'
               if(loop2.eq.2)char='CA'
               if(loop2.eq.1)char='C'
               if(loop2.eq.4)char='O'
               if(loop2.eq.5)char='O'
      
              
               write(1004,FMT2)spot,'TT',char,spot, x,y,z
               
            enddo
         enddo

         box_units = box/10.0d0

         write(1003,FMT3)box_units,box_units,box_units
         write(1004,FMT3)box_units,box_units,box_units

         close(1003)
      enddo
    end subroutine Grofile_template


    subroutine GroFile
      implicit none
      double precision :: x,y,z
      integer :: length_clus,type,cluster,numAt
      character(len=*),parameter :: FM1 = "(i5)"
      character(len=*),parameter :: FMT2 = "( i5 , 2a5 , i5 , 3f8.3)"
      character(len=*),parameter :: FMT3 = "(f10.5,f10.5,f10.5)"
      integer :: loop
      integer :: typ
      
      open(unit = 1003,file  = trim(grofile_path)//'/'//'LJ_Gro.dat')
      
      
      write(1003,*)'Lennard Jones'
      write(1003,FM1)np
      do loop = 1, np
         
         x = position(loop)%x
         y = position(loop)%y
         z = position(loop)%z
         write(1003,FMT2)loop,'SOlid','SO',loop, x,y,z
         
      enddo
      
      write(1003,FMT3)box,box,box
      close(1003)
      
    end subroutine GroFile


    
    subroutine GroFile_COM(filename,type)
      implicit none
      double precision :: x,y,z
      integer :: length_clus,type,cluster,numAt
      character(len=*),parameter :: FM1 = "(i5)"
      character(len=*),parameter :: FMT2 = "( i5 , 2a5 , i5 , 3f8.3)"
      character(len=*),parameter :: FMT3 = "(f10.5,f10.5,f10.5)"
      character(len=*) :: filename
      integer :: loop
      integer :: typ
      
      open(unit = 1003,file  = trim(grofile_path)//'/'//filename)
      
      
      write(1003,*)'Lennard Jones'
      write(1003,FM1)NumMolArray(type)
      do loop = 1, NumMolArray(type)
         
         x = com(loop)%x/10.0d0
         y = com(loop)%y/10.0d0
         z = com(loop)%z/10.0d0
         write(1003,FMT2)loop,'SOlid','SO',loop, x,y,z
         
      enddo
      
      write(1003,FMT3)box/10.0d0,box/10.0d0,box/10.0d0
      close(1003)
      
    end subroutine GroFile_COM

   subroutine mol_GroFile_cluster(filename,cluster)
      implicit none
      double precision :: x,y,z
      double precision :: box_units
      integer :: length_clus,type,cluster,numAt
      character(len=*),parameter :: FM1 = "(i5)"
      character(len=*),parameter :: FMT2 = "( i5 , 2a5 , i5 , 3f8.3)"
      character(len=*),parameter :: FMT3 = "(f10.5,f10.5,f10.5)"
      integer :: loop,j,k,ca,countmol
      integer :: typ
      character(len=3) :: char
      character(len=*) :: filename


      open(unit = 1003,file  = trim(grofile_path)//'/'//filename)
      ca = 0
      countmol = 0
      write(1003,*)'Molecular System'
      write(1003,FM1)np
      open(unit = 1004, file ='check.dat')
      do loop = 1,numMolType

         do j = 1,numMolArray(loop)
         
            countmol  = countmol+1
            do k = 1,numMolAtom(loop)
               ca = ca + 1
               typ  = position(ca)%type
               char = type2string(typ)
               char = adjustl(char)
               char = trim(char)              

               !---must convert to nm
               x = position(ca)%x/10.0d0
               y = position(ca)%y/10.0d0
               z = position(ca)%z/10.0d0

               if(j.ne.769)then
                  write(1003,FMT2)countmol,'GLY',char,ca, x,y,z
               else
                  write(1003,FMT2)countmol,'GlY','SOL',ca,x,y,z
               endif


              ! if(belongsto(j).ne.cluster)then
              !    write(1003,FMT2)countmol,'SOL','SOL',ca, x,y,z
              !    write(1004,*)j,denstrack(j),minrmsd(j)
              ! else
              !    write(1003,FMT2)countmol,'LIQ',char,ca, x,y,z
              ! endif
               
            enddo
         enddo
      enddo
      
      box_units = box/(10.0d0)
      write(1003,FMT3)box_units,box_units,box_units
      close(1003)
      close(1004)
      
    end subroutine Mol_GroFile_cluster

    subroutine mol_GroFile(filename)
      implicit none
      double precision :: x,y,z
      double precision :: box_units
      integer :: length_clus,type,cluster,numAt
      character(len=*),parameter :: FM1 = "(i5)"
      character(len=*),parameter :: FMT2 = "( i5 , 2a5 , i5 , 3f8.3)"
      character(len=*),parameter :: FMT3 = "(f10.5,f10.5,f10.5)"
      integer :: loop,j,k,ca,countmol
      integer :: typ
      character(len=3) :: char
      character(len=*) :: filename


      open(unit = 1003,file  = trim(grofile_path)//'/'//filename)
      ca = 0
      countmol = 0
      write(1003,*)'Molecular System'
      write(1003,FM1)np
      do loop = 1,numMolType

         do j = 1,numMolArray(loop)
         
            countmol  = countmol+1
            do k = 1,numMolAtom(loop)
               ca = ca + 1
               typ  = position(ca)%type
               char = type2string(typ)
               char = adjustl(char)
               char = trim(char)              

               !---must convert to nm
               x = position(ca)%x/10.0d0
               y = position(ca)%y/10.0d0
               z = position(ca)%z/10.0d0



               if(loop.eq.1)then
                  write(1003,FMT2)countmol,'GLY',char,ca, x,y,z
               else
                  write(1003,FMT2)countmol,'CTC',char,ca, x,y,z
               endif
            enddo
         enddo
      enddo
      
      box_units = box/(10.0d0)
      write(1003,FMT3)box_units,box_units,box_units
      close(1003)
      
    end subroutine Mol_GroFile


    subroutine timeseries_GroFile
      implicit none
      double precision :: x,y,z
      integer :: length_clus,type,cluster,numAt
      character(len=*),parameter :: FM1 = "(i5)"
      character(len=*),parameter :: FMT2 = "( i5 , 2a5 , i5 , 3f8.3)"
      character(len=*),parameter :: FMT3 = "(f10.5,f10.5,f10.5)"
      integer :: loop,typ,id,j,jflag
      logical :: exist
      character(len=3) :: char
      integer :: T1,T2,clock_rate,clock_max
      
      
      
      call system_clock(T1,clock_rate,clock_max)

      inquire(file=trim(grofile_path)//'/'//'time_series.gro',exist=exist)
      if(exist)then
         open(unit = 99999, file = trim(grofile_path)//'/'//'time_series.gro',status='old',position='append',action='write')
      else
         open(unit = 99999, file = trim(grofile_path)//'/'//'time_series.gro',status='new',action='write')
      endif

      write(99999,*)'Lennard Jones'
      write(99999,FM1)np
      do loop = 1, np
         

         do j = 1,np
            id = global_id(j)
            if(id.eq.loop)then
               jflag = j
               exit
            endif
         enddo

         typ = position(jflag)%type
         char = type2string(typ)
         char = adjustl(char)
         char = trim(char)

         typ  = ltype2ftype(typ)

         x = position(jflag)%x/10.0d0; y = position(jflag)%y/10.0d0; z = position(jflag)%z/10.0d0

         write(99999,FMT2)loop,'TTA',char,loop,x,y,z
         
      enddo
      
      write(99999,FMT3)box/10.d0,box/10.d0,box/10.d0
      call system_clock(T2,clock_rate,clock_max)
      
      time_video = time_video + real(T2-T1)/real(clock_rate)

    end subroutine timeseries_GroFile





    
    subroutine timeseries_ptta_GroFile
      implicit none
      double precision :: x,y,z
      integer :: length_clus,type,cluster,numAt
      character(len=*),parameter :: FM1 = "(i5)"
      character(len=*),parameter :: FMT2 = "( i5 , 2a5 , i5 , 3f8.3)"
      character(len=*),parameter :: FMT3 = "(f10.5,f10.5,f10.5)"
      integer :: loop,typ,id,j,jflag,i
      character(len=3) :: char
      integer :: T1,T2,clock_rate,clock_max
      logical :: exist
      integer :: atom
      call system_clock(T1,clock_rate,clock_max)

  
      inquire(file=trim(grofile_path)//'/'//'time_series.gro',exist=exist)
      if(exist)then
         open(unit = 99999, file = trim(grofile_path)//'/'//'time_series.gro',status='old',position='append',action='write')
      else
         open(unit = 99999, file = trim(grofile_path)//'/'//'time_series.gro',status='new',action='write')
      endif


      atom = nummolatom(1)
      do i = 1,np

         id = global_id(i)

         if(id.lt.nummolarray(1)*nummolatom(1))then
        
            x = position(i)%x
            y = position(i)%y
            z = position(i)%z

            do while (x > box) 
               x = x - box 
            enddo
            do while(x < 0.0)                     
               x = x + box     
            enddo
            
            
            do while (y > box) 
               y = y - box
            enddo
            do while (y < 0.0) 
               y = y + box     
            enddo
            
            do while (z > box) 
               z = z - box 
            enddo
            do while  (z < 0.0) 
               z = z + box     
            enddo


            ptta_grof(id)%x = x
            ptta_grof(id)%y = y
            ptta_grof(id)%z = z


         endif
      enddo
      
      
      
      
      write(99999,*)'Lennard Jones'
      write(99999,FM1)nummolarray(1)*nummolatom(1)
      do loop = 1, nummolarray(1)
         do j = 1, atom
            
            if(j.eq.1)char='C'
            if(j.eq.2)char='CA'
            if(j.eq.3)char='N'
            if(j.eq.4)char='O'
            if(j.eq.5)char='O'
            if(j.eq.6)char='H'
            if(j.eq.7)char='H'
            if(j.eq.8)char='H'
            if(j.eq.9)char='H'
            if(j.eq.10)char='H'

     
            jflag = (loop-1)*atom + j
            
            x = ptta_grof(jflag)%x; y = ptta_grof(jflag)%y; z = ptta_grof(jflag)%z
            
            x = x/10.0d0; y = y/10.0d0; z = z/10.0d0
            write(99999,FMT2)loop,'TTA',char,atom,x,y,z
         enddo
      enddo
      
      write(99999,FMT3)box/10.0d0,box/10.0d0,box/10.0d0
      close(99999)
      call system_clock(T2,clock_rate,clock_max)
      
      time_video = time_video + real(T2-T1)/real(clock_rate)

    end subroutine timeseries_ptta_GroFile


   subroutine timeseries_ptta_GroFile_seeded
      implicit none
      double precision :: x,y,z
      integer :: length_clus,type,cluster,numAt
      character(len=*),parameter :: FM1 = "(i5)"
      character(len=*),parameter :: FMT2 = "( i5 , 2a5 , i5 , 3f8.3)"
      character(len=*),parameter :: FMT3 = "(f10.5,f10.5,f10.5)"
      integer :: loop,typ,id,j,jflag,i
      character(len=3) :: char
      integer :: T1,T2,clock_rate,clock_max
      logical :: exist


      call system_clock(T1,clock_rate,clock_max)



      inquire(file=trim(grofile_path)//'/'//trim(sim_string)//'/'//'cluster_num.dat',exist=exist)
      if(exist)then
         open(unit = 99999, file = trim(grofile_path)//'/'//trim(sim_string)//'/'//'cluster_num.dat',status='old',position='append',action='write')
      else
         open(unit = 99999, file = trim(grofile_path)//'/'//trim(sim_string)//'/'//'cluster_num.dat',status='new',action='write')
      endif
   

      do i = 1,np
         if(mol(1,i).eq.1)then
            x = ptta(mol(2,i))%x
            y = ptta(mol(2,i))%y
            z = ptta(mol(2,i))%z

            if (x > box) then
               x = x - box 
            elseif  (x < 0.d0) then                    
               x = x + box     
            endif
            
            
            if (y > box) then
               y = y - box 
            elseif  (y < 0.d0) then
              y = y + box     
            endif
            
            if (z > box) then
               z = z - box 
            elseif  (z < 0.d0) then
               z = z + box     
            endif


            ptta_grof(mol(2,i))%x = x
            ptta_grof(mol(2,i))%y = y
            ptta_grof(mol(2,i))%z = z


         endif
      enddo
      
      
      
      
      write(201,*)'Lennard Jones'
      write(201,FM1)nummolarray(1)*numatns
      do loop = 1, nummolarray(1)
         do j = 1,numatns
            
            if(j.eq.1)char='C'
            if(j.eq.2)char='CA'
            if(j.eq.3)char='N'
            if(j.eq.4)char='O'
            if(j.eq.5)char='O'
     
            jflag = (loop-1)*numatns + j
            
            x = ptta_grof(jflag)%x; y = ptta_grof(jflag)%y; z = ptta_grof(jflag)%z
            
            x = x/10.0d0; y = y/10.0d0; z = z/10.0d0
            write(99999,FMT2)loop,'TTA',char,numatns,x,y,z
         enddo
      enddo
      
      write(99999,FMT3)box/10.0d0,box/10.0d0,box/10.0d0
      
      call system_clock(T2,clock_rate,clock_max)
      
      time_video = time_video + real(T2-T1)/real(clock_rate)

    end subroutine timeseries_ptta_GroFile_seeded











    subroutine timeseries_ptta_GroFile_clus(greatest)
      implicit none
      double precision :: x,y,z
      integer :: length_clus,type,cluster,numAt
      character(len=*),parameter :: FM1 = "(i5)"
      character(len=*),parameter :: FMT2 = "( i5 , 2a5 , i5 , 3f8.3)"
      character(len=*),parameter :: FMT3 = "(f10.5,f10.5,f10.5)"
      integer :: loop,typ,id,j,jflag,i
      integer :: greatest
      character(len=3) :: char


      do i = 1,np
         if(mol(1,i).eq.1)then
            x = ptta(mol(2,i))%x
            y = ptta(mol(2,i))%y
            z = ptta(mol(2,i))%z

            if (x > box) then
               x = x - box 
            elseif  (x < 0.d0) then                    
               x = x + box     
            endif
            
            
            if (y > box) then
               y = y - box 
            elseif  (y < 0.d0) then
              y = y + box     
            endif
            
            if (z > box) then
               z = z - box 
            elseif  (z < 0.d0) then
               z = z + box     
            endif


            ptta_grof(mol(2,i))%x = x
            ptta_grof(mol(2,i))%y = y
            ptta_grof(mol(2,i))%z = z


         endif
      enddo
      
      
      
      
      write(202,*)'Lennard Jones'
      write(202,FM1)nummolarray(1)*numatns
      do loop = 1, nummolarray(1)


         if(belongsto(loop).ne.greatest)then
            do j = 1,numatns
               
               if(j.eq.1)char='C'
               if(j.eq.2)char='CA'
               if(j.eq.3)char='N'
               if(j.eq.4)char='O'
               if(j.eq.5)char='O'
               
               jflag = (loop-1)*numatns + j
               
               x = ptta_grof(jflag)%x/10.0d0; y = ptta_grof(jflag)%y/10.0d0; z = ptta_grof(jflag)%z/10.0d0
               
               write(202,FMT2)loop,'TTA',char,numatns,x,y,z
            enddo
            
         else
             do j = 1,numatns
               
               if(j.eq.1)char='SO'
               if(j.eq.2)char='SO'
               if(j.eq.3)char='SO'
               if(j.eq.4)char='SO'
               if(j.eq.5)char='SO'
               
               jflag = (loop-1)*numatns + j
               
               x = ptta_grof(jflag)%x/10.0d0; y = ptta_grof(jflag)%y/10.0d0; z = ptta_grof(jflag)%z/10.0d0
               
               write(202,FMT2)loop,'TTA',char,numatns,x,y,z
            enddo
         endif
      end do



         
      write(202,FMT3)box/10.d0,box/10.d0,box/10.d0
      
    end subroutine timeseries_ptta_GroFile_clus












    subroutine XYZfile
      implicit none
      integer :: i

      write(200,*)np
      write(200,*)
      do i = 1,np
         write(200,*)'LJ',position(i)%x,position(i)%y,position(i)%z
      enddo
    end subroutine XYZfile

    subroutine tinker_XYZfile(filename)
      implicit none
      double precision :: x,y,z
      integer :: i,typ,J
      integer :: a1,a2,a3
      character(len=3) :: char
      character(len=*),parameter :: FM1 = "(i6,a5,f12.6,2f12.6,5i6)"
      character(len=*),parameter :: FM2 = "(i6)"
      character(len=*) :: filename


      open(unit = 500,file = trim(grofile_path)//'/'//filename )

      write(500,FM2)np

      do i = 1,np
         typ  = position(i)%type
         
         x = position(i)%x
         y = position(i)%y
         z = position(i)%z
         char = type2string(typ)
         char = adjustl(char)
         char = trim(char)

         typ  = ltype2ftype(typ)




         if(num1bond(i).eq.1)then
            write(500,FM1)i,char,x,y,z,typ,specbond(1,i)

         elseif(num1bond(i).eq.2)then
           write(500,FM1)i,char,x,y,z,typ,specbond(1,i),specbond(2,i)

         elseif(num1bond(i).eq.3)then
            write(500,FM1)i,char,x,y,z,typ,specbond(1,i),specbond(2,i),specbond(3,i)

         elseif(num1bond(i).eq.4)then
            write(500,FM1)i,char,x,y,z,typ,specbond(1,i),specbond(2,i),specbond(3,i),specbond(4,i)

         elseif(num1bond(i).eq.0)then
            write(500,FM1)i,char,x,y,z,typ
         endif
      end do
      print*,'we have left tinker:',i,np
      close(500)
    end subroutine tinker_XYZfile


    subroutine XYZ_2_TinkXYZfile(filename1,filename2,length)
      implicit none
      double precision :: x,y,z
      integer :: i,typ,J,length
      integer :: a1,a2,a3
      character(len=3) :: char
      character(len=*),parameter :: FM1 = "(i6,a5,f12.6,2f12.6,5i6)"
      character(len=*),parameter :: FM2 = "(i6)"
      character(len=*) :: filename1,filename2


      open(unit = 500,file = trim(grofile_path)//'/'//filename1)
      open(unit = 501,file = trim(grofile_path)//'/'//filename2)

      write(501,FM2)length

      do i = 1,length
         typ  = position(i)%type
         
         read(500,*)x,y,z


         char = type2string(typ)
         char = adjustl(char)
         char = trim(char)

         typ  = ltype2ftype(typ)


         if(num1bond(i).eq.1)then
            write(501,FM1)i,char,x,y,z,typ,specbond(1,i)

         elseif(num1bond(i).eq.2)then
           write(501,FM1)i,char,x,y,z,typ,specbond(1,i),specbond(2,i)

         elseif(num1bond(i).eq.3)then
            write(501,FM1)i,char,x,y,z,typ,specbond(1,i),specbond(2,i),specbond(3,i)

         elseif(num1bond(i).eq.4)then
            write(501,FM1)i,char,x,y,z,typ,specbond(1,i),specbond(2,i),specbond(3,i),specbond(4,i)

         elseif(num1bond(i).eq.0)then
            write(501,FM1)i,char,x,y,z,typ
         endif
      end do
      close(500)
      close(501)
    end subroutine XYZ_2_TinkXYZfile






    subroutine read_in_tinker(filename,flag,flag2)
      implicit none
      double precision :: x,y,z
      double precision :: dummy
      integer :: i,typ,j,k,l,localnp
      integer :: a1,a2,a3
      integer :: flag,flag2
      character(len=3) :: char
      character(len=*) :: filename

      open(unit = 3500, file = trim(grofile_path)//'/'//filename)

      
      if(flag2.eq.0)then
         read(3500,*)np
         
         if(flag.eq.1)then
            read(3500,*)dummy,dummy,dummy,dummy,dummy,dummy
         endif
         
         
         
         do i = 1 ,np
            
            if(num1bond(i).eq.1)then
               read(3500,*)l,char,x,y,z
               
               position(i)%x = x; position(i)%y = y; position(i)%z = z
               
               
            elseif(num1bond(i).eq.2)then
               read(3500,*)l,char,x,y,z
               
               position(i)%x = x; position(i)%y = y; position(i)%z = z
               
            elseif(num1bond(i).eq.3)then
               read(3500,*)l,char,x,y,z
               
               position(i)%x = x; position(i)%y = y; position(i)%z = z
               
            elseif(num1bond(i).eq.4)then
               read(3500,*)l,char,x,y,z
               
               position(i)%x = x; position(i)%y = y; position(i)%z = z
               
            elseif(num1bond(i).eq.0)then
               read(3500,*)l,char,x,y,z
               
               position(i)%x = x; position(i)%y = y; position(i)%z = z
               
            endif
         end do
         close(3500)
         
         

      elseif(flag2.eq.1)then
         read(3500,*)np
        
     
         numMolArray(2) = (np-NumMolArray(1)*NumMolAtom(1))/numMolAtom(2)
         print*,'NUMBER OF SOLVENT MOLECULES',numMolArray(2),np,NumMolArray(1)

         close(3500)
      end if
         
       
      
    end subroutine read_in_tinker

    subroutine tinker_KEYfile(flag)
      implicit none
      !character(len=*),parameter :: FM1 = "(i6,a5,f12.6,2f12.6,5i6)"

      character(len=*),parameter :: FM1 = "(a10,a21)"
      character(len=*),parameter :: FM2 = "(a5,i10,i5,i5,2f10.2)"
      character(len=*),parameter :: FM3 = "(a7,i8,i5,i5,i5,f11.3,f4.1,i2,f8.3,f6.1,i2,f8.3,f6.1,i2)"
      character(len=*),parameter :: FM4 = "(a6,f21.2)"
      character(len=*),parameter :: FM5 = "(a10,f16.2)"
      character(len=*),parameter :: FM6 = "(a12,f21.2)"
      character(len=*),parameter :: FM7 = "(a12,3i6)"
      character(len=*),parameter :: FM8 = "(a20,a20)"
      character(len=*),parameter :: FM9 = "(a13)"


      character(len=*),parameter :: FM20 = "(a5)"
      character(len=*),parameter :: FM21 = "(a8,f21.2)"
      character(len=*),parameter :: FM22 = "(a8,3i6)"
      character(len=*),parameter :: FM23 = "(a11,a12)"

      integer :: i,flag
      open(unit = 3500,file =trim(grofile_path)//'/'//'tink.key')
      
      if(flag.eq.3.or.flag.eq.4)then
         write(3500,FM1)'parameters','oplsaa_actual.prm'
      else
         write(3500,FM1)'parameters','oplsaa_nocharge.prm'
      endif


      write(3500,FM2)'angle', 5,3, 19, 70.00, 108.00
      write(3500,FM2)'angle', 4,3, 19, 80.00, 120.40
      write(3500,FM2)'angle', 3,19,19, 150.00,180.00

      write(3500,FM3)'torsion', 4,3,19,19, 0.000,0.0,1, 0.000,180.0,2, 0.000,0.0,3
      write(3500,FM3)'torsion', 5,3,19,19, 0.000,0.0,1, 0.000,180.0,2, 0.000,0.0,3
      write(3500,FM3)'torsion', 19,3,5,7,  1.500,0.0,1, 5.500,180.0,2, 0.000,0.0,3

      do i = 1,3
         write(3500,*)
      end do
      write(3500,FM9)'neighbor-list'


      do i = 1,3
         write(3500,*)
      end do

      

      write(3500,FM4)'a-axis',box
      write(3500,FM5)'vdw-cutoff',rcut


      if(flag.eq.3.or.flag.eq.4)then
         
         do i = 1,3
            write(3500,*)
         enddo
         write(3500,FM20)'ewald'
         write(3500,FM6)'ewald-cutoff',12.0
         !write(3500,FM22)'pme-grid',64,64,64
         !write(3500,FM23)'fft-package','FFTW'
      endif
      
      close(3500)

      print*,'DONE MAKING TINK KEY'
    end subroutine tinker_KEYfile

    subroutine lets_convert
      implicit none
      character(len=150) :: app
      character(len=150) :: files
      integer :: i,num
      
      do i = 1,400

         if(i.lt.10)then
            files = 'tta.00'
         elseif(i.ge.10.and.i.lt.100)then
            files = 'tta.0'
         else
            files = 'tta.'
         endif
         num = i
         write(app,'(i8)')num
         app = ADJUSTL(app)
         app = trim(app)
         
         files = trim(files)//trim(app)

         print*,'files',files
         call read_in_tinker(files,1,0)
      !   call pbc(0)
         call timeseries_GroFile
      enddo
    end subroutine lets_convert

    subroutine clear
        implicit none
      character(len=150) :: app
      character(len=150) :: files
      integer :: i,num
      
     
      call system('rm coords.001')
      call system('rm coords.002')
      call system('rm coords.003')
      call system('rm coords.004')
      call system('rm coords.005')
      call system('rm coords.006')
      call system('rm coords.007')
      call system('rm coords.008')
      call system('rm coords.009')
      call system('rm coords.010')
      call system('rm coords.011')
      call system('rm coords.012')
      call system('rm coords.013')
      call system('rm coords.014')
      call system('rm coords.015')
      call system('rm coords.016')
      call system('rm coords.017')
      call system('rm coords.018')
      call system('rm coords.019')
      call system('rm coords.020')
    end subroutine clear


    subroutine read_in_lammps(filename)
      implicit none
      integer :: i,j,junk,id
      double precision :: x,y,z
      character(len=*) :: filename
      
      open(unit = 500, file = filename)
      !read(500,*)np
      !read(500,*)

      do i = 1,np
         read(500,*)id,x,y,z
        
         position(id)%x = x
         position(id)%y = y
         position(id)%z = z
         
      end do
      close(500)
    end subroutine read_in_lammps

     subroutine read_in_lammps_beta(filename)
      implicit none
      integer :: i,j,junk,count
      double precision :: x,y,z
      double precision :: xx(10),yy(10),zz(10)
      character(len=*) :: filename
      
      open(unit = 500, file = filename)
      !read(500,*)np
      !read(500,*)

      count = 0

      print*,'WE ARE IN LAMMPS BETA'
      do i = 1,nummolarray(1)

         do j = 1,nummolatom(1)
            read(500,*)x,y,z
            xx(j) = x
            yy(j) = y
            zz(j) = z
         enddo
         do j = 1,nummolatom(1)
            junk = (i-1)*nummolatom(1)
            
            if(j.eq.1)then
               position(junk+1)%x = xx(1)
               position(junk+1)%y = yy(1)
               position(junk+1)%z = zz(1)
            endif
            
            if(j.eq.2)then
               position(junk+2)%x = xx(2)
               position(junk+2)%y = yy(2)
               position(junk+2)%z = zz(2)
            endif

            if(j.eq.3)then
               position(junk+6)%x = xx(3)
               position(junk+6)%y = yy(3)
               position(junk+6)%z = zz(3)
            endif
     
            if(j.eq.4)then
               position(junk+7)%x = xx(4)
               position(junk+7)%y = yy(4)
               position(junk+7)%z = zz(4)
            endif

            if(j.eq.5)then
               position(junk+8)%x = xx(5)
               position(junk+8)%y = yy(5)
               position(junk+8)%z = zz(5)
            endif
            
            if(j.eq.6)then
               position(junk+9)%x = xx(6)
               position(junk+9)%y = yy(6)
               position(junk+9)%z = zz(6)
            endif

            if(j.eq.7)then
               position(junk+10)%x = xx(7)
               position(junk+10)%y = yy(7)
               position(junk+10)%z = zz(7)
            endif
            
            if(j.eq.8)then
               position(junk+3)%x = xx(8)
               position(junk+3)%y = yy(8)
               position(junk+3)%z = zz(8)
            endif
            
            if(j.eq.9)then
               position(junk+4)%x = xx(9)
               position(junk+4)%y = yy(9)
               position(junk+4)%z = zz(9)
            endif

            if(j.eq.10)then
               position(junk+5)%x = xx(10)
               position(junk+5)%y = yy(10)
               position(junk+5)%z = zz(10)
            endif
         enddo
      enddo
      close(500)
    end subroutine read_in_lammps_beta


    subroutine read_in_lammps_velocity(filename)
      implicit none
      integer :: i,j,junk,id
      double precision :: x,y,z,vx,vy,vz
      character(len=*) :: filename
      
      open(unit = 500, file = filename)
      !read(500,*)np
      !read(500,*)

      do i = 1,np
         read(500,*)id,x,y,z,vx,vy,vz
        
         position(id)%x = x
         position(id)%y = y
         position(id)%z = z
         v(id)%x = vx
         v(id)%y = vy
         v(id)%z = vz

      end do
      close(500)
    end subroutine read_in_lammps_velocity


    subroutine shift_lammps
      implicit none
      integer :: i
      do i = 1,np
         position(i)%x = position(i)%x + 0.50d0*box
         position(i)%y = position(i)%y + 0.50d0*box
         position(i)%z = position(i)%z + 0.50d0*box
      enddo
    end subroutine shift_lammps


    subroutine center_coords
      integer :: i,j,count
      double precision :: x,y,z
      double precision :: COMx,COMy,COMz,boxm
      

      count = 0
 
      comx = 0.0d0; comy =0.0d0; comz=0.0d0
      do i = 1,np


         comx = comx + position(i)%x
         comy = comy + position(i)%y
         comz = comz + position(i)%z
         
      end do
      comx  = comx/np; comy=comy/np; comz=comz/np
      boxm = box/2.0d0
      do i=1,np
         position(i)%x = position(i)%x + (boxm-COMx)
         position(i)%y = position(i)%y + (boxm-COMy)
         position(i)%z = position(i)%z + (boxm-COMz)
      enddo
    end subroutine center_coords



    subroutine make_grofile_sol(length,filename)
      implicit none
      integer :: length,type
      double precision :: x,y,z
      character(len=50) :: rank_char
      character(len=*):: filename
      character(len=*),parameter :: FM1 = "(i5)"
      character(len=*),parameter :: FMT2 = "( i5 , 2a5 , i5 , 3f8.3)"
      character(len=*),parameter :: FMT3 = "(f10.5,f10.5,f10.5)"
      character(len=50):: solid,clus,char
      integer :: loop,typ


      open(unit = 7896, file = trim(grofile_path)//'/'//filename)
      write(7896,*)'Lennard Jones'
      write(7896,FM1)length
      do loop = 1, length
         
         typ  = position(loop)%type
         char = type2string(typ)
         char = adjustl(char)
         char = trim(char)  
         x = position(loop)%x; y = position(loop)%y; z = position(loop)%z
         
         if(belongsto(loop).ne.0)then
            
            write(7896,FMT2)loop,'SOli','SO',loop,x,y,z

         elseif(belongsto(loop).eq.0)then
            write(7896,FMT2)loop,'LIqui','LI',loop, x,y,z
         endif
         
      enddo
      
      write(1,FMT3)box,box,box
      close(1)
      
    end subroutine make_grofile_sol
    
    
    subroutine make_grofile_cluster(cluster,length_clus,filename)
      implicit none
      double precision :: x,y,z
      integer :: length_clus,type,cluster,numAt
      character(len=50) :: rank_char,char
      character(len=*):: filename
      character(len=*),parameter :: FM1 = "(i5)"
      character(len=*),parameter :: FMT2 = "( i5 , 2a5 , i5 , 3f8.3)"
      character(len=*),parameter :: FMT3 = "(f10.5,f10.5,f10.5)"
      integer :: loop,typ,j,offset,i,atom

      atom = nummolatom(1)
      open(unit = 1,file  = trim(grofile_path)//'/'//filename)
      
      write(1,*)'Lennard Jones'
      write(1,FM1)length_clus*nummolatom(1)
      do loop = 1, nummolarray(1)
         
         if(belongsto(loop).eq.cluster)then

            do j =1 ,atom
               offset =(loop-1)*atom + j

                  
               if(j.eq.1)char='C'
               if(j.eq.2)char='Ca'
               if(j.eq.3)char='N'
               if(j.eq.4)char='O'
               if(j.eq.5)char='O'
               if(j.eq.6)char='H'
               if(j.eq.7)char='H'
               if(j.eq.8)char='H'
               if(j.eq.9)char='H'
               if(j.eq.10)char='H'
               
            
               x = ptta_full(offset)%x; y = ptta_full(offset)%y; z = ptta_full(offset)%z
               x = x/10.0d0; y = y/10.0d0; z = z/10.0d0
               write(1,FMT2)loop,'SOlid',char,loop,x,y,z
            enddo
         endif
      enddo
      
      write(1,FMT3)box/10.0d0,box/10.0d0,box/10.0d0
      close(1)
      
    end subroutine make_grofile_cluster


        
    subroutine make_grofile_cluster_step(step,cluster,length_clus)
      implicit none
      double precision :: x,y,z
      integer :: length_clus,type,cluster,numAt
      character(len=50) :: rank_char,char
      character(len=1500):: filename
      character(len=*),parameter :: FM1 = "(i5)"
      character(len=*),parameter :: FMT2 = "( i5 , 2a5 , i5 , 3f8.3)"
      character(len=*),parameter :: FMT3 = "(f10.5,f10.5,f10.5)"
      integer :: loop,typ,j,offset,i,step,atom

      atom = nummolatom(1)
      write(filename,'(i8)')step
      filename = adjustl(filename)
      filename= trim(grofile_path)//'/'//trim(filename)//'.gro'
      open(unit = 99999,file  = trim(filename))
      
      write(99999,*)'Lennard Jones'
      write(99999,FM1)length_clus*atom

     
      do loop = 1, nummolarray(1)
         
         if(belongsto(loop).eq.cluster)then

            do j =1 ,atom
               offset =(loop-1)*atom + j

                  
               if(j.eq.1)char='C'
               if(j.eq.2)char='Ca'
               if(j.eq.3)char='N'
               if(j.eq.4)char='O'
               if(j.eq.5)char='H'
               if(j.eq.6)char='H'
               if(j.eq.7)char='H'        
               if(j.eq.8)char='H'
               if(j.eq.9)char='H'
               if(j.eq.10)char='H'
               x = ptta_full(offset)%x; y = ptta_full(offset)%y; z = ptta_full(offset)%z
               x = x/10.0d0; y = y/10.0d0; z = z/10.0d0
               write(99999,FMT2)loop,'SOlid',char,loop,x,y,z
            enddo
         endif
      enddo
      
      write(99999,FMT3)box/10.0d0,box/10.0d0,box/10.0d0
      close(99999)
      
    end subroutine make_grofile_cluster_step
    
    

   subroutine make_grofile_cluster_mol(cluster,length_clus,filename)
      implicit none
      double precision :: x,y,z
      real*4  :: local_array(3,length_clus)
      integer :: length_clus,type,cluster
      integer :: i,j
      integer :: count
      character(len=50) :: rank_char,char
      character(len=*):: filename
      character(len=*),parameter :: FM1 = "(i5)"
      character(len=*),parameter :: FMT2 = "( i5 , 2a5 , i5 , 3f8.3)"
      character(len=*),parameter :: FMT3 = "(f10.5,f10.5,f10.5)"
      integer :: loop,typ


      count = 0
      do loop = 1, nummolarray(1)
         
         if(belongsto(loop).eq.cluster)then
            
            !print*,'found molecule',loop,denstrack(loop),minrmsd(loop)
            do i = 1,numatns

               count = count + 1
               j = (loop-1)*numatns +i
               x = ptta(j)%x; y = ptta(j)%y; z = ptta(j)%z


               do while (x.gt.box)
                  x = x - box 
               end do
               do while  (x .lt.0.d0)               
                  x = x + box   
               end do
               
               
               do while (y.gt.box) 
                  y = y - box
               enddo
               do while  (y.lt.0.d0) 
                  y = y + box    
               enddo
               
               
               do while (z.gt.box) 
                  z = z - box 
               enddo
               do while  (z.lt.0.d0) 
                  z = z + box    
               enddo
               
               
               local_array(1,count) = x
               local_array(2,count) = y
               local_array(3,count) = z
               
               
            enddo
         endif
         
      enddo



      count = 0
      open(unit = 1,file  = filename)
      write(1,*)'Lennard Jones'
      write(1,FM1)length_clus
      do loop = 1, nummolarray(1)
         
         
         if(belongsto(loop).eq.cluster)then
            
            do i = 1,numatns
               
               count = count + 1

               if(denstrack(loop).ge.4)then
                  if(i.eq.1)char='C'
                  if(i.eq.2)char='CA'
                  if(i.eq.3)char='N'
                  if(i.eq.4)char='O'
                  if(i.eq.5)char='O'
               else
                  char ='L'
               endif

               j = (loop-1)*numatns+i
      
               char = trim(char)  
               char = adjustl(char)

               x = local_array(1,count); y = local_array(2,count); z = local_array(3,count)
               x = x/10.0d0; y = y/10.0d0; z = z/10.0d0

               write(1,FMT2)loop,'GLY',char,loop,x,y,z
               
               !if(denstrack(loop).ge.densmin)then
               !   write(1,FMT2)loop,'SO','SO',loop,x,y,z
               !else
               !   write(1,FMT2)loop,'SO','LI',loop,x,y,z
               !endif
            enddo
         endif
         
      enddo
      
      write(1,FMT3)box/10.0d0,box/10.0d0,box/10.0d0
      close(1)
      
    end subroutine make_grofile_cluster_mol





    subroutine xyzf_template_gro(xyzf,local_template,filename)
      implicit none
      double precision :: x,y,z,box_units
      double precision :: local_template(3,150),xyzf(3,150)
      character(len=*),parameter :: FM1 = "(i5)"
      character(len=*),parameter :: FMT2 = "( i5 , 2a5 , i5 , 3f8.3)"
      character(len=*),parameter :: FMT3 = "(f10.5,f10.5,f10.5)"
      character(len=3) :: char2
      character(len=*) :: filename
      integer :: i,flag,spot,j,spot1
      integer :: typ,type,count
      integer :: numolgro
      character(len=3) :: char

      open(unit = 7000, file = trim(grofile_path)//'/'//filename)

      write(7000,*)'xyz'

      write(7000,FM1)2*uniquepairs*numatompertemp
 
      count = 0
      do i= 1,uniquepairs
         do j = 1,numatompertemp
            spot = (i-1)*10+j
            typ  = position(spot)%type
            char = type2string(typ)
            char = adjustl(char)
            char = trim(char)   

            count = count + 1

            
            if(j.eq.1)char='C'
            if(j.eq.2)char='Ca'
            if(j.eq.3)char='N'
            if(j.eq.4)char='O'
            if(j.eq.5)char='O'
            

            x    = xyzf(1,count)/10.0d0
            y    = xyzf(2,count)/10.0d0
            z    = xyzf(3,count)/10.0d0
            write(7000,FMT2)i,'TT',char,i,x,y,z

         enddo
      enddo

      count = 0
      do i= 1,uniquepairs
         do j = 1,numatompertemp
          
            
          ! if(j.eq.1)char='C'
          ! if(j.eq.2)char='Ca'
          ! if(j.eq.3)char='N'
          ! if(j.eq.4)char='O'

            char ='T'
            count = count + 1
            x    = tempf(1,count)/10.0d0
            y    = tempf(2,count)/10.0d0
            z    = tempf(3,count)/10.0d0
            write(7000,FMT2)i,'TT',char,i,x,y,z

         enddo
      enddo

      !do i =1,numatomtemp
      !   write(7000,FMT2)numatomtemp+i,'LJ','ljsol',numatomtemp+i, local_template(1,i),local_template(2,i),local_template(3,i)
      !enddo
   
      box_units = box/(10.0d0)
      
      write(7000,FMT3)box_units,box_units,box_units
      close(7000)

    end subroutine xyzf_template_gro


    subroutine xyz_template_centroid(local_template,len,filename)
      implicit none
      double precision :: local_template(3,len)
      double precision :: box_units
      double precision :: x,y,z
      character(len=*) :: filename
      character(len=*),parameter :: FM1 = "(i5)"
      character(len=*),parameter :: FMT2 = "( i5 , 2a5 , i5 , 3f8.3)"
      character(len=*),parameter :: FMT3 = "(f10.5,f10.5,f10.5)"
      integer :: i,flag
      integer :: typ,type,len
      character(len=3) :: char
      
      open(unit = 7000, file = trim(grofile_path)//'/'//filename)
      write(7000,*)'xyz'
      
      write(7000,FM1)2*numatompertemp
      
      
      do i = 1,numatompertemp
         typ  = position(i)%type
         char = type2string(typ)
         char = adjustl(char)
         char = trim(char)          
         
         
         if(i.eq.1)char='C'
         if(i.eq.2)char='Ca'
         if(i.eq.3)char='N'
         if(i.eq.4)char='O'
         

         x = xyz(1,i)/10.0d0
         y = xyz(2,i)/10.0d0
         z = xyz(3,i)/10.0d0
         
         write(7000,FMT2)i,'TTA',char,i, x,y,z
      enddo
      
      do i = 1,numatompertemp
         if(i.eq.1)char='C'
         if(i.eq.2)char='Ca'
         if(i.eq.3)char='N'
         if(i.eq.4)char='O'
         
         
         x = local_template(1,i)/10.0d0
         y = local_template(2,i)/10.0d0
         z = local_template(3,i)/10.0d0
         
         write(7000,FMT2)i,'TTA',char,i, x,y,z
      enddo
      
      
      box_units = box/(10.0d0)
      
      write(7000,FMT3)box_units,box_units,box_units
      close(7000)
    end subroutine xyz_template_centroid


    subroutine xyzf_template_centroid(local_template,len,filename)
      implicit none
      double precision :: local_template(3,len)
      double precision :: box_units
      double precision :: x,y,z
      character(len=*) :: filename
      character(len=*),parameter :: FM1 = "(i5)"
      character(len=*),parameter :: FMT2 = "( i5 , 2a5 , i5 , 3f8.3)"
      character(len=*),parameter :: FMT3 = "(f10.5,f10.5,f10.5)"
      integer :: i,flag
      integer :: typ,type,len,count
      character(len=3) :: char
      
      open(unit = 7000, file = trim(grofile_path)//'/'//filename)
      write(7000,*)'xyz'
      
      write(7000,FM1)2*numatompertemp
      
      count = 0
      do i = 1,numatompertemp
         typ  = position(i)%type
         char = type2string(typ)
         char = adjustl(char)
         char = trim(char)          
         
         if(mod(count,numatompertemp).ne.0)then
            count = count + 1
         else
            count = 1
         endif

         if(count.le.4)then
            char = 'C'
         else
            char ='O'
         endif

         
         x = xyzf(1,i)/10.0d0
         y = xyzf(2,i)/10.0d0
         z = xyzf(3,i)/10.0d0
         
         write(7000,FMT2)i,'TTA',char,i, x,y,z
      enddo
      
      count = 0
      do i = 1,numatompertemp
         typ  = position(i)%type
         char = type2string(typ)
         char = adjustl(char)
         char = trim(char)          
         
         if(mod(count,numatompertemp).ne.0)then
            count = count + 1
         else
            count = 1
         endif
         
         if(count.le.4)then
            char = 'C'
         else
            char ='O'
         endif

         
         x = local_template(1,i)/10.0d0
         y = local_template(2,i)/10.0d0
         z = local_template(3,i)/10.0d0
         
         write(7000,FMT2)i,'TTA',char,i, x,y,z
      enddo
      
      
      box_units = box/(10.0d0)
      
      write(7000,FMT3)box_units,box_units,box_units
      close(7000)
    end subroutine xyzf_template_centroid


     subroutine template_Gro(local_template,len,filename,flag)
      implicit none
      double precision :: local_template(3,len)
      double precision :: box_units
      double precision :: x,y,z
      character(len=*) :: filename
      character(len=*),parameter :: FM1 = "(i5)"
      character(len=*),parameter :: FMT2 = "( i5 , 2a5 , i5 , 3f8.3)"
      character(len=*),parameter :: FMT3 = "(f10.5,f10.5,f10.5)"
      integer :: i,flag
      integer :: typ,type,len,count
      character(len=3) :: char

      open(unit = 7000, file = trim(grofile_path)//'/'//filename)
      write(7000,*)'xyz'
      write(7000,FM1)len
      


      if(flag.eq.0)then

         count = 0
         do i = 1,len
            typ  = position(i)%type
            char = type2string(typ)
            char = adjustl(char)
            char = trim(char)          


      
            if(mod(count,numatompertemp).ne.0)then
               count = count + 1
            else
               count = 1
            endif
            
            if(count.eq.1)char='C'
            if(count.eq.2)char='Ca'
            if(count.eq.3)char='N'
            if(count.eq.4)char='O'

            
            x = local_template(1,i)/10.0d0
            y = local_template(2,i)/10.0d0
            z = local_template(3,i)/10.0d0
            
            write(7000,FMT2)i,'TTA',char,i, x,y,z
         enddo
      endif


      box_units = box/(10.0d0)
      
      write(7000,FMT3)box_units,box_units,box_units
      close(7000)
    end subroutine template_Gro
    
   

    
end module mod_GroFile
