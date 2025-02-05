module mod_init_posit
  use ifport
  use global
  use get_started
  
  

  implicit none
  


contains

  subroutine init_position(situation,file_name)
    
    implicit none
    integer:: counter,junk
    integer:: numpd
    double precision ::gsize
    integer :: i,j,k
    integer :: flag
    character(len=*) :: situation,file_name
    character(len=*),parameter :: FMT10 = "(f7.5,f8.5,f8.5)"
    character(len=*),parameter :: FMT11 = "(f7.5,f8.5,f9.5)"
    character(len=*),parameter :: FMT12 = "(f7.5,f9.5,f8.5)"
    character(len=*),parameter :: FMT13 = "(f7.5,f9.5,f9.5)"

    character(len=*),parameter :: FMT14 = "(f8.5,f9.5,f8.5)"
    character(len=*),parameter :: FMT15 = "(f8.5,f8.5,f9.5)"
    character(len=*),parameter :: FMT16 = "(f8.5,f9.5,f9.5)"
    character(len=*),parameter :: FMT17 = "(f8.5,f8.5,f8.5)"
    character(len=*),parameter :: FMT3 = "(f8.5,f10.5,f10.5)"
    character(len=*),parameter :: FMT4 = "(f9.5,f9.5,f9.5)"
    double precision :: half,deltax,remainder,spot

    print*,'what is np',np
    print*,'what is dens and box:',dens,box

    
    select case(situation)

    case('cubic grid')

       numpd = ceiling(np**(1.0d0/3.0d0))
       gsize = box/(1.0d0*numpd) 
       
       print*,'----------------initializing positions ----------------------'
       print*,'              INITIALIZING PARTICLES ON CUBIC GRID'
       print*,'              GRID SPACING:',gsize
       print*,'              NUMBER PER DIMENSION',numpd
       counter =1
       
       do i=1,numpd
          do j=1,numpd
             do k=1,numpd
                
                if(counter.le.np)then
                   

                   position(counter)%x = (k-1)*gsize
                   position(counter)%y = (j-1)*gsize
                   position(counter)%z = (i-1)*gsize


                   if(mod(counter,2).eq.0)then
                      type(counter) = 1
                   else
                      type(counter) = 2
                   endif
                   
                endif
                counter = counter+1
                
             enddo
          enddo
       enddo


    case('random')
       
       do i = 1,np
          position(i)%x = box*rand();position(i)%y = box*rand(); position(i)%z = box*rand()
          type(i) = 1
       enddo


    case('read from file')
       open(unit = 500, file = file_name)

       read(500,*)box,box,box
       read(500,*)np
       
       print*,'reading in values:',box,box,box,np

       do i= 1,np
          read(500,*)junk,position(i)%x,position(i)%y,position(i)%z
       enddo


       close(500)
    end select





    open(unit = 500,file = 'posit_start.dat')
    write(500,*)box,box,box
    write(500,*)np
    do i = 1,np
       write(500,*)position(i)%x,position(i)%y,position(i)%z
    enddo

    vol = box*box*box
    dens = real(np)/vol

    print*,'new density and volume:',dens,vol
    
   
    
    
    print*,'-------------end initializing particles-------------'
    
    return
  end subroutine init_position

   

 

  subroutine init_velocity
    implicit none
    double precision :: sumv(3), sumv2
    double precision :: fs,fs1
    double precision :: mm,factor
    integer :: i
    integer :: type

    print*,'========================================'
    PRINT*,'initializing velocity'
    print*
    do i = 1 ,np
       v(i)%x = rand()-0.50
       v(i)%y = rand()-0.50
       v(i)%z = rand()-0.50
    enddo
    

    sumv(:) = 0.0d0; sumv2 = 0.0d0

    tot_mass =0.0d0
    do i = 1,np
       type    = position(i)%type
       mm      = mass(type)

       sumv(1) = sumv(1) + mm*v(i)%x
       sumv(2) = sumv(2) + mm*v(i)%y
       sumv(3) = sumv(3) + mm*v(i)%z

       sumv2 = sumv2 + 0.50*mm*(v(i)%x**2 + v(i)%y**2 + v(i)%z**2)
       tot_mass = tot_mass + mm
    enddo
    ncons = ncons + 3

    !---velocity center of mass
    sumv(:) = sumv(:)/tot_mass

    !---mean squared velocity

    !---
    sumv2 = amu2e*sumv2/(0.50*kboltz*(3.0*real(np)-ncons))    
    fs  = sqrt(temp/sumv2)

 
    do i = 1,np
       v(i)%x = (v(i)%x-sumv(1))*fs
       v(i)%y = (v(i)%y-sumv(2))*fs
       v(i)%z = (v(i)%z-sumv(3))*fs
       
    enddo


    !---calcuate initialial kinetic energy
    sumv(:) = 0.0d0; sumv2 = 0.0d0

    sumv2 = 0.0d0
    do i = 1,np
       sumv2 = sumv2 + mass(position(i)%type)*(v(i)%x**2 +v(i)%y**2 +v(i)%z**2)
    enddo
       
    ke = 0.5d0*sumv2*amu2e
    
    print*,'ncons',ncons
    print*,'INITIAL TEMPERATURE:',ke/(0.50*kboltz*(3.0d0*real(np)-ncons))
    

  end subroutine init_velocity


end module mod_init_posit


