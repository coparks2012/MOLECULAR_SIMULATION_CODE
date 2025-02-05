module mod_ensemble
  use global
  use mod_Gaussian
  use mod_Grofile
 
   
contains
  
  ! -------------------------------------------------------------------------
  ! constant NVE ensemble
  subroutine nve_v(dt_partial)
    implicit none
    real*4 :: dt_partial
    double precision :: dtmass
    double precision :: mm
    integer :: i
    integer :: T1,T2,clock_rate,clock_max
    
    call system_clock(T1,clock_rate,clock_max)
    ke = 0.0d0

    do i=1,np


       dtmass = dt_partial/mass(position(i)%type)
       
       
       v(i)%x = v(i)%x + dtmass*ff(i)%x
       v(i)%y = v(i)%y + dtmass*ff(i)%y
       v(i)%z = v(i)%z + dtmass*ff(i)%z
       
       
       ke = ke+ mass(position(i)%type)*(v(i)%x*v(i)%x+v(i)%y*v(i)%y+v(i)%z*v(i)%z) 

    enddo
    ke = ke*0.50*amu2e
    call system_clock(T2,clock_rate,clock_max)
    time_ver = time_ver + real(T2-T1)/real(clock_rate)

  end subroutine nve_v
  



  ! -------------------------------------------------------------------------
  ! constant NVE ensemble
  
  subroutine nve_x
    implicit none
    double precision :: dx,dy,dz
    double precision :: dtmass2
    double precision :: mm
    integer:: i
    integer:: T1,T2,clock_rate,clock_max

    call system_clock(T1,clock_rate,clock_max)
  
    do i = 1,np


       dtmass2 = dthalf2/mass(position(i)%type)
     

       dx = dt*v(i)%x
       dy = dt*v(i)%y
       dz = dt*v(i)%z

       drx(i) = drx(i) + dx
       dry(i) = dry(i) + dy
       drz(i) = drz(i) + dz
       
       position(i)%x = position(i)%x + dx
       position(i)%y = position(i)%y + dy
       position(i)%z = position(i)%z + dz

      ! if(mol(1,i).eq.1)then
      !    ptta( mol(2,i) )%x = ptta( mol(2,i) )%x + dx 
      !    ptta( mol(2,i) )%y = ptta( mol(2,i) )%y + dy 
      !    ptta( mol(2,i) )%z = ptta( mol(2,i) )%z + dz
      ! endif
    enddo

    call system_clock(T2,clock_rate,clock_max)

    time_ver = time_ver + real(T2-T1)/real(clock_rate)
    
  end subroutine nve_x


     ! -------------------------------------------------------------------------
  ! constant NVE ensemble
  subroutine nve_v_integrate_flag(dt_partial)
    implicit none
    double precision :: dtmass
    real*4 :: dt_partial
    double precision :: mm
    integer :: i
    integer :: T1,T2,clock_rate,clock_max
    
    call system_clock(T1,clock_rate,clock_max)
    ke = 0.0d0
    do i=1,np

       if(integrate_flag(i).eq.1)then
          dtmass = dt_partial/mass(position(i)%type)
          
          
          v(i)%x = v(i)%x + dtmass*ff(i)%x
          v(i)%y = v(i)%y + dtmass*ff(i)%y
          v(i)%z = v(i)%z + dtmass*ff(i)%z
          
          
          ke = ke+ 0.50d0*mass(position(i)%type)*(v(i)%x*v(i)%x+v(i)%y*v(i)%y+v(i)%z*v(i)%z) 
       endif
    enddo
    
    call system_clock(T2,clock_rate,clock_max)
    time_ver = time_ver + real(T2-T1)/real(clock_rate)
    
  end subroutine nve_v_integrate_flag
  

  ! -------------------------------------------------------------------------
  ! constant NVE ensemble
  
  subroutine nve_x_integrate_flag
    implicit none
    double precision :: dx,dy,dz
    double precision :: dtmass2
    double precision :: mm
    integer:: i
    integer:: T1,T2,clock_rate,clock_max

    call system_clock(T1,clock_rate,clock_max)
  
    do i = 1,np
     
       if(integrate_flag(i).eq.1)then
          dx = dt*v(i)%x
          dy = dt*v(i)%y
          dz = dt*v(i)%z
          
          
          drx(i) = drx(i) + dx
          dry(i) = dry(i) + dy
          drz(i) = drz(i) + dz
          
          position(i)%x = position(i)%x + dx
          position(i)%y = position(i)%y + dy
          position(i)%z = position(i)%z + dz
          
       !   if(mol(1,i).eq.1)then
       !      ptta( mol(2,i) )%x = ptta( mol(2,i) )%x + dx 
       !      ptta( mol(2,i) )%y = ptta( mol(2,i) )%y + dy 
       !      ptta( mol(2,i) )%z = ptta( mol(2,i) )%z + dz
       !   endif
       endif

    enddo
    
    call system_clock(T2,clock_rate,clock_max)

    time_ver = time_ver + real(T2-T1)/real(clock_rate)
    
  end subroutine nve_x_integrate_flag
  
  subroutine berendsen
    implicit none
    double precision :: mu,currx,curry,currz,dx,dy,dz
    integer :: step
    integer :: i
 
    mu = 1 -(dt/tau)*(ptarget-pressure)
    mu = mu**(1.0d0/3.0d0)
    box = box*mu
    vol = box*box*box
    hbox = 0.50*box
    ibox = 1.0/box

    do i =1,np
       currx = position(i)%x
       curry = position(i)%y
       currz = position(i)%z

       position(i)%x = position(i)%x*mu
       position(i)%y = position(i)%y*mu
       position(i)%z = position(i)%z*mu

       dx = position(i)%x-currx
       dy = position(i)%y-curry
       dz = position(i)%z-currz


       drx(i) = drx(i) + dx
       dry(i) = dry(i) + dy
       drz(i) = drz(i) + dz
       !if(mol(1,i).eq.1)then
       !   ptta( mol(2,i) )%x = ptta( mol(2,i) )%x + dx 
       !   ptta( mol(2,i) )%y = ptta( mol(2,i) )%y + dy 
       !   ptta( mol(2,i) )%z = ptta( mol(2,i) )%z + dz
       !endif

    end do

  end subroutine berendsen

 


  subroutine thermalize
    implicit none
    double precision :: x1
    double precision :: MassInv
    integer :: i
    do i = 1 ,np
       MassInv = 1.0d0/mass(position(i)%type)
       MassInv = sqrt(MassInv)
       
       v(i)%x = Gaussian()*STDtemp*MassInv
       v(i)%y = Gaussian()*STDtemp*MassInv
       v(i)%z = Gaussian()*STDtemp*MassInv
    enddo
  end subroutine thermalize

  subroutine anderson_thermo
    implicit none
    double precision :: x1
    double precision :: MassInv
    integer :: i


    do i = 1 ,np
       if(rand().lt.0.001d0)then
          MassInv = 1.0d0/mass(position(i)%type)
          MassInv = sqrt(MassInv)
          
          v(i)%x = Gaussian()*STDtemp*MassInv
          v(i)%y = Gaussian()*STDtemp*MassInv
          v(i)%z = Gaussian()*STDtemp*MassInv
       endif
    enddo
  end subroutine anderson_thermo


  subroutine scale_velocity
    implicit none
    double precision :: sum,calcTemp
    double precision :: scale
    integer :: i
    
    sum = 0.0d0
    do i = 1, np
       sum = sum + mass(position(i)%type)*(v(i)%x*v(i)%x+v(i)%y*v(i)%y+v(i)%z*v(i)%z)
    enddo
    
    sum = sum*0.50d0
    calcTemp = 2.0d0*sum/(kboltz*(3.0d0*np-ncons))
    ! calcTemp = calcTemp/(kboltz*r)
    
    ! scale = sqrt(calcTemp/temp)
    scale  = sqrt(temp/calcTemp)
    !print*,'sum kboltz',sum,kboltz,temp,calcTemp,scale
    
    do i = 1,np
       v(i)%x = v(i)%x*scale
       v(i)%y = v(i)%y*scale
       v(i)%z = v(i)%z*scale
    enddo
  end subroutine scale_velocity
  

 
    

  
end module mod_ensemble
