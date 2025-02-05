module mod_pbc
use global

contains 
  
  subroutine pbc(step)
    implicit none
    integer :: step,i 
    double precision :: box_low
    
    box_low = 0.0d0
    do i=1,np
       
       do while (position(i)%x.gt.box) 
          position(i)%x = position(i)%x - box 
       enddo
       do while(position(i)%x.le.box_low)                     
          position(i)%x = position(i)%x + box     
       enddo

 
       do while (position(i)%y.gt.box) 
          position(i)%y = position(i)%y - box
       enddo
       do while (position(i)%y.le.box_low) 
          position(i)%y = position(i)%y + box     
       enddo
 
       do while (position(i)%z.gt.box) 
          position(i)%z = position(i)%z - box 
       enddo
       do while  (position(i)%z.le.box_low) 
          position(i)%z = position(i)%z + box     
       enddo
      
    enddo
   

  end subroutine pbc


  subroutine pbc_particle(i)
    implicit none
    integer ::i 
    
    if (position(i)%x > box) then
       position(i)%x = position(i)%x - box 
    elseif  (position(i)%x < 0.d0) then
       position(i)%x = position(i)%x + box     
    endif
    
    
    if (position(i)%y > box) then
       position(i)%y = position(i)%y - box 
    elseif  (position(i)%y < 0.d0) then
       position(i)%y = position(i)%y + box     
    endif
    
    if (position(i)%z > box) then
       position(i)%z = position(i)%z - box 
    elseif  (position(i)%z < 0.d0) then
       position(i)%z = position(i)%z + box     
    endif
    
  end subroutine pbc_particle

  
   subroutine pbc_triclinic
    implicit none
    double precision :: dx,dy,dz
    double precision :: box_x,box_y,box_z
    double precision :: x1,y1,z1
    double precision :: xbound_lo,xbound_hi
    integer :: i


    do i = 1,np
       x1 =position(i)%x
       y1 =position(i)%y
       z1 =position(i)%z
       
       if (z1.lt.zlo)then 
          z1 = z1 + zprd
          y1 = y1 + yz
          x1 = x1 + xz
        
       elseif(z1.gt.zhi)then
          z1 = z1 - zprd
          y1 = y1 - yz
          x1 = x1 - xz
         
       endif
    
       if (y1.lt.ylo)then
          y1 = y1 + yprd
          x1 = x1 + xy
       elseif(y1.gt.yhi)then
          y1 = y1 - yprd
          x1 = x1 - xy
       endif

       xbound_hi = (xz/zprd)*z1+xprd
       xbound_lo = (xz/zprd)*z1

       if(x1.lt.xbound_lo)then
          x1 = x1 + xprd
       elseif(x1.gt.xbound_hi)then
          x1 = x1 -xprd
       endif


    

       position(i)%x=x1
       position(i)%y=y1
       position(i)%z=z1
       
    enddo
  end subroutine pbc_triclinic
  
  subroutine minimum_image_triclinic(dx,dy,dz)
    implicit none
    double precision :: dx,dy,dz
    double precision :: xy,xz,yz
    double precision :: box_x,box_y,box_z
    double precision :: xprd_half,yprd_half,zprd_half


    !---
    ! for 1024 bulk glycine crystal
   !  box_x = 40.8160000
   !  box_y = 47.88036
   !  box_z = 40.5659281
     
   !  xy = -8.0715716 
   !  xz = 0.000000 
   !  yz = 0.00000

    xprd_half = 0.5*xprd
    yprd_half = 0.5*yprd
    zprd_half = 0.5*zprd


    
    if (abs(dz) > zprd_half) then
       if (dz < 0.0)then 
          dz = dz + zprd
          dy = dy + yz
          dx = dx + xz
       else 
          dz = dz - zprd
          dy = dy - yz
          dx = dx - xz
       endif
    endif
    
    if (abs(dy) > yprd_half)then
       if (dy < 0.0)then
          dy = dy + yprd
          dx = dx + xy
       else 
          dy = dy - yprd
          dx = dx - xy
       endif
    endif
    
    if (abs(dx) > xprd_half) then
       if (dx < 0.0) then
          dx = dx + xprd
       else 
          dx = dx - xprd
       endif
    endif
  end subroutine minimum_image_triclinic

    subroutine minimum_image_triclinic_cluster(dx, dy,dz, shiftx,shifty,shiftz)
    implicit none
    double precision :: dx,dy,dz
    double precision :: box_x,box_y,box_z
    double precision :: shiftx,shifty,shiftz
    double precision :: xprd_half,yprd_half,zprd_half

    !---
    ! for 1024 bulk glycine crystal
   !  box_x = 40.8160000
   !  box_y = 47.88036
   !  box_z = 40.5659281
     
   !  xy = -8.0715716 
   !  xz = 0.000000 
   !  yz = 0.00000

    xprd_half = 0.5*xprd
    yprd_half = 0.5*yprd
    zprd_half = 0.5*zprd

    shiftx=0.0d0
    shifty=0.0d0
    shiftz=0.0d0
    
    if (abs(dz) > zprd_half) then
       if (dz < 0.0)then 
          shiftz = shiftz + zprd
          shifty = shifty + yz
          shiftx = shiftx + xz
       else 
          shiftz = shiftz - zprd
          shifty = shifty - yz
          shiftx = shiftx - xz
       endif
    endif
    
    if (abs(dy) > yprd_half)then
       if (dy < 0.0)then
          shifty = shifty + yprd
          shiftx = shiftx + xy
       else 
          shifty = shifty - yprd
          shiftx = shiftx - xy
       endif
    endif
    
    if (abs(dx) > xprd_half) then
       if (dx < 0.0) then
          shiftx = shiftx + xprd
       else 
          shiftx = shiftx - xprd
       endif
    endif
  end subroutine minimum_image_triclinic_cluster

  
  subroutine closest_image_triclinic2(i,j, xnew,ynew,znew)
    implicit none
    double precision :: dx,dy,dz
    double precision :: xy,xz,yz
    double precision  :: xnew,ynew,znew,x0,y0,z0
    double precision :: xprd_half,yprd_half,zprd_half
    integer :: i,j


    dx = position(j)%x-position(i)%x
    dy = position(j)%y-position(i)%y
    dz = position(j)%z-position(i)%z

    x0 = position(j)%x
    y0 = position(j)%y
    z0 = position(j)%z
    
    

    xprd_half = 0.5*xprd
    yprd_half = 0.5*yprd
    zprd_half = 0.5*zprd
    
    if(abs(dz).gt.zprd_half.or.abs(dy).gt.yprd_half.or.abs(dx).gt.xprd_half)then
       
       if (abs(dz).gt.zprd_half) then
          if (dz.lt.0.0)then 
             z0 = z0 + zprd
             y0 = y0 + yz
             x0 = x0 + xz
          else 
             z0 = z0 - zprd
             y0 = y0 - yz
             x0 = x0 - xz
          endif
       endif
       
       
       if (abs(dy).gt.yprd_half)then
          if (dy.lt.0.0)then
             y0 = y0 + yprd
             x0 = x0 + xy
          else 
             y0 = y0 - yprd
             x0 = x0 - xy
          endif
       endif
       
       if (abs(dx).gt.xprd_half) then
          if (dx.lt.0.0) then
             x0 = x0 + xprd
          else 
             x0 = x0 - xprd
          endif
       endif
    else
       dx =0.0d0
       dy =0.0d0
       dz =0.0d0
    endif


    xnew = x0
    ynew = y0
    znew = z0

    !xnew = xnew + dx
    !ynew = ynew + dy
    !znew = znew + dz
    
  end subroutine closest_image_triclinic2
  



end module mod_pbc
