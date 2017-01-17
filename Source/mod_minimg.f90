module mod_minimg
  use global
  contains

    subroutine minimg(dx,dy,dz)
      implicit none
      double precision :: dx,dy,dz


      dx = dx-box*nint(ibox*dx)
      dy = dy-box*nint(ibox*dy)
      dz = dz-box*nint(ibox*dz)


    end subroutine minimg

  end module mod_minimg
