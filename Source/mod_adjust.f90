module mod_adjust

use global
implicit none

contains 
  
  subroutine adjust
    
    implicit none
    double precision :: ratio
    integer :: att
    att = naccept + nreject

    print*
    print*,'what the hell is att',att,naccept,nreject
    print*,'initial kn',kn

    if (att .ne. 0) then
       ratio=(naccept*1.d0)/(att*1.0d0)
       print*,'what the hell is ratio',ratio,naccept,att
       if (ratio.gt.0.50) then
          kn = kn * 1.05
       endif
       if(ratio.lt.0.50)then
          kn = kn * 0.95
       endif
    endif
    
    
    naccept = 0
    nreject = 0

    print*,'we are now leaving adjust',kn
    print*
  end subroutine adjust






end module mod_adjust
