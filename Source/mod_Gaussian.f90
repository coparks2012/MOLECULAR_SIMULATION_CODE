module mod_Gaussian
  use ifport
  implicit none
contains
  
  double precision function Gaussian()
    implicit none
    double precision :: pi
    double precision :: x1,x2
    double precision :: r

    r = 2.0d0
    do while(r.ge.1.)
       x1 = 2.0d0*rand() -1.0d0
       x2 = 2.0d0*rand() -1.0d0
       
       r =x1*x1 + x2*x2
    enddo

    Gaussian = x1*sqrt(-2.0d0*log(r)/r)
    
    
  end function Gaussian


    ! -------------------------------------------------------------------------
  ! Marsaglia RNG
  !  compute in double precision, return single or double
  !  MUST be initialized by calling with iseed > 0
  !  call to return a random # is by calling with iseed = 0
  
  double precision function ranmars(iseed)
    implicit none
    
    ! argument variables
    
    integer iseed
    
    ! local variables
    
    integer ij,kl,i,j,k,l,ii,jj,m,i97,j97,iflag
    double precision u(97),s,t,c,cd,cm,uni
    data iflag /0/
    save
    
    if (iseed.gt.0.or.iflag.eq.0) then
       if (iseed.eq.0)then
          print*,'Uninitialized Marsaglia RNG'
          stop
       endif
       ij = (iseed-1)/30082
       kl = (iseed-1) - 30082*ij
       i = mod(ij/177,177) + 2
       j = mod(ij,177) + 2
       k = mod(kl/169,178) + 1
       l = mod(kl,169)
       do ii = 1,97
          s = 0.0D0
          t = 0.5D0
          do jj = 1,24
             m = mod(mod(i*j,179)*k,179)
             i = j
             j = k
             k = m
             l = mod(53*l+1,169)
             if (mod(l*m,64).ge.32) s = s + t
             t = 0.5*t
          enddo
          u(ii) = s
       enddo
       c = 362436.0D0/16777216.0D0
       cd = 7654321.0D0/16777216.0D0
       cm = 16777213.0D0/16777216.0D0
       i97 = 97
       j97 = 33
       iflag = 1
    endif
    
    uni = u(i97) - u(j97)
    if (uni.lt.0.0) uni = uni + 1.0
    u(i97) = uni
    i97 = i97 - 1
    if (i97.eq.0) i97 = 97
    j97 = j97 - 1
    if (j97.eq.0) j97 = 97
    c = c - cd
    if (c.lt.0.0) c = c + cm
    uni = uni - c
    if (uni.lt.0.0) uni = uni + 1.0
    ranmars = uni
    
    
  end function ranmars
  


  
end module mod_Gaussian
