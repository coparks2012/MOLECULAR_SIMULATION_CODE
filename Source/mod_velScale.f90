module mod_velScale
  use global
  
  implicit none

  contains

    subroutine velScale
      implicit none
      double precision :: sumv2,fs
      integer :: i

      sumv2 = 0.0d0
      do i =1,np
         sumv2 = sumv2 + mass(type(i))*(v(i)%x*v(i)%x+v(i)%y*v(i)%y+v(i)%z*v(i)%z)
      enddo

      sumv2 = sumv2/real(np)

      fs = sqrt(3.0d0*temp/sumv2)

      do i=1,np
         v(i)%z = v(i)%z*fs; v(i)%y = v(i)%y*fs; v(i)%x = v(i)%x*fs
      enddo
         
    end subroutine velScale

  end module mod_velScale
