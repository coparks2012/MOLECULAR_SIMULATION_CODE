include 'mkl_vsl.f90'
module mod_langevin
  USE MKL_VSL_TYPE
  use global  
  use ifport
  TYPE (VSL_STREAM_STATE) :: stream
  integer(kind=4) brng,method,seeds,T1,clock_rate,clock_max
contains
  ! -------------------------------------------------------------------------
  ! constant NVE ensemble
  subroutine langevin(dt_partial)
    
  end subroutine langevin

  subroutine langevin_debug_intel
    real(kind=4)    :: a,b
    integer(kind=4) :: remainder,errcode
    integer(kind=4) :: i,j,n,nn
    integer(kind=4) :: start,end
    integer(kind=4) :: T1,T2,clock_rate,clock_max
    real(kind=8)    :: rmass,gamma1,gamma2
    call system_clock(T1,clock_rate,clock_max)
    a = -0.50000; b = 0.50000; n = np
    errcode = vsRngUniform(method,stream,n,r1,a,b)
    errcode = vsRngUniform(method,stream,n,r2,a,b)
    errcode = vsRngUniform(method,stream,n,r3,a,b)

 
    do j =1 ,np
       rmass =  mass(position(j)%type)
       gamma1 = -fric*rmass
       gamma2 = gconst*sqrt(rmass)
       
       
       ff_lang(j)%x =  gamma1*v(j)%x + gamma2*r1(j)
       ff_lang(j)%y =  gamma1*v(j)%y + gamma2*r2(j)
       ff_lang(j)%z =  gamma1*v(j)%z + gamma2*r3(j)

    end do
 
    call system_clock(T2,clock_rate,clock_max)
    time_lan = time_lan + real(T2-T1)/real(clock_rate)
  end subroutine langevin_debug_intel

  subroutine langevin_debug
    implicit none
    integer :: i,j,n,nn
    integer :: T1,T2,clock_rate,clock_max
    double precision   :: rmass,gamma1,gamma2
    double precision   :: rand1,rand2,rand3
    call system_clock(T1,clock_rate,clock_max)

    !stop
    do j = 1,np
       rmass =  mass(position(j)%type)
       gamma1 = -fric*rmass
       gamma2 = gconst*sqrt(rmass)
       
  
       ff_lang(j)%x =  gamma1*v(j)%x + gamma2*(rand()-0.50d0)
       ff_lang(j)%y =  gamma1*v(j)%y + gamma2*(rand()-0.50d0)
       ff_lang(j)%z =  gamma1*v(j)%z + gamma2*(rand()-0.50d0)
       
    end do
    call system_clock(T2,clock_rate,clock_max)
    time_lan = time_lan + real(T2-T1)/real(clock_rate)
  end subroutine langevin_debug

  subroutine init_rand_stream
    brng=VSL_BRNG_MT19937
    !brng=VSL_BRNG_MCG31
    method=VSL_RNG_METHOD_UNIFORM_STD

    call system_clock(T1,clock_rate,clock_max)
    seeds=T1*(rank+1)
    
    print*,'============================================================'
    print*
    print*,'printing seed information'
    print*,'what is seeds',seeds
    print*
    print*,'============================================================'
    
    errcode=vslnewstream( stream, brng,  seeds)
  end subroutine init_rand_stream

end module mod_langevin
