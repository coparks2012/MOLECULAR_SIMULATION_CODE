module mod_LRC

  use global
  
  contains

    !--- use this value for shifted potential calculations
    subroutine cutoff_energy
      implicit none
      double precision :: rcut2i,rcut6i,rcut12i,rcut3

      rcut3   = rcut*rcut*rcut
      rcut2i  = 1.0d0/rcut2
      rcut6i  = rcut2i*rcut2i*rcut2i
      rcut12i = rcut6i*rcut6i

      ecut = 4.0d0*(rcut12i-rcut6i)
      ecutAvr = 0.50d0*(4.0d0/3.0d0)*acos(-1.0d0)*rcut3*dens*ecut*np


      print*
      print*,'INITIATING VALUE OF CUT OFF ENERGY'
      print*,'average cutoff energy:',ecutAvr/real(np)
      print*

    end subroutine cutoff_energy


    subroutine long_range_corr_LJ
      implicit none
      double precision :: pi
      double precision :: vol2,vol2i
      double precision :: ri,ri3,ri9

      ecorr = 0.0d0
      pcorr = 0.0d0
      pi    = acos(-1.0d0)
      
      vol2  = vol*vol
      vol2i = 1.0d0/vol2
      
      ri = 1.0d0/rcut
      ri3 = ri**3
      ri9 = ri3**3

      ecorr =   np*(8.0d0/3.0d0)*pi*dens*((1.0d0/3.0d0)*ri9 -ri3)
      pcorr =  (16.0d0/3.0d0)*pi*vol2i*np*np*((2.0d0/3.0d0)*ri9-ri3)

      print*
      print*,'INITIATING VALUE OF CUT OFF CORRECTIONS'
      print*,'energy cut off correction:',ecorr
      print*,'pressure cut off correction:',pcorr

    end subroutine long_range_corr_LJ


    subroutine long_range_corr_init
      implicit none
      double precision :: sig3,ri,ri3,ri9,pi
      double precision :: vol2,vol2i,local_eps
      double precision :: sig2
      double precision :: sig6
      double precision :: rcut3,rcut6,rcut9
      integer :: i,j,np1,np2,k,m
      integer :: alli,allj
      integer :: ioffset
      
      ecorr = 0.0d0
      pcorr = 0.0d0
      pi    = acos(-1.0d0)
   
 
     do i = 1, numatomtype
        
        alli =0
        do j = 1,np
           if(position(j)%type.eq.i)then
              alli = alli + 1
           endif
        enddo

        do j = 1,numatomtype
           
           allj=0
           do m = 1,np
              if(position(m)%type.eq.j)then
                 allj = allj+1
              endif
           enddo
           
           local_eps = sqrt(eps(i)*eps(j))
           !sig2 = 0.5*(sig(i)+sig(j))*0.5*(sig(i)+sig(j))

           sig2 = ((sig(i)+sig(j))*(sig(i)+sig(j)))/(2.0**(1.0/3.0))
           sig6 = sig2*sig2*sig2
           rcut3 = rcut*rcut*rcut
           rcut6 = rcut3*rcut3
           rcut9 = rcut3*rcut6
           ecorr = ecorr + 8.0*pi*alli*allj*local_eps*&
               sig6*(sig6-3.0*rcut6)/(9.0*rcut9)
           pcorr = pcorr + 16.0*pi*alli*allj*local_eps *&
                sig6 * (2.0*sig6 - 3.0*rcut6) / (9.0*rcut9)
        enddo
     enddo
   
  end subroutine long_range_corr_init

  subroutine calc_long_range_correction
    implicit none
    vdw_corr  = ecorr/vol
    pres_corr = pcorr/vol
  end subroutine calc_long_range_correction

end module mod_LRC
