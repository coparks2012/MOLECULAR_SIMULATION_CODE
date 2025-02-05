module get_started

  use global
  use mod_memory
  implicit none
contains
  

  subroutine init_variables_wrapper
    implicit none
    call calc_number_of_particles
    call calculate_gewald
    call calculate_total_number_bonds
    call calculate_total_number_angles
    call calculate_total_number_dihed
    call calculate_total_number_impro
    call calculate_total_number_shake
    call init

    call print_input

  end subroutine init_variables_wrapper


  subroutine init
    implicit none
      double precision :: erfcc,erfcd
    integer :: i

    !----LETS DEFINE CONSTANTS
    !====kcal/mol  to AMU
    kcal2amu  = 0.0004184 
    amu2e     = 2390.057361
    coulpre   = 332.06371
    nktv2p    = 68568.415

    !---units of kboltz are kcal/mol/k
    kboltz    = 0.001987204181D0

    !---map kboltz(kcal/mol) to AMU
    !kboltz    = kboltz*dt_factor 
    !kboltz = 0.00000008314


    rv = neighbor_pad+cut_coul
    rv2 = rv*rv
    rvrc =  neighbor_pad
    rvrc2 = rvrc*rvrc
    !rvrc2 = 0.250d0*rvrc2
      
    
    dthalf  = 0.50d0*dt*kcal2amu
    dthalf2 = kcal2amu*dt*dt*0.50d0

    print*,'HELLO FROM RANK TEMP',rank,temp
    STDtemp = sqrt(temp)
    beta   = 1.0d0/(kboltz*temp)
    nu = 10

  
    dens = real(np)/(box**3)
    mdens = numMolArray(1)/(box**3)

    vol = box**3
    hbox   = box/(2.0d0)
    ibox   = 1.0d0/box
    
    
    !---shake constaint values
    maxiter   = 25
    tolerance = 0.0001


    rcut2      = rcut*rcut
    cut_coulsq = cut_coul*cut_coul
    mrcut2     = mrcut*mrcut
    denscut2    = denscut*denscut


    !---thermostat variables
    fric   = 0.0005
    gconst =  sqrt(12.0*kboltz*kcal2amu*temp*fric/(0.50*dt))
    fric   = fric*amu2e
    gconst = gconst*amu2e
    viscous_damp = 0.005
    viscous_damp = viscous_damp*amu2e

    !---barostat variables
    vmax = 0.05
    tau = 1000.0
    !---multiply tau by bulk modulus
    tau = tau*21454.85352


    !--density histogram
    num_dens_ob = num_dens_ob + 1




    !---DSF parameters
    EWALD_F =  1.12837917
    EWALD_P =  0.3275911
    A1      =  0.254829592
    A2      = -0.284496736
    A3      =  1.421413741
    A4      = -1.453152027
    A5      =  1.061405429
    MY_PIS =   1.77245385090551602729
    MY_PIS_INV = 1.0d0/MY_PIS

    erfcc = erfc(alpha*cut_coul)
    erfcd = exp(-alpha*alpha*cut_coul*cut_coul)

    f_shift = -(erfcc/cut_coulsq + 2.0/MY_PIS*alpha*erfcd/cut_coul)
    e_shift =  erfcc/cut_coul - f_shift*cut_coul


    !---initialize scale value for pruning
    CObond = 0.615
    

    !---spring constant
    kn = 0.025*kboltz*temp
    
    print*,'============================================='
    print*,'total number of atoms',np

    do i = 1,nummoltype
       print*,'number of molecules type',i,nummolarray(i)
    enddo
   
    print*,'box length',box
    print*,'IN GET STARTED:',e_shift,f_shift
    print*,'what is cut_coul:',cut_coul
    print*,'what is cut_coulsq:',cut_coulsq
    print*,'what is MY_PIS:',MY_PIS
    print*,'what is alpha:',alpha
    print*,'what is erfcc:',erfcc
    print*,'what is erfcd:',erfcd
    PRINT*,'what is gewald',gewald
    print*
    print*
    print*,'========================================'
  end subroutine init

  subroutine calculate_gewald
    implicit none
    double precision :: ewaldcof
    double precision :: ewaldcof_lo, ewaldcof_hi
    integer :: i
    double precision :: tolerance
    
    tolerance = 0.000001
    ewaldcof = 1.0
    print*,erfc(ewaldcof*cut_coul)/cut_coul,cut_coul
    
    do while ( erfc(ewaldcof*cut_coul)/cut_coul.gt.tolerance ) 
       gewald = gewald* 2.0
    enddo
    ewaldcof_lo = 0.0
    ewaldcof_hi = ewaldcof
    do i =1,100
       ewaldcof = 0.5 * ( ewaldcof_lo + ewaldcof_hi )
       if ( erfc(ewaldcof*cut_coul)/cut_coul.gt.tolerance ) then
          ewaldcof_lo = ewaldcof
          
       else 
          ewaldcof_hi = ewaldcof
       endif
    enddo
    gewald =ewaldcof
  end subroutine calculate_gewald

  subroutine calc_number_of_particles
    implicit none
    integer :: i,j

    !---calculate the number of atoms total in simulation
    np = 0
    
    do i = 1,nummoltype
       do j = 1,numMolArray(i)
          np = np + numMolAtom(i)
       enddo
    enddo
  end subroutine calc_number_of_particles


  
  !-----------------------------------------------------------------------------------------!
    !
    ! calculates the total number of bonds in the system 
    subroutine calculate_total_number_bonds
      implicit none
      integer :: i,j
      
      !---zero numbond and numangle
      numbond  = 0

      do i = 1,numMolType
         do j = 1,numMolAtom(i)

            numbond = numbond + numMolArray(i)*Mol2NumBonds(i)%numbonds(1,j)


         enddo
      enddo

      !---newtons 3rd law
      numbond = numbond/2

      print*,'TOTAL NUMBER OF 12 BONDS:',numbond
      call totalbond_memory

    end subroutine calculate_total_number_bonds


    !-----------------------------------------------------------------------------------------!
    !
    ! calculates the total number of bonds in the system 
    subroutine calculate_total_number_angles
      implicit none
      integer :: i,j
      
      !---zero numbond and numangle
      numangle = 0

      do i = 1,numMolType
         do j = 1,numMolAtom(i)

            numangle = numangle + numMolArray(i)*Mol2NumBonds(i)%numbonds(2,j)

         enddo
      enddo

      !---newtons 3rd law
      numangle = numangle/2
      print*,'TOTAL NUMBER OF 1-3 BONDS:',numangle

      call totalangle_memory
    end subroutine calculate_total_number_angles


    subroutine calculate_total_number_dihed
      implicit none
      integer :: i,j
      
      !---zero numbond and numangle
      numdihed = 0

      do i = 1,numMolType
         do j = 1,numMolAtom(i)

            numdihed = numdihed + numMolArray(i)*Mol2NumBonds(i)%numbonds(3,j)

         enddo
      enddo

      !---newtons 3rd law
      numdihed = numdihed/2
      print*,'TOTAL NUMBER OF 1-4 BONDS',numdihed

      call totaldihed_memory

    end subroutine calculate_total_number_dihed

    subroutine calculate_total_number_impro
      implicit none
      integer :: i,j
      
      !---zero numbond and numangle
      numimpro = 0

      do i = 1,numMolType
         do j = 1,numMolAtom(i)

            numimpro = numimpro + numMolArray(i)*Mol2NumBonds(i)%numbonds(4,j)

         enddo
      enddo

      !---newtons 3rd law
      print*,'TOTAL NUMBER OF IMPRO BONDS',numimpro

      call improdihed_memory

    end subroutine calculate_total_number_impro
  

    
    subroutine calculate_total_number_shake
      implicit none
      integer :: i,j
      
      !---zero numbond and numangle
      nshake = 0

      do i = 1,nummoltype
         nshake = nshake + nummolarray(i)*nummolshake(i)
      end do

      !---newtons 3rd law
      write(101,*)'TOTAL NUMBER OF shake constraints',nshake

      call total_shake_memory

    end subroutine calculate_total_number_shake


  !-------------------------------------------!
  ! prints input parameters to screen
  !-------------------------------------------!
  subroutine print_input
    
    implicit none
    
  
    return
  end subroutine print_input


 

  


  
  
end module get_started
