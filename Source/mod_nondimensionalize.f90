module mod_nondimensionalize

  use global
  contains

    subroutine calc_ref_values
      implicit none
      double precision :: boltzz
      double precision :: local
      integer :: i 
  
      reps = eps(1); rmass = mass(1); rsig = sig(1)
      
      rdeltat = 48.88821D0*rsig*sqrt(rmass/reps)

      !rq = sqrt(rsig*reps/332.0681883)
      rq = sqrt(rsig*reps/332.0636)
      
      boltzz = 0.00198720418d0
      rtemp = reps/boltzz

      print*,'REFERENCE VALUES:'
      print*,'eps:',reps
      print*,'mass:',rmass
      print*,'length:',rsig
      print*,'time:',rdeltat
      print*,'charge:',rq
      print*,'temp:',rtemp


    end subroutine calc_ref_values


    subroutine nondimensionalize_forcefield
      implicit none
      double precision :: x,y,z
      double precision :: sig6,sig12
      integer :: i,k,j
      integer :: atom1
      integer :: offset
      integer :: num,total

  
      

      !---nondimensionalize sig,eps,and mass
      do i = 1,numAtomType
         eps(i)   =  eps(i)/reps
         mass(i)  =  mass(i)/rmass
         sig(i)   =  sig(i)/rsig
         qtype(i) =  qtype(i)/rq
      end do


      do i = 1,nummoltype
         x = FC(i)%coords(1)%x; y = FC(i)%coords(1)%y; z = FC(i)%coords(1)%z;
         do j = 1,numMolAtom(i)
            FC(i)%coords(j)%x = FC(i)%coords(j)%x - x
            FC(i)%coords(j)%y = FC(i)%coords(j)%y - y
            FC(i)%coords(j)%z = FC(i)%coords(j)%z - z

            FC(i)%coords(j)%x = FC(i)%coords(j)%x/rsig
            FC(i)%coords(j)%y = FC(i)%coords(j)%y/rsig
            FC(i)%coords(j)%z = FC(i)%coords(j)%z/rsig
         enddo
      enddo


      total =0 
      do i = 1,nummoltype
         total = total + nummolshake(i)
      enddo

      do i = 1,total
         num = frac_shake(i)%num
         if(num.eq.1)then
            frac_shake(i)%bond = frac_shake(i)%bond/rsig
         elseif(num.eq.3)then
            frac_shake(i)%bond  = frac_shake(i)%bond/rsig
            frac_shake(i)%bond2 = frac_shake(i)%bond2/rsig
         elseif(num.eq.4)then
            frac_shake(i)%bond  = frac_shake(i)%bond/rsig
            frac_shake(i)%bond2 = frac_shake(i)%bond2/rsig
            frac_shake(i)%bond3 = frac_shake(i)%bond3/rsig
         elseif(num.eq.-3)then
            frac_shake(i)%bond  = frac_shake(i)%bond/rsig
            frac_shake(i)%bond2 = frac_shake(i)%bond2/rsig
            frac_shake(i)%bond3 = frac_shake(i)%bond3/rsig
         endif
      end do



      do i = 1 , numAtomType
         do j = 1,numAtomType
            
            offset = numAtomType*(i-1)+j
            
            sig6  = rsig**6
            sig12 = rsig**12 

            lj1(offset)   = lj1(offset)/reps/sig12
            lj2(offset)   = lj2(offset)/reps/sig6

            lj3(offset) = lj3(offset)/reps/sig12
            lj4(offset) = lj4(offset)/reps/sig6

            lj14_1(i,j)   = lj14_1(i,j)/reps/sig12
            lj14_2(i,j)   = lj14_2(i,j)/reps/sig6

            lj14_3(i,j) = lj14_3(i,j)/reps/sig12
            lj14_4(i,j) = lj14_4(i,j)/reps/sig6            
         enddo
      enddo



      do i = 1,numbondtype
         bondcoeff(1,i) = (bondcoeff(1,i)/reps)*rsig*rsig
         bondcoeff(2,i) = bondcoeff(2,i)/rsig
      end do
      

      do i = 1,numangletype
         anglecoeff(1,i) = anglecoeff(1,i)/reps
      enddo


      do i = 1,numdihedtype
         dihedcoeff(1,i) = dihedcoeff(1,i)/reps
         dihedcoeff(2,i) = dihedcoeff(2,i)/reps
         dihedcoeff(3,i) = dihedcoeff(3,i)/reps
         dihedcoeff(4,i) = dihedcoeff(4,i)/reps
      end do

      do i = 1,numimprotype
         improcoeff(1,i) = improcoeff(1,i)/reps
      end do
      

     
    end subroutine nondimensionalize_forcefield


    subroutine non_dimensionalize_input
      implicit none
      double precision :: erfcc,erfcd
      integer :: i,j

      dt = dt/rdeltat
      dthalf = dthalf/rdeltat
      dthalf2 = dthalf2/rdeltat/rdeltat
      dt2 = dt*dt
      rcut = rcut/rsig
      rcut2 = rcut*rcut
      mrcut = mrcut/rsig
      mrcut2 = mrcut*mrcut
      denscut = denscut/rsig
      denscut2 = denscut*denscut

      rv = rv/rsig
      rv2 = rv2/rsig/rsig
      rvrc = rvrc/rsig
      rvrc2 = rvrc2/rsig/rsig

 
      print*,'THIS IS BOX LENGTH BEFORE',box/rsig
      box  = dble(box)/dble(rsig)
      hbox = dble(hbox)/dble(rsig)
      vol  = dble(vol)/dble((rsig*rsig*rsig))
      ibox = dble(ibox) *dble(rsig)
      print*,'after',box



      dens = dens*rsig*rsig*rsig
      mdens = mdens*rsig*rsig*rsig

      temp = temp/rtemp
      STDtemp = STDtemp/sqrt(rtemp)
      anneal_temp = anneal_temp/rtemp

      cut_coul = cut_coul/rsig
      cut_coulsq = cut_coul * cut_coul
      alpha = alpha*rsig
      gewald = gewald*rsig

      tempa_noprune(1)%temp(:,:) = tempa_noprune(1)%temp(:,:)/rsig
      tempa_noprune(2)%temp(:,:) = tempa_noprune(2)%temp(:,:)/rsig
      tempa_noprune(3)%temp(:,:) = tempa_noprune(3)%temp(:,:)/rsig
      tempa_noprune(4)%temp(:,:) = tempa_noprune(4)%temp(:,:)/rsig
      tempa_noprune(5)%temp(:,:) = tempa_noprune(5)%temp(:,:)/rsig
      tempa_noprune(6)%temp(:,:) = tempa_noprune(6)%temp(:,:)/rsig

      !---nondimensionalize scale for molecular pruning
      CObond = CObond/rsig


      !---nondimensionalize shake
      tolerance = tolerance/rsig

      !---nondimensionalize andersen thermostat
      fric   = 0.0005*rdeltat
      gconst =  sqrt(12.0*kboltz*temp*fric/dthalf)

      print*,'what the fuck is gconst',gconst
      


      print*,'JUST BEFORE E_SHIFT F_SHIFT:',e_shift,f_shift
      print*,'REPS AND RSIG:',reps,rsig
      !e_shift = e_shift/reps
      !f_shift = f_shift*rsig/reps
      
      e_shift = e_shift*rsig
      f_shift = f_shift*rsig*rsig

      print*,'IN NONDIM:',e_shift,f_shift
      print*,'what is cut_coul:',cut_coul
      print*,'what is cut_coulsq:',cut_coulsq
      print*,'what is MY_PIS:',MY_PIS
      print*,'what is alpha:',alpha
  

      print*,'-----------------------------------------------------'
      print*,'NON DIMENSIONALZIED INPUT:'
      print*,'dt:',dt
      print*,'dthalf2:',dthalf2
      print*,'dthalf:',dthalf
      print*,'rcut:',rcut
      print*,'rv:',rv,rv2
      print*,'rcrv2',rvrc2
      print*,'cut_coul',cut_coul
      print*,'mrcut',mrcut
      print*,'mrcut2',mrcut2
      print*,'denscut',denscut
      print*,'denscut2',denscut2
      print*,'dens DID NOT NONDIMENSIONALIZE:',dens
      print*,'temp',temp
      print*,'kn',kn
      print*,'anneal_temp',anneal_temp
      print*,'E_SHIFT AND F_SHIFT:',e_shift,f_shift
      print*,'WE ARE NOW DONE NONDIMENSIONALIZING'
    end subroutine non_dimensionalize_input


  end module mod_nondimensionalize
