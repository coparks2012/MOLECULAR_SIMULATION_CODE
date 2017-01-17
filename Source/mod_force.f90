module mod_force
  use global
  use mod_minimg
  use mod_force_bond
  use mod_force_angle
  use mod_force_dihed
  use mod_force_dihed_amber
  use mod_force_impro_amber
  use mod_force_correction
  use mod_force_nonbond
  use mod_langevin
  use omp_lib

    interface
     subroutine cal_e_pme_kspace(lx, ly, lz, x, y, z, cg, e_elec_kspace, fx, fy, fz,virial) bind (c)
       use iso_c_binding
       real ( c_float ), VALUE :: lx,ly,lz
       real ( c_float ) ::  x(*),y(*),z(*),cg(*),fx(*),fy(*),fz(*)
       real ( c_double ) ::  e_elec_kspace
       real ( c_double ) ::  virial(*)

     end subroutine cal_e_pme_kspace
  end interface


contains

  subroutine force_driver(step,neighbor_flag)
    implicit none
    integer :: step,neighbor_flag
    
    !========PME with MIC
    if(intel_flag.eq.1.and.coul_flag.eq.2)then
       call force_wrapper_pme_mic(step,neighbor_flag)

    !========PME with no MIC
    elseif(intel_flag.eq.0.and.coul_flag.eq.2)then
       call force_wrapper_PME_NOMIC(step,neighbor_flag)

    !========DSF with MIC
    elseif(intel_flag.eq.1.and.coul_flag.eq.1)then
       call force_wrapper_DSF_MIC(step,neighbor_flag)


    !========DSF no MIC
    elseif(intel_flag.eq.0.and.coul_flag.eq.1)then
       call force_wrapper_DSF_NOMIC(step,neighbor_flag)
    endif

  end subroutine force_driver


  subroutine force_wrapper_pme_mic(step,neighbor_flag)
    implicit none
    double precision :: virial(9)
    integer :: step
    integer :: neighbor_flag
    integer :: i
    integer :: sig=1
    integer :: tid
    integer :: T1,T2,T3,T4,T5,T6,T7,clock_rate,clock_max
    call zero_values(step)
    
    call system_clock(T1,clock_rate,clock_max)
    !========== begin of mic region ===================
    !dir$ offload  target(mic:ourmic), signal(sig)                                &
    !dir$ in(position: alloc_if(.false.),free_if(.false.)),                       &
    !dir$ inout(ff: alloc_if(.false.), free_if(.false.)),                         &
    !dir$ nocopy(nlist: alloc_if(.false.) free_if(.false.)),                      &
    !dir$ nocopy(numneigh: alloc_if(.false.) free_if(.false.)),                   &
    !dir$ in(q: alloc_if(.false.) free_if(.false.)),                              &
    !dir$ nocopy(lj1: alloc_if(.false.) free_if(.false.)),                        &
    !dir$ nocopy(lj2: alloc_if(.false.) free_if(.false.)),                        &
    !dir$ nocopy(lj3: alloc_if(.false.) free_if(.false.)),                        &
    !dir$ nocopy(lj4: alloc_if(.false.) free_if(.false.)),                        &
    !dir$ in(rcut2),in(cut_coulsq),in(box),in(ibox),in(np),                       &
    !dir$ in(gewald),in(EWALD_P),in(MY_PIS_INV),in(EWALD_F),in(coulpre)           &
    !dir$ inout(potential), inout(e_coul),in(step),                               &
    !dir$ in(numAtomType),in(neigh_alloc),                                        &
    !dir$ inout(virialx,virialy,virialz),                                         &
    !dir$ in(A1),in(A2),in(A3),in(A4),in(A5)                                   
    call lj_cut_coul_long_nonewton_MIC
    call system_clock(T2,clock_rate,clock_max)
    time_nonbond = time_nonbond + real(T2-T1)/real(clock_rate)

    
    !========begin of CPU region===================
    !======== begin PME region=====================
    
    call system_clock(T3,clock_rate,clock_max)
    !---convert to PME arrays
    do i = 1,np
       pmex(i) = position(i)%x
       pmey(i) = position(i)%y
       pmez(i) = position(i)%z

    enddo
    call cal_e_pme_kspace(box, box, box, pmex, pmey, pmez, q, e_coul_long, pmefx, pmefy, pmefz,virial)
    virialx_host = virial(1)
    virialy_host = virial(5)
    virialz_host = virial(9)
    !==== the next do loop is commented out currently
    !==== when you call the PME driver, uncomment this out
    do i = 1,np
       ffh(i)%x = ffh(i)%x + pmefx(i)
       ffh(i)%y = ffh(i)%y + pmefy(i)
       ffh(i)%z = ffh(i)%z + pmefz(i)
    end do
    call system_clock(T4,clock_rate,clock_max)
    time_pme = time_pme + real(T4-T3)/real(clock_rate)


    !======= end PME region =================
    call system_clock(T5,clock_rate,clock_max)
    !$omp parallel
    !$omp single

    !$omp task 
    if(numbond.gt.0)then
       call bond_harmonic(step)
       call correct_12_bonds
    endif
    !$omp end task

    
    !$omp task
    if(numangle.gt.0)then
       call angle_harmonic(step)
       call correct_13_bonds
    endif
    !$omp end task
    
    !$omp task
    if(numdihed.gt.0)then
       call dihed_PME_amber
       call correct_14_bonds
    endif
    !$omp end task

    !$omp task
    if(numimpro.gt.0)then
       call impro_amber_force
    endif
    !$omp end task

    !$omp task
    call langevin_debug_intel
    !$omp end task
    
    !$omp end single
    !$omp end parallel
    call system_clock(T6,clock_rate,clock_max)
    time_bond_total = time_bond_total + real(T6-T5)/real(clock_rate)

    !=======end of region============================
    !dir$ offload_wait target(mic:ourmic) wait(sig)


    !==========merge CPU and MIC
    do i = 1,np
       ff(i)%x = ff(i)%x + ffh(i)%x +ff_bond(i)%x + ff_angle(i)%x + ff_dihed(i)%x+&
            ff_impro(i)%x + ff_lang(i)%x

       ff(i)%y = ff(i)%y + ffh(i)%y +ff_bond(i)%y + ff_angle(i)%y + ff_dihed(i)%y+&
            ff_impro(i)%y + ff_lang(i)%y

       ff(i)%z = ff(i)%z + ffh(i)%z +ff_bond(i)%z + ff_angle(i)%z + ff_dihed(i)%z+&
            ff_impro(i)%z + ff_lang(i)%z
    enddo
    e_coul = e_coul+e_coul_host+e_coul_bond+e_coul_angle+e_coul_dihed

    virialx = virialx + virialx_host + v_bondx + v_anglex + v_dihedx + v_improx
    virialy = virialy + virialy_host + v_bondy + v_angley + v_dihedy + v_improy
    virialz = virialz + virialz_host + v_bondz + v_anglez + v_dihedz + v_improz

    call system_clock(T7,clock_rate,clock_max)
    time_force = time_force + real(T7-T1)/real(clock_rate)

  end subroutine force_wrapper_pme_mic

 subroutine force_wrapper_pme_mic_noasync(step,neighbor_flag)
    implicit none
    double precision :: virial(9)
    integer :: step
    integer :: neighbor_flag
    integer :: i
    integer :: sig=1
    integer :: tid
    integer :: T1,T2,T3,T4,T5,T6,T7,clock_rate,clock_max
    call zero_values(step)
    
    call system_clock(T1,clock_rate,clock_max)
    !========== begin of mic region ===================
    !dir$ offload  target(mic:ourmic),                                            &
    !dir$ in(position: alloc_if(.false.),free_if(.false.)),                       &
    !dir$ inout(ff: alloc_if(.false.), free_if(.false.)),                         &
    !dir$ nocopy(nlist: alloc_if(.false.) free_if(.false.)),                      &
    !dir$ nocopy(numneigh: alloc_if(.false.) free_if(.false.)),                   &
    !dir$ in(q: alloc_if(.false.) free_if(.false.)),                              &
    !dir$ nocopy(lj1: alloc_if(.false.) free_if(.false.)),                        &
    !dir$ nocopy(lj2: alloc_if(.false.) free_if(.false.)),                        &
    !dir$ nocopy(lj3: alloc_if(.false.) free_if(.false.)),                        &
    !dir$ nocopy(lj4: alloc_if(.false.) free_if(.false.)),                        &
    !dir$ in(rcut2),in(cut_coulsq),in(box),in(ibox),in(np),                       &
    !dir$ in(gewald),in(EWALD_P),in(MY_PIS_INV),in(EWALD_F),in(coulpre)           &
    !dir$ inout(potential), inout(e_coul),in(step),                               &
    !dir$ in(numAtomType),in(neigh_alloc),                                        &
    !dir$ inout(virialx,virialy,virialz),                                         &
    !dir$ in(A1),in(A2),in(A3),in(A4),in(A5)                                   
    call lj_cut_coul_long_nonewton_MIC
    call system_clock(T2,clock_rate,clock_max)
    time_nonbond = time_nonbond + real(T2-T1)/real(clock_rate)



    time_pme = time_pme + real(T2-T1)/real(clock_rate)
    
    !========begin of CPU region===================
    !======== begin PME region=====================
    
    call system_clock(T3,clock_rate,clock_max)
    !---convert to PME arrays
    do i = 1,np
       pmex(i) = position(i)%x
       pmey(i) = position(i)%y
       pmez(i) = position(i)%z

    enddo
    call cal_e_pme_kspace(box, box, box, pmex, pmey, pmez, q, e_coul_long, pmefx, pmefy, pmefz,virial)
    virialx_host = virial(1)
    virialy_host = virial(5)
    virialz_host = virial(9)
    !==== the next do loop is commented out currently
    !==== when you call the PME driver, uncomment this out
    do i = 1,np
       ffh(i)%x = ffh(i)%x + pmefx(i)
       ffh(i)%y = ffh(i)%y + pmefy(i)
       ffh(i)%z = ffh(i)%z + pmefz(i)
    end do
    call system_clock(T4,clock_rate,clock_max)
    time_pme = time_pme + real(T4-T3)/real(clock_rate)


    !======= end PME region =================
    call system_clock(T5,clock_rate,clock_max)
    !$omp parallel
    !$omp single

    !$omp task 
    if(numbond.gt.0)then
       call bond_harmonic(step)
       call correct_12_bonds
    endif
    !$omp end task

    
    !$omp task
    if(numangle.gt.0)then
       call angle_harmonic(step)
       call correct_13_bonds
    endif
    !$omp end task
    
    !$omp task
    if(numdihed.gt.0)then
       call dihed_PME_amber
       call correct_14_bonds
    endif
    !$omp end task

    !$omp task
    if(numimpro.gt.0)then
       call impro_amber_force
    endif
    !$omp end task

    !$omp task
    call langevin_debug_intel
    !$omp end task
    
    !$omp end single
    !$omp end parallel
    call system_clock(T6,clock_rate,clock_max)
    time_bond = time_bond + real(T6-T5)/real(clock_rate)
    call system_clock(T6,clock_rate,clock_max)
    time_bond = time_bond + real(T6-T5)/real(clock_rate)

    !=======end of region============================


    !==========merge CPU and MIC
    do i = 1,np
       ff(i)%x = ff(i)%x + ffh(i)%x +ff_bond(i)%x + ff_angle(i)%x + ff_dihed(i)%x+&
            ff_impro(i)%x + ff_lang(i)%x

       ff(i)%y = ff(i)%y + ffh(i)%y +ff_bond(i)%y + ff_angle(i)%y + ff_dihed(i)%y+&
            ff_impro(i)%y + ff_lang(i)%y

       ff(i)%z = ff(i)%z + ffh(i)%z +ff_bond(i)%z + ff_angle(i)%z + ff_dihed(i)%z+&
            ff_impro(i)%z + ff_lang(i)%z
    enddo
    e_coul = e_coul+e_coul_host+e_coul_bond+e_coul_angle+e_coul_dihed

    virialx = virialx + virialx_host + v_bondx + v_anglex + v_dihedx + v_improx
    virialy = virialy + virialy_host + v_bondy + v_angley + v_dihedy + v_improy
    virialz = virialz + virialz_host + v_bondz + v_anglez + v_dihedz + v_improz

    call system_clock(T7,clock_rate,clock_max)
    time_force = time_force + real(T7-T1)/real(clock_rate)

  end subroutine force_wrapper_pme_mic_noasync


  subroutine force_wrapper_pme_mic_noasync_notask(step,neighbor_flag)
    implicit none
    double precision :: virial(9)
    integer :: step
    integer :: neighbor_flag
    integer :: i
    integer :: sig=1
    integer :: tid
    integer :: T1,T2,T3,T4,T5,T6,T7,clock_rate,clock_max
    call zero_values(step)
    
    call system_clock(T1,clock_rate,clock_max)
    !========== begin of mic region ===================
    !dir$ offload  target(mic:ourmic),                                            &
    !dir$ in(position: alloc_if(.false.),free_if(.false.)),                       &
    !dir$ inout(ff: alloc_if(.false.), free_if(.false.)),                         &
    !dir$ nocopy(nlist: alloc_if(.false.) free_if(.false.)),                      &
    !dir$ nocopy(numneigh: alloc_if(.false.) free_if(.false.)),                   &
    !dir$ in(q: alloc_if(.false.) free_if(.false.)),                              &
    !dir$ nocopy(lj1: alloc_if(.false.) free_if(.false.)),                        &
    !dir$ nocopy(lj2: alloc_if(.false.) free_if(.false.)),                        &
    !dir$ nocopy(lj3: alloc_if(.false.) free_if(.false.)),                        &
    !dir$ nocopy(lj4: alloc_if(.false.) free_if(.false.)),                        &
    !dir$ in(rcut2),in(cut_coulsq),in(box),in(ibox),in(np),                       &
    !dir$ in(gewald),in(EWALD_P),in(MY_PIS_INV),in(EWALD_F),in(coulpre)           &
    !dir$ inout(potential), inout(e_coul),in(step),                               &
    !dir$ in(numAtomType),in(neigh_alloc),                                        &
    !dir$ inout(virialx,virialy,virialz),                                         &
    !dir$ in(A1),in(A2),in(A3),in(A4),in(A5)                                   
    call lj_cut_coul_long_nonewton_MIC
    call system_clock(T2,clock_rate,clock_max)
    time_nonbond = time_nonbond + real(T2-T1)/real(clock_rate)

    
    !========begin of CPU region===================
    !======== begin PME region=====================
    
    call system_clock(T3,clock_rate,clock_max)
    !---convert to PME arrays
    do i = 1,np
       pmex(i) = position(i)%x
       pmey(i) = position(i)%y
       pmez(i) = position(i)%z

    enddo
    call cal_e_pme_kspace(box, box, box, pmex, pmey, pmez, q, e_coul_long, pmefx, pmefy, pmefz,virial)
    virialx_host = virial(1)
    virialy_host = virial(5)
    virialz_host = virial(9)
    !==== the next do loop is commented out currently
    !==== when you call the PME driver, uncomment this out
    do i = 1,np
       ffh(i)%x = ffh(i)%x + pmefx(i)
       ffh(i)%y = ffh(i)%y + pmefy(i)
       ffh(i)%z = ffh(i)%z + pmefz(i)
    end do
    call system_clock(T4,clock_rate,clock_max)
    time_pme = time_pme + real(T4-T3)/real(clock_rate)


    !======= end PME region =================
    call system_clock(T5,clock_rate,clock_max)

    if(numbond.gt.0)then
       call bond_harmonic(step)
       call correct_12_bonds
    endif
    
    if(numangle.gt.0)then
       call angle_harmonic(step)
       call correct_13_bonds
    endif
    
    if(numdihed.gt.0)then
       call dihed_PME_amber
       call correct_14_bonds
    endif

    if(numimpro.gt.0)then
       call impro_amber_force
    endif

    call langevin_debug_intel
    

    call system_clock(T6,clock_rate,clock_max)
    time_bond_total = time_bond_total + real(T6-T5)/real(clock_rate)
    !=======end of region============================


    !==========merge CPU and MIC
    do i = 1,np
       ff(i)%x = ff(i)%x + ffh(i)%x +ff_bond(i)%x + ff_angle(i)%x + ff_dihed(i)%x+&
            ff_impro(i)%x + ff_lang(i)%x

       ff(i)%y = ff(i)%y + ffh(i)%y +ff_bond(i)%y + ff_angle(i)%y + ff_dihed(i)%y+&
            ff_impro(i)%y + ff_lang(i)%y

       ff(i)%z = ff(i)%z + ffh(i)%z +ff_bond(i)%z + ff_angle(i)%z + ff_dihed(i)%z+&
            ff_impro(i)%z + ff_lang(i)%z
    enddo
    e_coul = e_coul+e_coul_host+e_coul_bond+e_coul_angle+e_coul_dihed

    virialx = virialx + virialx_host + v_bondx + v_anglex + v_dihedx + v_improx
    virialy = virialy + virialy_host + v_bondy + v_angley + v_dihedy + v_improy
    virialz = virialz + virialz_host + v_bondz + v_anglez + v_dihedz + v_improz

    call system_clock(T7,clock_rate,clock_max)
    time_force = time_force + real(T7-T1)/real(clock_rate)

  end subroutine force_wrapper_pme_mic_noasync_notask











  subroutine force_wrapper_DSF_MIC(step,neighbor_flag)
    implicit none
    integer :: step
    integer :: neighbor_flag
    integer :: i
    integer :: sig=1
    call zero_values(step)


    !========== begin of mic region ===================
    !dir$ offload  target(mic:0) signal(sig),                                     &
    !dir$ in(position: alloc_if(.false.),free_if(.false.)),                       &
    !dir$ inout(ff: alloc_if(.false.), free_if(.false.)),                         &
    !dir$ nocopy(nlist: alloc_if(.false.) free_if(.false.)),                      &
    !dir$ nocopy(numneigh: alloc_if(.false.) free_if(.false.)),                   &
    !dir$ in(q: alloc_if(.false.) free_if(.false.)),                              &
    !dir$ nocopy(lj1: alloc_if(.false.) free_if(.false.)),                        &
    !dir$ nocopy(lj2: alloc_if(.false.) free_if(.false.)),                        &
    !dir$ nocopy(lj3: alloc_if(.false.) free_if(.false.)),                        &
    !dir$ nocopy(lj4: alloc_if(.false.) free_if(.false.)),                        &
    !dir$ in(rcut2),in(cut_coulsq),in(box),in(ibox),in(np),                       &
    !dir$ in(alpha),in(EWALD_P),in(MY_PIS_INV),                                   &
    !dir$ inout(potential), inout(e_coul),in(step),                               &
    !dir$ in(numAtomType),in(neigh_alloc),                                        &
    !dir$ in(A1),in(A2),in(A3),in(A4),in(A5),                                     &
    !dir$ in(e_shift),in(f_shift)              
    call lj_cut_coul_dsf_nonewton_MIC(step)
   


    !========begin of CPU region===================
    if(numbond.gt.0)then
      call bond_harmonic(step)
    endif
 
    
    if(numangle.gt.0)then
       call angle_harmonic(step)
    endif
    if(numdihed.gt.0.and.neighbor_flag.eq.1)then
       call dihedral_opls(step)
       
    elseif(numdihed.gt.0.and.neighbor_flag.eq.0)then
       call dihedral_opls_n2(step)

    elseif(numdihed.gt.0.and.neighbor_flag.eq.3)then
       call dihedral_opls_DSF(step)
    elseif(numdihed.gt.0.and.neighbor_flag.eq.4)then
       call dihedral_opls_n2(step)
    endif
     
    if(numimpro.gt.0)then
       call improper_harmonic_lammps14(step)

    endif

 

    !=======end of region============================

    !dir$ offload_wait target(mic:0) wait(sig)


    !==========merge CPU and MIC
    do i = 1,np
       ff(i)%x = ff(i)%x + ffh(i)%x
       ff(i)%y = ff(i)%y + ffh(i)%y
       ff(i)%z = ff(i)%z + ffh(i)%z
    enddo


  end subroutine force_wrapper_DSF_MIC
  
  
  subroutine force_wrapper_DSF_NOMIC(step,neighbor_flag)
    implicit none
    integer :: step
    integer :: neighbor_flag
    integer :: i
    integer :: sig=1
    call zero_values(step)

      
    call lj_cut_coul_dsf_nonewton_NOMIC(step)
   

    if(numbond.gt.0)then
      call bond_harmonic(step)
    endif
 
    
    if(numangle.gt.0)then
       call angle_harmonic(step)
    endif
    if(numdihed.gt.0.and.neighbor_flag.eq.1)then
       call dihedral_opls(step)
       
    elseif(numdihed.gt.0.and.neighbor_flag.eq.0)then
       call dihedral_opls_n2(step)

    elseif(numdihed.gt.0.and.neighbor_flag.eq.3)then
       call dihedral_opls_DSF(step)
    elseif(numdihed.gt.0.and.neighbor_flag.eq.4)then
       call dihedral_opls_n2(step)
    endif
     
    if(numimpro.gt.0)then
       call improper_harmonic_lammps14(step)

    endif

 
    !==========merge CPU and MIC
    do i = 1,np
       ff(i)%x = ff(i)%x + ffh(i)%x
       ff(i)%y = ff(i)%y + ffh(i)%y
       ff(i)%z = ff(i)%z + ffh(i)%z
    enddo
  end subroutine force_wrapper_DSF_NOMIC


  
  subroutine force_wrapper_PME_NOMIC(step,neighbor_flag)
    implicit none
    double precision :: virial(9)
    integer :: step
    integer :: neighbor_flag
    integer :: i
    integer :: sig=1
    integer :: tid
  
    call zero_values(step)
    

                                  
    call lj_cut_coul_long_nonewton_NOMIC

    !call force_N2_14
    !========begin of CPU region===================
    !======== begin PME region=====================
    

    !---convert to PME arrays
    do i = 1,np
       pmex(i) = position(i)%x
       pmey(i) = position(i)%y
       pmez(i) = position(i)%z

    enddo
    call cal_e_pme_kspace(box, box, box, pmex, pmey, pmez, q, e_coul_long, pmefx, pmefy, pmefz,virial)
    virialx_host = virial(1)
    virialy_host = virial(5)
    virialz_host = virial(9)
   !==== the next do loop is commented out currently
   !==== when you call the PME driver, uncomment this out
    do i = 1,np
       ffh(i)%x = ffh(i)%x + pmefx(i)
       ffh(i)%y = ffh(i)%y + pmefy(i)
       ffh(i)%z = ffh(i)%z + pmefz(i)
    end do
    !======= end PME region =================



    if(numbond.gt.0)then
       call bond_harmonic(step)
       call correct_12_bonds
    endif
    
    if(numangle.gt.0)then
       call angle_harmonic(step)
       call correct_13_bonds
    endif
    
    
    if(numdihed.gt.0)then
       call dihed_PME_amber
       call correct_14_bonds
    endif
    
    if(numimpro.gt.0)then
       call impro_amber_force
    endif
    
    
    !call langevin(dthalf)
    


    !==========merge CPU and MIC
    do i = 1,np
       ff(i)%x = ff(i)%x + ffh(i)%x +ff_bond(i)%x + ff_angle(i)%x + ff_dihed(i)%x+&
            ff_impro(i)%x + ff_lang(i)%x

       ff(i)%y = ff(i)%y + ffh(i)%y +ff_bond(i)%y + ff_angle(i)%y + ff_dihed(i)%y+&
           ff_impro(i)%y + ff_lang(i)%y

       ff(i)%z = ff(i)%z + ffh(i)%z +ff_bond(i)%z + ff_angle(i)%z + ff_dihed(i)%z+&
            ff_impro(i)%z+ff_lang(i)%z
    enddo


    e_coul = e_coul+e_coul_host+e_coul_bond+e_coul_angle+e_coul_dihed
 

    virialx = virialx + virialx_host + v_bondx + v_anglex + v_dihedx + v_improx
    virialy = virialy + virialy_host + v_bondy + v_angley + v_dihedy + v_improy
    virialz = virialz + virialz_host + v_bondz + v_anglez + v_dihedz + v_improz

  end subroutine force_wrapper_PME_NOMIC


  subroutine zero_values(step)
    implicit none
    integer :: step
    integer :: i

    do i =1 ,np
       ff(i)%x = 0.0d0
       ff(i)%y = 0.0d0
       ff(i)%z = 0.0d0

       ff_bond(i)%x = 0.0d0
       ff_bond(i)%y = 0.0d0
       ff_bond(i)%z = 0.0d0


       ff_angle(i)%x = 0.0d0
       ff_angle(i)%y = 0.0d0
       ff_angle(i)%z = 0.0d0

       ff_dihed(i)%x = 0.0d0
       ff_dihed(i)%y = 0.0d0
       ff_dihed(i)%z = 0.0d0

       ff_impro(i)%x = 0.0d0
       ff_impro(i)%y = 0.0d0
       ff_impro(i)%z = 0.0d0

       ffh(i)%x = 0.0d0
       ffh(i)%y = 0.0d0
       ffh(i)%z = 0.0d0

       ff_lang(i)%x = 0.0d0
       ff_lang(i)%y = 0.0d0
       ff_lang(i)%z =0.0d0
    end do

    potential = 0.0d0
    pot_14 = 0.0d0
    e_coul = 0.0d0
    e_coul_14 = 0.0d0
    e_coul_long = 0.0d0
    e_coul_host = 0.0d0
    e_coul_bond =0.0d0
    e_coul_angle=0.0d0
    e_coul_dihed=0.0d0
    e_bond = 0.0d0
    e_angle = 0.0d0
    e_dihedral = 0.0d0
    e_improper = 0.0d0

    virialx = 0.0d0;virialy=0.0d0;virialz=0.0d0
    virialx_host=0.0d0;virialy_host=0.0d0;virialz_host=0.0d0
    v_bondx=0.d0 ; v_bondy=0.0d0 ; v_bondz=0.0d0
    v_anglex=0.d0; v_angley=0.0d0; v_anglez=0.0d0
    v_dihedx=0.d0; v_dihedy=0.0d0; v_dihedz=0.0d0
    v_improx=0.d0; v_improy=0.0d0; v_improz=0.0d0



  end subroutine zero_values



  
  subroutine calculate_pressure
    implicit none
    
    virialx = virialx + pres_corr
    virialy = virialy + pres_corr
    virialz = virialz + pres_corr

    pressure =  nktv2p*( (2.0/3.0)*ke + (1.0/3.0)*(virialx+virialy+virialz))

    pressure = pressure/vol
  

  end subroutine calculate_pressure


  subroutine calculate_ke
    implicit none
    integer :: i
    double precision sum

    sum =0.0d0
    
    do i = 1,np
       sum  = sum + mass(position(i)%type)*(v(i)%x*v(i)%x + &
            v(i)%y*v(i)%y + v(i)%z*v(i)%z)
    enddo
    ke = 0.50d0*sum*amu2e

  end subroutine calculate_ke



 


end module mod_force
