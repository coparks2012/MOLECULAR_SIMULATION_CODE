module mod_memory

  use global
  contains

    subroutine memory
      implicit none
      integer :: totalchange,i
      integer :: nummol,numat
      integer :: alloc
      integer :: mol_max
      

      !---lets perform some initial allocation
      allocate(drx(np),dry(np),drz(np))
      if(nshake.gt.0)then
         allocate(xshake(np))
      endif
      allocate(r1(np),r2(np),r3(np))
      allocate(ff(np),ff0(np),ffh(np))
      allocate(ff_bond(np))
      allocate(ff_angle(np))
      allocate(ff_dihed(np))
      allocate(ff_impro(np))
      allocate(ff_lang(np))
      allocate(pmex(np),pmey(np),pmez(np))
      allocate(pmefx(np),pmefy(np),pmefz(np))
      allocate(atombin(np))
      allocate(v(np),v0(np),vs(np))
      allocate(position(np),p0(np))
      allocate(q(np),qs(np),q0(np))
      allocate(ps(np))
      allocate(specbond(12,np),specbond_ss(12,np),specbond0(12,np))
      allocate(num1bond(np),num2bond(np),num3bond(np))
      allocate(num1bondss(np),num2bondss(np),num3bondss(np))
      allocate(num1bond0(np),num2bond0(np),num3bond0(np))

      !----************NPT swap arrays
      allocate(p0p(np))
      allocate(v0p(np))
      allocate(global_id0p(np))
      allocate(mol0p(2,np))
      allocate(q0p(np))
      allocate(num1bond0p(np),num2bond0p(np),num3bond0p(np))
    
      


      allocate(tt(np))
      allocate(ttb(np))
      allocate(mol(2,np),mol0(2,np),molsort(2,np))
      allocate(global_id(np),global_id_ss(np),global_id0(np))
      allocate(integrate_flag(np))
      allocate(integrate_flag_ss(np))
      

      nummol = numMolArray(1)
      numat  = numMolAtom(1)
      densvlalloc = 12*nummol
      allocate(com(nummol))
      allocate(scom(nummol))
      allocate(mttb(numMolArray(1)))
      allocate(pbcshift(3,nummol))
      allocate(xyzf(3,nummol),tempf(3,numatomtemp))
      allocate(pruned(numatompertemp*nummol))
      allocate(denstrack(nummol))
      allocate(minrmsd(nummol))
      allocate(belongsto(nummol))
      allocate(surfaceclus(nummol))
      allocate(densvl(densvlalloc*nummol))
      allocate(densnn(nummol))
      !allocate(ptta(nummol*numatns),ptta0(nummol*numatns),ptta_grof(nummol*numatns))
      allocate(ptta(nummol*numatns))
      allocate(ptta_full(nummol*numat))
      allocate(ptta_grof(nummol*numat))
      allocate(clusterhisto(nummolarray(1)))
      allocate(denshisto(nummolarray(1)))
      allocate(dihed(nummolarray(1)))


      !---dens hist variables
      allocate(histogram_dens(0:30))
    




      !---allocate clustering neighbor arrays
      !mol_alloc    = numMolArray(1)
      !xyz_alloc    = mol_alloc*numatompertemp
      !xyzcom_alloc = mol_alloc
      !allocate( mvl(mol_alloc*numMolArray(1)   ) )
      !allocate( shift(mol_alloc*numMolArray(1)) )
      !allocate( mnum( numMolArray(1)) )
      !allocate( xyz(3, xyz_alloc))
      !allocate( xyzcom(3, xyzcom_alloc))
      !allocate( xyzcom2match(2,xyzcom_alloc))
      !allocate( sortcom(xyzcom_alloc))
      !allocate( mll(numMolArray(1)))
      !allocate( mHOC(0:mncellT-1))
      !allocate( mstart(0:mncellT-1))
      !allocate( mendposit(0:mncellT-1))
      !allocate( mcnum(0:26*mncellT))
      mol_alloc    = numMolArray(1)
      xyz_alloc    = mol_alloc*numatompertemp
      xyzcom_alloc = mol_alloc
      allocate( mvl(mol_alloc*numMolArray(1)   ) )
      allocate( shift(mol_alloc*numMolArray(1)) )
      allocate( mnum( numMolArray(1)) )
      allocate( xyz(3, xyz_alloc))
      allocate( xyzcom(3, xyzcom_alloc))
      allocate( xyzcom2match(2,xyzcom_alloc))
      allocate( sortcom(xyzcom_alloc))
      allocate( mll(numMolArray(1)))
      allocate( mHOC(0:mncell_alloc))
      allocate( mstart(0:mncell_alloc))
      allocate( mendposit(0:mncell_alloc))
      allocate( mcnum(0:26*mncell_alloc))
      

      !---allocate verlet lists 
      if(vlist_flag.ne.0)then
         allocate(ll(np))
         allocate(HOC(0:ncell_alloc))
         allocate(start(0:ncell_alloc))
         allocate(endposit(0:ncell_alloc))


         if(newton3rd_flag.eq.1)then
            allocate(cnum(0:13*ncell_alloc))
         elseif(newton3rd_flag.eq.0)then
            
            allocate(cnum(0:27*ncell_alloc))
         endif
         allocate(numneigh(np))

         neigh_alloc = 512
         nlistsize   = neigh_alloc*np

         cache_size = 2048    

         if(intel_flag.eq.0)then
       !     allocate(nlist(nlistsize))
       !     allocate(pptr_cache(ncellT*cache_size))
       !     allocate(cache(3*cache_size*ncellT))
       !     allocate(xcache(cache_size*ncellT))
       !     allocate(ycache(cache_size*ncellT))
       !     allocate(zcache(cache_size*ncellT))
       !     allocate(dr2array(np*cache_size))
       !     allocate(nncache(0:ncellT*cache_size))


            allocate(nlist(nlistsize))
            allocate(dr2array(np*cache_size))
         endif

         if(intel_flag.eq.1)then
            call xeon_phi_memory_neighbor_alloc_full
            !call xeon_phi_memory_backup
         endif


      endif


    end subroutine memory
    


   


    subroutine memory_santiso
      implicit none
      integer :: totalchange,i
      integer :: nummol,numat
      integer :: alloc
      integer :: mol_max
      
      !---lets perform some initial allocation

      allocate(position(np),p0(np))
      allocate(ttb(np))
      allocate(mol(2,np))
      

      nummol = numMolArray(1)
      numat  = numMolAtom(1)
      densvlalloc = 12*nummol
      allocate(com(nummol))
      allocate(scom(nummol))
      allocate(mttb(numMolArray(1)))
      allocate(pbcshift(3,nummol))
      allocate(xyzf(3,nummol),tempf(3,numatomtemp))
      allocate(pruned(numatompertemp*nummol))
      allocate(denstrack(nummol))
      allocate(minrmsd(nummol))
      allocate(belongsto(nummol))
      allocate(densvl(densvlalloc*nummol))
      allocate(densnn(nummol))
      allocate(ptta(nummol*numatns),ptta0(nummol*numatns),ptta_grof(nummol*numatns))
      allocate(clusterhisto(nummolarray(1)))
      allocate(denshisto(nummolarray(1)))
      allocate(dihed(nummolarray(1)))


      !---dens hist variables
      allocate(histogram_dens(0:30))
  
      mol_alloc    = numMolArray(1)
      xyz_alloc    = mol_alloc*numatompertemp
      xyzcom_alloc = mol_alloc
      allocate( mvl(mol_alloc*numMolArray(1)   ) )
      allocate( shift(mol_alloc*numMolArray(1)) )
      allocate( mnum( numMolArray(1)) )
      allocate( xyz(3, xyz_alloc))
      allocate( xyzcom(3, xyzcom_alloc))
      allocate( xyzcom2match(2,xyzcom_alloc))
      allocate( sortcom(xyzcom_alloc))
      allocate( mll(numMolArray(1)))
      allocate( mHOC(0:mncell_alloc))
      allocate( mstart(0:mncell_alloc))
      allocate( mendposit(0:mncell_alloc))
      allocate( mcnum(0:26*mncell_alloc))
      

 
 

    end subroutine memory_santiso
    




    subroutine xeon_phi_memory_neighbor_alloc_full
      implicit none

      !dir$ offload begin target(mic:ourmic) in(nlistsize,cache_size,np,ncellT),    &
      !dir$ nocopy(nlist:length(0) alloc_if(.false.) free_if(.false.)),        &
      !dir$ nocopy(numneigh:length(0) alloc_if(.false.) free_if(.false.)),     &
      !dir$ nocopy(dr2array:length(0) alloc_if(.false.) free_if(.false.))

      write(1234,*)'hello from mic'
      allocate(nlist(nlistsize))
      allocate(numneigh(np))
      allocate(dr2array(np*cache_size))

   
     ! write(1234,*)'nlist',allocated(nlist)
     ! write(1234,*)'numneigh',allocated(numneigh)
     ! write(1234,*)'dr2array',allocated(dr2array(np*cache_size))

     ! print*,'mic memory from rank',rank
     ! print*,'nlist',allocated(nlist)
     ! print*,'numneigh',allocated(numneigh)
     ! print*,'dr2array',allocated(dr2array(np*cache_size))
     ! print*
    
      !dir$ end offload



     
      print*,'we are now allocating mic memory'

      !dir$ offload_transfer target(mic:ourmic) in(position: alloc_if(.true.) free_if(.false.)), &
      !dir$ in(specbond: alloc_if(.true.) free_if(.false.)),                                &
      !dir$ in(num3bond: alloc_if(.true.) free_if(.false.)),                                &
      !dir$ in(atombin: alloc_if(.true.) free_if(.false.)),                                 &
      !dir$ in(start: alloc_if(.true.) free_if(.false.)),                                   &
      !dir$ in(endposit: alloc_if(.true.) free_if(.false.)),                                &
      !dir$ in(cnum: alloc_if(.true.) free_if(.false.)),                                    &           
      !dir$ in(ff: alloc_if(.true.), free_if(.false.)),                                     &
      !dir$ in(q: alloc_if(.true.) free_if(.false.)),                                       &
      !dir$ in(lj1: alloc_if(.true.) free_if(.false.)),                                     &
      !dir$ in(lj2: alloc_if(.true.) free_if(.false.)),                                     &
      !dir$ in(lj3: alloc_if(.true.) free_if(.false.)),                                     &
      !dir$ in(lj4: alloc_if(.true.) free_if(.false.))                                  
    

 
      print*,'we managed to allocate that'


    end subroutine xeon_phi_memory_neighbor_alloc_full



    
    subroutine moltype_memory
      implicit none
      integer :: i,j
      allocate(FC(numMolType))
      allocate(numMolArray(numMolType))
      allocate(numMolAtom(numMolType))
      allocate(Mol2NumBonds(numMolType))
      allocate(Mol2SpecTotal(numMolType))
      allocate(nummolshake(nummoltype))
            
      
    end subroutine moltype_memory


    subroutine moltype_datatype_memory(i)
      integer :: i

      allocate(FC(i)%coords(numMolAtom(i)))
      allocate(Mol2NumBonds(i)%numbonds(4,numMolAtom(i)))
      allocate(Mol2specTotal(i)%specbond(16,numMolAtom(i)))

    end subroutine moltype_datatype_memory

    
    subroutine atomtype_memory
      implicit none
      !---allocate arrays that are atom type look ups
      allocate(eps(numAtomType),sig(numAtomType),qtype(numAtomType))
      allocate(ltype2ftype(numAtomType))
      allocate(type1(numAtomType*numAtomType))
      allocate(type2(numAtomType*numAtomType))
      allocate(mass(numAtomType))
      allocate(type2string(numAtomType))
      allocate(lj1(numAtomType*numAtomType))
      allocate(lj2(numAtomType*numAtomType))
      allocate(lj3(numAtomType*numAtomType))
      allocate(lj4(numAtomType*numAtomType))
      allocate(lj14_1(numAtomType,numAtomType))
      allocate(lj14_2(numAtomType,numAtomType))
      allocate(lj14_3(numAtomType,numAtomType))
      allocate(lj14_4(numAtomType,numAtomType))
      
    end subroutine atomtype_memory
    
    
    subroutine bondtype_memory
      implicit none
      
      allocate(bondcoeff(2,numbondtype))
      allocate(Bond2atypes(numbondtype))
      
    end subroutine bondtype_memory
    
    
    subroutine angletype_memory
      implicit none
      
      allocate(anglecoeff(4,numangletype))
      allocate(angle2atypes(numangletype))
    end subroutine angletype_memory
    
    
    subroutine dihedtype_memory
      implicit none
      allocate(dihedcoeff(5,numdihedtype))

      !---new dihedral arrays
      allocate(shift_dihed(numdihedtype))
      allocate(cos_shift(numdihedtype))
      allocate(sin_shift(numdihedtype))
      allocate(multiplicity(numdihedtype))
      allocate(k_dihed(numdihedtype))
    end subroutine dihedtype_memory

    subroutine improdihedtype_memory
      implicit none
      allocate(improcoeff(2,numimprotype))

      allocate(k_impro(numimprotype))
      allocate(multiplicity_impro(numimprotype))
      allocate(sign_impro(numimprotype))
    end subroutine improdihedtype_memory



    subroutine totalbond_memory
      implicit none
      allocate(bondlist(3,numbond),bondlist0(3,numbond))

      allocate(bondlist0p(3,numbond))

    end subroutine totalbond_memory


    subroutine totalangle_memory
      implicit none
      allocate(anglelist(4,numangle),anglelist0(4,numangle))

      allocate(anglelist0p(4,numangle))

    end subroutine totalangle_memory

    subroutine totaldihed_memory
      implicit none
      allocate(dihedlist(5,numdihed),dihedlist0(5,numdihed))
      
      allocate(dihedlist0p(5,numdihed))
    end subroutine totaldihed_memory


    subroutine improdihed_memory
      implicit none
      allocate(improlist(5,numimpro),improlist0(5,numimpro))
      allocate(improlist0p(5,numimpro))

      

    end subroutine improdihed_memory


    subroutine total_shake_memory
      implicit none

      print*,'WE ARE ALLOCATING SHAKE ARRAYS'
      print*,'SHAKE WILL BE USED WITH A TIMSTEP:',dt

      allocate(shakeatom(nshake),shakeatom0(nshake))

      allocate(shakeatom0p(nshake))
    end subroutine total_shake_memory

  end module mod_memory
  
