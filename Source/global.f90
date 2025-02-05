module global
  implicit none

!---------------------------------------------------------------------
!*** simulation parameters
!*** read at run time


  !======offload variables
  !dir$ attributes offload:mic::position
  !dir$ attributes offload:mic::ff,ffh
  !dir$ attributes offload:mic::ff_bond
  !dir$ attributes offload:mic::ff_angle
  !dir$ attributes offload:mic::ff_dihed
  !dir$ attributes offload:mic::ff_impro
  !dir$ attributes offload:mic::ff_langevin
  !dir$ attributes offload:mic::q
  !dir$ attributes offload:mic::dthalf
  !dir$ attributes offload:mic::rcut
  !dir$ attributes offload:mic::rcut2,cut_coulsq
  !dir$ attributes offload:mic::box,ibox,boxlei
  !dir$ attributes offload:mic::np
  !dir$ attributes offload:mic::intel_flag,coul_flag
  !dir$ attributes offload:mic::rdf,rdf_delta,num_rdf,numob_rdf
  !dir$ attributes offload:mic::potential,potential_host
  !dir$ attributes offload:mic::pot_14
  !dir$ attributes offload:mic::e_coul,e_coul_host,e_coul_long,e_coul_bond,e_coul_angle,e_coul_dihed
  !dir$ attributes offload:mic::e_coul_14
  !dir$ attributes offload:mic::e_bond
  !dir$ attributes offload:mic:: e_angle
  !dir$ attributes offload:mic:: e_dihedral
  !dir$ attributes offload:mic:: e_improper
  !dir$ attributes offload:mic::virialx,virialy,virialz
  !dir$ attributes offload:mic::virialx_host,virialy_host,virialz_host
  !dir$ attributes offload:mic::time_nonbond
  !dir$ attributes offload:mic::time_single_update
  !dir$ attributes offload:mic::shake_time
  !dir$ attributes offload:mic::lj1,lj2,lj3,lj4
  !dir$ attributes offload:mic::alpha,gewald
  !dir$ attributes offload:mic::EWALD_F,EWALD_P,MY_PIS_INV,MY_PIS
  !dir$ attributes offload:mic::A1,A2,A3,A4,A5
  !dir$ attributes offload:mic::e_shift,f_shift
  !dir$ attributes offload:mic::coulpre
  !dir$ attributes offload:mic::rv2
  !dir$ attributes offload:mic::pptr_cache  
  !dir$ attributes offload:mic::cache
  !dir$ attributes offload:mic::xcache,ycache,zcache
  !dir$ attributes offload:mic::dr2array
  !dir$ attributes offload:mic::dr2AOS
  !dir$ attributes offload:mic::atombin
  !dir$ attributes offload:mic::nlist
  !dir$ attributes offload:mic::numneigh
  !dir$ attributes offload:mic::cache_size
  !dir$ attributes offload:mic::neigh_alloc,nlistsize
  !dir$ attributes offload:mic::start,endposit
  !dir$ attributes offload:mic::cnum
  !dir$ attributes offload:mic::nncache
  !dir$ attributes offload:mic::rn
  !dir$ attributes offload:mic::ncellT

   !---file IDs
  integer :: pressure_f_id, box_f_id, vol_f_id, clus_f_id,PE_f_id,total_energy_f_id,KE_f_id,temp_f_id


  
  !-----------------------------------------------------------------------
  !*** DSF variables and PME VARIABLES
  !***
  double precision :: alpha,gewald
  double precision :: EWALD_F,EWALD_P,MY_PIS,MY_PIS_INV
  double precision :: A1,A2,A3,A4,A5
  double precision :: e_shift,f_shift
  real*4  :: coulpre
  integer :: ngrid_x,ngrid_y,ngrid_z
  real*4,allocatable :: pmex(:),pmey(:),pmez(:)
  real*4,allocatable :: pmefx(:),pmefy(:),pmefz(:)




  !dir$ attributes align:64::ainfo
  type ainfo
     real*4  :: x,y,z
     integer :: type
  end type ainfo

  type coords
     real*4 :: x,y,z
  end type coords

  !---pme variables
  

  !--- nondimensional groups
  double precision :: reps,rmass,rsig,rdeltat,rtemp,rq

  integer :: restart_mult

  type(ainfo),allocatable :: position(:)
  type(ainfo),allocatable :: ps(:)
  type(ainfo),allocatable :: p0(:)
  type(coords),allocatable :: ff(:),ffh(:),ff_bond(:),ff_angle(:),ff_dihed(:),ff_impro(:),ff_lang(:)
  type(coords),allocatable :: ff0(:)

  type(coords),allocatable :: v(:),vs(:),v0(:)
  real*4, allocatable :: drx(:),dry(:),drz(:)
  real*4, allocatable :: q(:)
  double precision, allocatable :: qs(:)
  double precision, allocatable :: q0(:)
  integer, allocatable :: integrate_flag(:),integrate_flag_ss(:)
  integer :: ncell_alloc


  !----total mass of system
  double precision :: tot_mass


  !---12 13 14 interaction corrections
  double precision :: fudge_12, fudge_13, fudge_14,fudge_14_coul


  !---amu to kcal/mol
  double precision :: amu2e
  !---amu to pressure
  double precision :: nktv2p

  integer, allocatable :: global_id(:),global_id_ss(:),global_id0(:)
  real*4 :: dt,dthalf,dthalf2,dt2,kcal2amu
  real*4 :: rcut,neighbor_pad,triggersq
  double precision :: rcut2, cut_coulsq
  real*4          :: cut_coul
 
  real*4 :: kboltz
  real*4 :: vol,hbox
  real*4 :: box,ibox
  real*4 :: dens
  real*4 :: temp,anneal_temp,P,beta,STDtemp,nu
  real*4 :: kn
  real*4 :: timeclock
  real*4 :: ecorr, pcorr,ecut,ecutAvr
  real*4 :: totKE,totPE,totE
  integer :: np
  integer :: newton3rd_flag,intel_flag,coul_flag
  integer :: istart
  integer :: mdeq,mdsteps
  integer :: nsample, nadjust
  integer :: numSweep,numStep,numSweepEquil
  character(len=3) :: ensemble
  character(len=3) :: STA



  !----AMBER DIHEDRAL VARIABLES
  double precision,allocatable :: k_dihed(:),shift_dihed(:),cos_shift(:),sin_shift(:),multiplicity(:)


  !---AMBER IMPROPER VARIABLES
  double precision,allocatable :: k_impro(:)
  integer,allocatable :: multiplicity_impro(:),sign_impro(:)
  
   !---------------------------------------------------------------------
  !*** langevin thermostat variables
  !***
  real*4 :: gconst
  real*4 :: fric,muT,gammaT,hT
  real*4 :: viscous_damp
  real*4,allocatable :: r1(:),r2(:),r3(:)

  !--------
  !*** monte carlo 
  integer :: naccept,nreject
  integer :: ar_vol
  integer :: acc_vol
  integer :: att_vol
  integer :: cc
  real*4 :: vmax


  !--------------
  !***  barostat variables
  real*4 :: vdw_corr,pres_corr
  real*4 :: pdamp,ptarget
  real*4 :: vol_vel,gammaP
  real*4 :: rconst
  real*4 :: piston_mass,inv_piston_mass
  real*4 :: xi, eta,muP,hP,tau
  integer,allocatable :: com_part(:,:),com_part0(:,:)

 !---msd variables
  integer :: msd_record

  !---radial distribut
  real*4, allocatable :: rdf(:)
  real*4 :: rdf_delta
  integer :: num_rdf,numob_rdf

  !---vhisto
  integer :: num_vdist
  real*4 :: vdist_delta
  real*4,allocatable :: vhisto(:)
  


  !---density histogramcs
  integer :: num_dens_ob
  integer,allocatable :: histogram_dens(:)


  !---triclinic box
  real*4 :: xlo,xhi,ylo,yhi,zlo,zhi
  real*4 :: xprd,yprd,zprd,xy,xz,yz
  
  
  !--------------
  !*** tip4p
  integer :: num_tip4p_pairs
  integer,allocatable :: tip4p_list(:)
  type(coords),allocatable :: nopbc(:)
  
  
  !-----------------------------------------------------------------------
  !*** measured quantities
  !*** 
  real*4 :: avgE,varE,avgP,varP,avgV,varV,avgRHO,varRHO,avgT,varT

  double precision :: e_self
  double precision :: pressure,ke,tot_e
  double precision :: potential,potential_host
  double precision :: pot_14
  double precision :: e_coul,e_coul_host,e_coul_long,e_coul_bond,e_coul_angle,e_coul_dihed
  double precision :: e_coul_14
  double precision :: e_bond
  double precision :: e_angle
  double precision :: e_dihedral
  double precision :: e_improper
  double precision :: virialx,virialy,virialz
  double precision :: virialx_host,virialy_host,virialz_host
  double precision :: v_bondx,v_bondy,v_bondz
  double precision :: v_anglex,v_angley,v_anglez
  double precision :: v_dihedx,v_dihedy,v_dihedz
  double precision :: v_improx,v_improy,v_improz
  double precision :: virial(9)


  !==== record timing data from run
  !------------------------------------------------------------------
  double precision :: time_bond, time_nonbond, time_angle,time_dihed,time_data
  double precision :: time_impro,time_bond_total
  double precision :: time_neigh
  double precision :: time_video  
  double precision :: time_pme
  double precision :: time_single_update
  double precision :: shake_time
  double precision :: time_initialize
  double precision :: time_integrate
  double precision :: time_cache
  double precision :: time_sort
  double precision :: time_lan
  double precision :: time_ver
  double precision :: time_clus
  double precision :: time_force
  double precision :: time_vol
  double precision :: time_12corr,time_13corr,time_14corr
  integer :: nn_builds
  


 

  !-----------------------------------------------------------------------
  !*** LL and Vlist 
  !*** 
  double precision:: rv,rv2,rvrc,rvrc2
  real*4, allocatable :: min_x(:),min_y(:),min_z(:)
  real*4, allocatable :: drxs(:),drys(:),drzs(:)
  integer, allocatable :: pptr_cache(:)
  real*4, allocatable :: cache(:)
  real*4, allocatable :: xcache(:),ycache(:),zcache(:)
  real*4, allocatable :: dr2array(:)
  type dr2
     real*4  :: dr2
     integer :: atom
  end type dr2
  type(dr2),allocatable :: dr2AOS(:)
  integer, allocatable :: atombin(:)
  integer,allocatable :: ll(:),HOC(:)
  integer,allocatable :: vlist(:,:)
  integer,allocatable :: nlist(:)
  integer,allocatable :: nn(:)
  integer,allocatable :: numneigh(:)
  integer,allocatable :: nnfirst(:), nnlast(:)
  integer  :: cache_size 
  integer :: maxneigh
  integer :: neigh_alloc,nlistsize
  integer :: vlist_flag
  integer,allocatable :: start(:),endposit(:)
  integer,allocatable :: cnum(:)
  integer,allocatable :: nncache(:)
  !-----------------------------------------------------------------------
  !*** position/list arrays and list parameters
  !*** 
  
  type atom_cell
     integer :: cell
     integer :: loc
  end type atom_cell

  type cell_atom_data
     type(ainfo),allocatable :: loc_r(:)
  end type cell_atom_data

  type cell_data
     integer :: num
     integer :: size
  end type cell_data

  type minim
     real*4 :: min_z,min_y,min_x
     integer :: cnum
  end type minim

  !dir$ attributes offload:mic::rn
  real*4 :: rn
  !dir$ attributes offload:mic::ncellT
  integer :: ncellT
  integer :: ncellD
  integer :: alloc_size
  
  !----------------------------------------------------------------------
  !***
  !*** sorting data

  integer, allocatable :: tt(:)
  integer, allocatable  :: ttb(:)

  type sort_data_bond
     integer :: ba,bi
  end type sort_data_bond

  type SOA_sort
     type(sort_data_bond), allocatable :: ptr(:)
     integer :: num
  end type SOA_sort
  
  type(SOA_sort),allocatable :: spec_ptr(:)
  type(SOA_sort),allocatable :: spec_ptr_ss(:)

  !----------------------------------------------------------------------
  !***
  !*** shake/rattle variables
  type shake_info
     double precision :: bond,bond2,bond3
     integer :: atom1,atom2,atom3,atom4
     integer :: num
  end type shake_info

  type xshake_data
     double precision :: x,y,z
  end type xshake_data

  !dir$ attributes offload:mic::nshake
  integer :: nshake, maxiter,ncons
  double precision :: tolerance
  type(xshake_data), allocatable :: xshake(:)
  type(shake_info),allocatable :: frac_shake(:)
  type(shake_info),allocatable :: shakeatom(:),shakeatom0(:)
  integer, allocatable :: nummolshake(:)
  



  !----***NPT swap arrays
  type(ainfo),allocatable :: p0p(:)
  type(coords),allocatable :: v0p(:)
  real*4,allocatable :: q0p(:)
  integer,allocatable :: num1bond0p(:),num2bond0p(:),num3bond0p(:)
  integer,allocatable :: bondlist0p(:,:),anglelist0p(:,:),dihedlist0p(:,:)
  integer,allocatable :: improlist0p(:,:)
  integer,allocatable :: specbond0p(:,:)
  integer,allocatable :: global_id0p(:)
  integer,allocatable :: mol0p(:,:)
  type(shake_info),allocatable :: shakeatom0p(:)









  
  !-----------------------------------------------------------------------
  !***
  !*** variables related to molecular geometries and types
  !*** numMol is number of molecules present in simulation
  !*** numMoltype is number of different molecular types present
  !*** numAtomType is number of different atom types
  !*** numbond is total number of 1-2 bonds
  !*** numbondtype is number of 1-2 bond types
  !*** numangle is total number of 1-3 bonds
  !*** numangletype is total number of 1-3 bond types
  !*** type contains type for each atom
  !*** numMolArray contains number of molecules for each molecular type
  !*** numMolAtom contains number of atoms in each molecular type
  !*** num1bond, num2bond, num3bond, are number of bonds for each atom. THIS IS CUMMULATIVE
  !*** specbond contains the global tags that each atom is bonded to
  !*** bondlist contains for each 12 bond the atoms present (newtons 3rd i.lt.j) and the bond type
  !*** anglelist contains for each 12 bond the atoms present (newtons 3rd i.lt.j) and the bond type
  !*** atom2bond returns bond type given atom types input
  !*** bondcoef contains k and r0 value for each 1-2 bond type
  !*** anglecoeff contains k and r0 value for each 1-2 bond type
  !*** Mol2AtomType contains atom types present for each molecular type
  !*** FC array contains fractional input coordinates for each molecular type
  !*** Mol2NumBonds contains for each molecular type,numbe of 12,13,14 bonds for each contained atom
  !*** Mol2SpecInfo contains, for each atom in each molecular type, the atoms to which it interacts 12 13 14


  integer :: numMol
  integer :: numMolType
  !dir$ attributes offload:mic::numAtomType
  integer :: numAtomType
  !dir$ attributes offload:mic:: numbond
  integer :: numbond
  integer :: numbondtype
  !dir$ attributes offload:mic:: numangle
  integer :: numangle
  integer :: numangletype
  !dir$ attributes offload:mic:: numdihed
  integer :: numdihed
  integer :: numdihedtype
  !dir$ attributes offload:mic:: numimpro
  integer :: numimpro
  integer :: numimprotype
  character(len=40),allocatable :: type2string(:)
  integer,allocatable :: ltype2ftype(:)
  !dir$ attributes offload:mic::type
  integer,allocatable :: type(:)
  integer,allocatable :: numMolArray(:)
  integer,allocatable :: numMolAtom(:)
  !dir$ attributes offload:mic::num3bond,num1bond,num2bond
  integer,allocatable :: num1bond(:),num2bond(:),num3bond(:)
  integer,allocatable :: num1bondss(:),num2bondss(:),num3bondss(:)
  integer,allocatable :: num1bond0(:), num2bond0(:),num3bond0(:)
  !dir$ attributes offload:mic::specbond
  integer,allocatable :: specbond(:,:),specbond_ss(:,:),specbond0(:,:)
  integer,allocatable :: bondlist(:,:),bondlist0(:,:)
  integer,allocatable :: anglelist(:,:), anglelist0(:,:)
  integer,allocatable :: dihedlist(:,:), dihedlist0(:,:)
  integer,allocatable :: improlist(:,:), improlist0(:,:)
  integer,allocatable :: atom2bond(:,:)
  double precision,allocatable :: bondcoeff(:,:)
  double precision,allocatable :: anglecoeff(:,:)
  double precision,allocatable :: dihedcoeff(:,:)
  double precision,allocatable :: improcoeff(:,:)
  double precision, allocatable :: qtype(:)
  double precision, allocatable :: eps(:)
  !dir$ attributes offload:mic::sig
  double precision, allocatable :: sig(:)
  double precision, allocatable :: type1(:)
  double precision, allocatable :: type2(:)
  real*4, allocatable :: lj1(:),lj2(:),lj3(:),lj4(:)
  double precision, allocatable :: lj14_1(:,:),lj14_2(:,:)
  double precision, allocatable :: lj14_3(:,:),lj14_4(:,:)
  double precision, allocatable :: mass(:)

  type MolAtomType
     integer,allocatable :: types(:)
  end type MolAtomType

  type(MolAtomType),allocatable :: Mol2AtomType(:)


  type CollecBondAtomTypes
     integer :: tag1,tag2
  end type CollecBondAtomTypes
  type(CollecBondAtomTypes),allocatable :: Bond2atypes(:)

  type CollecAngleAtomTypes
     integer :: tag1,tag2,tag3
  end type CollecAngleAtomTypes
  type(CollecAngleAtomTypes),allocatable :: angle2atypes(:)



  type molinfo
     double precision :: x,y,z
     integer :: type
  end type molinfo

  type fracCoords
     type(molinfo),allocatable :: coords(:)
  end type fracCoords

  type(fracCoords),allocatable :: FC(:)


  type Mol2BondNumType
     integer,allocatable :: numbonds(:,:)
  end type Mol2BondNumType

  type(Mol2BondNumType),allocatable :: Mol2NumBonds(:) 


  type spec
     integer :: atom,atom2,atom3,type
  end type spec

  type Mol2MolSpec
     type(spec),allocatable :: spec12(:,:)
     type(spec),allocatable :: spec13(:,:)
     type(spec),allocatable :: spec14(:,:)
  end type Mol2MolSpec

  type(Mol2MolSpec),allocatable :: Mol2MolSpecInfo(:)

  type Mol2MolSpecSum
     type(spec),allocatable :: specbond(:,:)
  end type Mol2MolSpecSum


  type(Mol2MolSpecSum),allocatable :: Mol2SpecTotal(:)

  !-----------------------------------------------------------------------
  !*** clustering algorithm variables. Streinhardt order parameters
  !*** 
  !*** center = positin in nuchisto array where n0 is stored
  type BOP
     double precision :: BOPcrit
     double precision :: rn
     double precision :: rcut
     double precision :: rcut2
     integer :: minfrenkel
     integer :: ncellT
     integer :: ncellD
     integer :: alloc_size
     integer :: center
  end type BOP

  type(BOP) :: BOPinfo
  type(atom_cell),allocatable :: mol_tr_cl(:)
  double precision,allocatable :: q4(:)
  double precision,allocatable :: q6(:)
  double precision,allocatable :: rQ4(:,:)
  double precision,allocatable :: rQ6(:,:)
  double precision,allocatable :: iQ6(:,:)
  double precision,allocatable :: iQ4(:,:)
  double precision,allocatable :: min_x_cl(:)
  double precision,allocatable :: min_y_cl(:)
  double precision,allocatable :: min_z_cl(:)
  double precision,allocatable :: min_x2(:)
  double precision,allocatable :: min_y2(:)
  double precision,allocatable :: min_z2(:)
  double precision,allocatable :: x_cl(:)
  double precision,allocatable :: y_cl(:)
  double precision,allocatable :: z_cl(:)

  integer,allocatable :: frenkelnumber(:)
  integer,allocatable :: numNeigh_cl(:)
  !integer,allocatable :: nlist(:,:)
  integer,allocatable :: clusterhisto(:)
  integer,allocatable :: tempnum_cl(:)
  integer,allocatable :: start_cl(:)
  integer,allocatable :: end_cl(:)
  integer,allocatable :: cnum2(:)
  integer,allocatable :: cnum_cl(:)
  integer,allocatable :: nuchisto(:)
  


  !---template matching order variables
  !*** 
  !*** 
  !*** ppta is array containing coordinates of atoms in TTA molecules
  !*** mol is an array of length np that contains what molecular type atom belongs to and ptr to ptta if needed
  !*** xyz(3,0:numneigh) is array to which we are comparing against template (tagged molecule + enviroment coordinates)
  !*** xyz_alloc is current allocation size of xyz, xyz_number is number in array
  !*** temp array will store template coordinates
  !*** numatomtemp is total number of atoms in template file
  !*** numoltemp   is total number of molecules in template file
  !*** numatompermol is number of atoms in one molecule in template
  !*** scale is the oxygen carbon bond length ( 0.615A)
  !*** CalcCluster is how many steps till we cluster
  !*** densmin is the number of minimum neighbors to be considered solid


  type xyz_type
     double precision :: x,y,z
  end type xyz_type

  type tempz
     double precision, allocatable :: temp(:,:)
  end type tempz

  type template_num
     integer :: numatom
     integer :: nummol
  end type template_num


  type(xyz_type), allocatable  :: ptta(:),ptta0(:),ptta_grof(:),ptta_full(:)
  type(xyz_type), allocatable  :: pruned(:)
  type(xyz_type), allocatable  :: com(:),scom(:)
  type(xyz_type), allocatable  :: shift(:)
  type(tempz), allocatable :: tempa(:)
  type(tempz), allocatable :: tempa_noprune(:)
  type(template_num), allocatable :: numtemp(:)
  double precision, allocatable :: dihed(:)
  double precision, allocatable :: xyzcom2match(:,:) 
  double precision, allocatable :: xyz(:,:),xyzcom(:,:),xyzf(:,:)
  double precision, allocatable :: pbcshift(:,:)
  double precision, allocatable :: tempf(:,:)
  double precision, allocatable :: comtemp(:,:)
  type(xyz_type) :: BOPtemp_prune(4)

  double precision :: scale
  double precision :: rmsdcut
  integer, allocatable :: mol(:,:),mol0(:,:)
  integer, allocatable :: molsort(:,:)
  integer, allocatable :: mvl(:)
  integer, allocatable :: mll(:),mHOC(:)
  integer, allocatable :: densvl(:),densnn(:)
  integer, allocatable :: mstart(:),mendposit(:),mcnum(:)
  integer, allocatable :: mnum(:)
  integer, allocatable :: mttb(:)
  integer, allocatable :: NumPolyTemp(:)
  integer, allocatable :: sortcom(:)
  integer, allocatable :: denstrack(:)
  integer,allocatable :: denshisto(:)
  integer, allocatable :: belongsto(:)
  integer, allocatable :: surfaceClus(:)
  integer, allocatable :: minrmsd(:)

 
  double precision :: mrcut,mrcut2,mrn
  double precision :: mdens
  double precision :: CObond
  double precision :: denscut,denscut2

  integer :: cluster_type
  integer :: comnumber
  integer :: mncellT,mncellD
  integer :: mol_alloc
  integer :: numatomtemp, nummoltemp,numatompertemp
  integer :: xyz_alloc,xyz_number,xyzcom_alloc
  integer :: numPolyMorph, numTemPlate,numatns
  integer :: densmin
  integer :: densvlalloc
  integer :: uniquepairs
  integer :: mncell_alloc


  !-----------------------------------------------------------------------
  !*** MPI VARIABLES
  !*** 
  !*** 
  integer,allocatable :: p2n(:)
  integer,allocatable :: p2w(:)
  !dir$ attributes offload:mic::size,ourmic,rank
  integer :: size,rank,ourmic
  character(len=150) :: rank_char,restart_files_path,restart_files_path_final
  character(len=150) :: grofile_path, data_files_path
  character(len=150) :: sim_string
  character(len=150) :: path_string


  !======================================================================
  !*** FILENAME VARIABLES
  !***
  character(len=500) :: filename_input
  character(len=500) :: filename_mol
  character(len=500) :: filename_topol
  character(len=500) :: filename_lammps
  character(len=5000) :: filename_gro
  character(len=5000) :: filename_restart
  integer :: file_type_flag
end module global
