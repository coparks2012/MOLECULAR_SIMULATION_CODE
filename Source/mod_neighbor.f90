module mod_neighbor
  use global
  use omp_lib
  use ifport
  implicit none
contains
  
  subroutine neighbor_wrapper(step)
    implicit none
    integer :: step,signal_flag
    integer :: T1,T2,clock_rate,clock_max
    integer :: i
  
   
    nn_builds = nn_builds + 1
    drx(:) = 0.0d0
    dry(:) = 0.0d0
    drz(:) = 0.0d0


    call system_clock(T1,clock_rate,clock_max)
    call check_allocation

    call new_nlist(step)
    call sort_posit_charge(step)

    if(intel_flag.eq.1)then
       call check_cache_size
 

    else
       !call check_cache_size_backup_NOMIC

       call check_cache_size_nomic
    endif

    if(numbond.gt.0)then
       call sort_specbond_numbond(step)
       call sort_bondlist(step)
    endif


    
    if(intel_flag.eq.1)then
       signal_flag =1

    endif

    if(numangle.gt.0)then
       call sort_anglelist(step)
    endif

    if(numdihed.gt.0)then
       call sort_dihedlist(step)
    endif

    
    if(numimpro.gt.0)then
       call sort_improlist(step)
    endif


    if(nshake.gt.0)then
       call sort_shakeatom(step)
    endif


    call system_clock(T2,clock_rate,clock_max)
    time_sort = time_sort + real(T2-T1)/real(clock_rate)



    if(intel_flag.eq.1)then
       call build_neighbor_nonewton_MIC(step,0)

       !call build_neighbor_nonewton_MIC_ISPFULLPAR(step,0)
       !call build_neighbor_nonewton_MIC_backup(step,0)

    else
       call build_neighbor_nonewton_nomic(step)
       !call build_neighbor_nonewton_nomic_backup(step,0)
    endif

  
  
  end subroutine neighbor_wrapper



  subroutine build_neighbor_nonewton_MIC(step,signal_flag)
    implicit none
    real*4  :: x1,y1,z1,x2,y2,z2
    real*4  :: dx,dy,dz,dr2
    real*4  :: boxdx,boxdy,boxdz
    integer :: bonds(12)
    integer :: start_stencil,end_stencil
    integer :: check
    integer :: i,j,k,z,m
    integer :: c1s,c1e,c2s,c2e
    integer :: cell_neigh,step
    integer :: num2tmp,num3tmp
    integer :: neigh_flag
    integer :: nnpt,n0,neigh
    integer :: bin,count
    integer :: T1,T2,clock_rate,clock_max
    integer :: cached_bin,ncache,l
    integer :: off1,off2,off3,signal_flag
    integer :: redo,first,max, threadmax,flag_cycle
    integer :: flag1,flag2,data_flag
    
    redo = 1
    data_flag = 1
    do while(redo.ne.0)
 
       call system_clock(T1,clock_rate,clock_max)
   
       !dir$ offload begin target(mic:ourmic) inout(max),                       &
       !dir$ in(position: alloc_if(.false.),free_if(.false.)),                  &
       !dir$ nocopy(nlist: alloc_if(.false.) free_if(.false.)),                 &
       !dir$ nocopy(numneigh: alloc_if(.false.) free_if(.false.)),              &
       !dir$ in(atombin: alloc_if(.false.) free_if(.false.)),                   &
       !dir$ in(start: alloc_if(.false.) free_if(.false.)),                     &
       !dir$ in(endposit: alloc_if(.false.) free_if(.false.)),                  &
       !dir$ in(cnum: alloc_if(.false.) free_if(.false.)),                      &
       !dir$ in(num3bond: alloc_if(.false.) free_if(.false.)),                  &
       !dir$ in(specbond: alloc_if(.false.) free_if(.false.)),                  &
       !dir$ nocopy(dr2array: alloc_if(.false.) free_if(.false.))

       
       max =0
       ncache = 0
       first =0
       numneigh(:) = 0
       
       !$omp parallel default(firstprivate),               &
       !$omp& shared(position,atombin),                    &
       !$omp& shared(num3bond,specbond),                   &
       !$omp& shared(dr2array),                            &
       !$omp& shared(cnum,start,endposit),                 &
       !$omp& shared(nlist,numneigh),                      &
       !$omp& shared(first,max,redo,neigh_alloc)
       threadmax = 0

 
       
       !$omp do schedule(dynamic) 

       do i = 1,np
          x1 = position(i)%x
          y1 = position(i)%y
          z1 = position(i)%z
          
          bin      = atombin(i)
          num3tmp  = num3bond(i)
          bonds(:) = specbond(:,i)          
          n0   = (i-1)*neigh_alloc

  

          off1 = bin*cache_size
          off2 = 3*off1
          off3 = (i-1)*cache_size


          start_stencil = bin*27
          end_stencil   = 27*bin+26
          
          ncache = 0
          do k = start_stencil,end_stencil
             cell_neigh = cnum(k)            
             c2s = start(cell_neigh); c2e = endposit(cell_neigh)
             flag1 = off3 + ncache + 1 - c2s
       
             !dir$ simd
             do z= c2s,c2e                                  
                dx = x1-position(z)%x
                dy = y1-position(z)%y
                dz = z1-position(z)%z
                
                
                boxdx = dx*ibox; boxdy = dy*ibox; boxdz = dz*ibox
                boxdx = (boxdx+sign(1/(epsilon(boxdx)),boxdx)) -sign(1/epsilon(boxdx),dx)
                boxdy = (boxdy+sign(1/(epsilon(boxdy)),boxdy)) -sign(1/epsilon(boxdy),dy)
                boxdz = (boxdz+sign(1/(epsilon(boxdz)),boxdz)) -sign(1/epsilon(boxdz),dz)
                dx = dx-boxdx*box
                dy = dy-boxdy*box
                dz = dz-boxdz*box

                dr2 = dx*dx + dy*dy + dz*dz

                dr2array(z+flag1) = dr2
             enddo
             ncache = ncache + (c2e-c2s)+1
          enddo
          

          nnpt   = 0
          ncache = 0
          
          do k = start_stencil,end_stencil
             cell_neigh = cnum(k)            
             c2s = start(cell_neigh); c2e = endposit(cell_neigh)
             
             do z= c2s,c2e                                  
                
                ncache = ncache + 1
                dr2    = dr2array(off3+ncache)
                
                if(z.ne.i.and.dr2.lt.rv2)then
                   !---check to see if atoms are bonded
                   neigh_flag = 1
                   do l = 1,num3tmp
                      if(bonds(l).eq.z)then
                         neigh_flag = 0
                      endif
                   enddo
                   
                   if(neigh_flag.eq.1)then
                      nnpt = nnpt+1
                      
                      if(n0+nnpt.le.nlistsize)then
                         nlist(n0+nnpt) = z
                      endif
                      
                   endif
                endif
             enddo
          enddo      
    
    
          numneigh(i) = nnpt



          if(nnpt.gt.neigh_alloc)then        
             if(nnpt.gt.threadmax)then
                threadmax = nnpt           
             endif
          endif

         
         
       enddo
    
       !$omp end  do
       !$omp critical
       if(threadmax.gt.max)then
          max= threadmax
       endif
       !$omp end critical
       !$omp end parallel
       !dir$ end offload

       call system_clock(T2,clock_rate,clock_max)
       time_neigh = time_neigh + real(T2-T1)/real(clock_rate)


       if(max.gt.neigh_alloc)then
          write(1234,*)'we are in here',max,step
          neigh_alloc = max*1.2
          do while(mod(neigh_alloc,16).ne.0)
             neigh_alloc = neigh_alloc + 1
          enddo
          write(1234,*)'we are in with neigh_alloc',neigh_alloc
          print*,'DEALLOCATING NEIGHBOR LIST FOR RANK',rank
          print*,'NEW NEIGH ALLOC SIZE',neigh_alloc

          nlistsize = neigh_alloc*np
          !---deallocate on host

          if(intel_flag.eq.1)then
             !---dealocate on MIc
             !dir$ offload begin target(mic:ourmic),&
             !dir$ nocopy(nlist: alloc_if(.false.) free_if(.false.))
             
             deallocate(nlist)
             
             !dir$ end offload
          else
             deallocate(nlist)
          endif



          if(intel_flag.eq.1)then
             !---allocate on MIC
             !dir$ offload begin target(mic:ourmic) in(nlistsize) nocopy(nlist:length(0) alloc_if(.false.) free_if(.false.))
             allocate(nlist(nlistsize))
             !dir$ end offload
          else
             allocate(nlist(nlistsize))
          endif

          print*,'WE HAVE SUCCESFULLY REALLOCATED FOR RANK',rank
         
       else
          redo = 0
       end if

    enddo

    
  end subroutine build_neighbor_nonewton_MIC









  subroutine check_cache_size
    implicit none
    integer :: i,z,ncache,start_stencil,end_stencil
    integer :: off1,off2,cell_neigh,c2s,c2e,max,k,threadmax
    integer :: T1,T2,clock_rate,clock_max

    call system_clock(T1,clock_rate,clock_max)
    !dir$ offload begin target(mic:ourmic),                          &
    !dir$ inout(cache_size),                                         &
    !dir$ in(ncellT,np),                                             &
    !dir$ in(cnum: alloc_if(.false.) free_if(.false.)),              &
    !dir$ in(start: alloc_if(.false.) free_if(.false.)),             &
    !dir$ in(endposit: alloc_if(.false.) free_if(.false.)),          &
    !dir$ nocopy(dr2array: alloc_if(.false.) free_if(.false.))         
   
    max = 0
    threadmax = 0

    !$omp parallel default(firstprivate),&
    !$omp& shared(start,endposit,max)
    !$omp do schedule(dynamic)
    do i =0,ncellT-1
       ncache = 0
       start_stencil = i*27
       end_stencil   = 27*i+26
       
       
       do k = start_stencil,end_stencil
          cell_neigh = cnum(k)            
          c2s = start(cell_neigh); c2e = endposit(cell_neigh)
    
          do z= c2s,c2e                                  
             ncache = ncache + 1        
          enddo
       enddo

       if(ncache.gt.cache_size)then
          if(ncache.gt.threadmax)then
             threadmax = ncache
          endif
       endif
    enddo
    !$omp end  do

    !$omp critical
    if(threadmax.gt.max)then
       max = threadmax
    endif
    !$omp end critical
    !$omp end parallel


    if(max.ne.0)then
  
       !deallocate(dr2AOS)
       deallocate(dr2array)
       cache_size = max
       do while(mod(cache_size,16).ne.0)
          cache_size = cache_size+1
       enddo
       write(1234,*)'say hello to the new cache size on MIC',cache_size,mod(cache_size,16)
       write(1234,*)'what is neigh alloc',neigh_alloc
       write(1234,*)'what is nlistsize',nlistsize
       write(1234,*)'don reallocatinge with cache'
       allocate(dr2array(np*cache_size))
       write(1234,*)'we have allocated dr2array'

     

    endif
    !dir$ end offload
    call system_clock(T2,clock_rate,clock_max)
    
    time_cache = time_cache + real(T2-T1)/real(clock_rate)

  end subroutine check_cache_size




  
  subroutine check_cache_size_backup
    implicit none
    integer :: i,z,ncache,start_stencil,end_stencil
    integer :: off1,off2,cell_neigh,c2s,c2e,max,k,threadmax
    integer :: T1,T2,clock_rate,clock_max
    
    call system_clock(T1,clock_rate,clock_max)
    !dir$ offload begin target(mic:ourmic),&
    !dir$ inout(cache_size),&
    !dir$ in(ncellT,np),&
    !dir$ in(cnum: alloc_if(.false.) free_if(.false.)),&
    !dir$ in(start: alloc_if(.false.) free_if(.false.)),&
    !dir$ in(endposit: alloc_if(.false.) free_if(.false.)),&
    !dir$ nocopy(pptr_cache: alloc_if(.false.) free_if(.false.)),&
    !dir$ nocopy(cache: alloc_if(.false.) free_if(.false.)),&
    !dir$ nocopy(dr2array: alloc_if(.false.) free_if(.false.)),&
    !dir$ nocopy(nncache: alloc_if(.false.) free_if(.false.))
    
    max = 0
    threadmax = 0
    
    !$omp parallel default(firstprivate),&
    !$omp& shared(start,endposit,max)
    !$omp do schedule(dynamic)
    do i =0,ncellT-1
       ncache = 0
       start_stencil = i*27
       end_stencil   = 27*i+26
       
       
       do k = start_stencil,end_stencil
          cell_neigh = cnum(k)            
          c2s = start(cell_neigh); c2e = endposit(cell_neigh)
          !ncache = (c2e-c2s) + 1
    
          do z= c2s,c2e                                  
             ncache = ncache + 1        
          enddo
       enddo

       if(ncache.gt.cache_size)then
          if(ncache.gt.threadmax)then
             threadmax = ncache
          endif
       endif
    enddo
    !$omp end  do

    !$omp critical
    if(threadmax.gt.max)then
       max = threadmax
    endif
    !$omp end critical
    !$omp end parallel


    if(max.ne.0)then
  
       deallocate(pptr_cache,cache,dr2array,nncache)

       cache_size = max
       do while(mod(cache_size,16).ne.0)
          cache_size = cache_size+1
       enddo

       print*,'say hello to the new cache size on MIC',cache_size,mod(cache_size,16)

       allocate(pptr_cache(ncellT*cache_size))
       allocate(cache(3*ncellT*cache_size))
       allocate(dr2array(np*cache_size))
       allocate(nncache(0:ncellT*cache_size))
    endif
    !dir$ end offload
    call system_clock(T2,clock_rate,clock_max)

    time_cache = time_cache + real(T2-T1)/real(clock_rate)
  !  print*,'elapsed time in cache',real(T2-T1)/real(clock_rate)




  end subroutine check_cache_size_backup


    subroutine check_cache_size_backup_NOMIC
    implicit none
    integer :: i,z,ncache,start_stencil,end_stencil
    integer :: off1,off2,cell_neigh,c2s,c2e,max,k,threadmax
    integer :: T1,T2,clock_rate,clock_max
    
    call system_clock(T1,clock_rate,clock_max)
    
    max = 0
    threadmax = 0
    

    do i =0,ncellT-1
       ncache = 0
       start_stencil = i*27
       end_stencil   = 27*i+26
       
       
       do k = start_stencil,end_stencil
          cell_neigh = cnum(k)            
          c2s = start(cell_neigh); c2e = endposit(cell_neigh)
          !ncache = (c2e-c2s) + 1
    
          do z= c2s,c2e                                  
             ncache = ncache + 1        
          enddo
       enddo

       if(ncache.gt.cache_size)then
          if(ncache.gt.threadmax)then
             threadmax = ncache
          endif
       endif
    enddo

    if(threadmax.gt.max)then
       max = threadmax
    endif



    if(max.ne.0)then
  
       deallocate(pptr_cache,cache,dr2array,nncache)

       cache_size = max
       do while(mod(cache_size,16).ne.0)
          cache_size = cache_size+1
       enddo

       print*,'say hello to the new cache size on MIC',cache_size,mod(cache_size,16)

       allocate(pptr_cache(ncellT*cache_size))
       allocate(cache(3*ncellT*cache_size))
       allocate(dr2array(np*cache_size))
       allocate(nncache(0:ncellT*cache_size))
    endif
    call system_clock(T2,clock_rate,clock_max)

    time_cache = time_cache + real(T2-T1)/real(clock_rate)
  !  print*,'elapsed time in cache',real(T2-T1)/real(clock_rate)

    print*,'done with cache size'


  end subroutine check_cache_size_backup_NOMIC





  subroutine check_cache_size_nomic
    implicit none
    real*4  :: max_mass,cur_mass
    integer :: i,z,ncache,start_stencil,end_stencil
    integer :: off1,off2,cell_neigh,c2s,c2e,max,k,threadmax
    integer :: T1,T2,clock_rate,clock_max

    call system_clock(T1,clock_rate,clock_max)
 
    
    max = 0
    threadmax = 0
    max_mass = 0.0
    do i =0,ncellT-1
       ncache = 0
       start_stencil = i*27
       end_stencil   = 27*i+26
       
       cur_mass = 0.0

       do k = start_stencil,end_stencil
          cell_neigh = cnum(k)            
          c2s = start(cell_neigh); c2e = endposit(cell_neigh)
    
          do z= c2s,c2e                                  
             ncache = ncache + 1        

             cur_mass = cur_mass + mass(position(z)%type)

          enddo
       enddo

       if(cur_mass.gt.max_mass)then
          max_mass = cur_mass
       endif

       if(ncache.gt.cache_size)then
          if(ncache.gt.threadmax)then
             threadmax = ncache
          endif
       endif

    enddo

    if(threadmax.gt.max)then
       max = threadmax
    endif


    if(max.ne.0)then
  
       !deallocate(dr2AOS)
       deallocate(dr2array)

       cache_size = max
       do while(mod(cache_size,16).ne.0)
          cache_size = cache_size+1
       enddo

       
       allocate(dr2array(np*cache_size))

       !allocate(dr2AOS(np*cache_size))
       
    endif
    call system_clock(T2,clock_rate,clock_max)

    time_cache = time_cache + real(T2-T1)/real(clock_rate)

  end subroutine check_cache_size_nomic



  subroutine build_neighbor_nonewton_nomic(step)
    implicit none
    real*4 :: x1,y1,z1,x2,y2,z2
    real*4:: dx,dy,dz,dr2
    real*4 :: boxdx,boxdy,boxdz
    integer :: bonds(12)
    integer :: start_stencil,end_stencil
    integer :: check
    integer :: i,j,k,z,m
    integer :: c1s,c1e,c2s,c2e
    integer :: cell_neigh,step
    integer :: num1tmp,num2tmp,num3tmp
    integer :: neigh_flag
    integer :: nnpt,n0,neigh
    integer :: bin,count
    integer :: T1,T2,clock_rate,clock_max
    integer :: cached_bin,ncache,l
    integer :: off1,off2,off3,signal_flag
    integer :: redo,first,max, threadmax,flag_cycle
    integer :: flag1,flag2
    


    redo = 1
    do while(redo.ne.0)
 
       call system_clock(T1,clock_rate,clock_max)
   
     
       max =0
       ncache = 0
       first =0
       numneigh(:) = 0


    
       threadmax = 0

 
       do i = 1,np
          x1 = position(i)%x
          y1 = position(i)%y
          z1 = position(i)%z
          
          bin      = atombin(i)
          num3tmp  = num3bond(i)
          bonds(:) = specbond(:,i)          
          n0   = (i-1)*neigh_alloc

  

          off1 = bin*cache_size
          off2 = 3*off1
          off3 = (i-1)*cache_size


          start_stencil = bin*27
          end_stencil   = 27*bin+26
          
          ncache = 0
          do k = start_stencil,end_stencil
             cell_neigh = cnum(k)            
             c2s = start(cell_neigh); c2e = endposit(cell_neigh)
             flag1 = off3 + ncache + 1 - c2s
             
             do z= c2s,c2e                                  
                dx = x1-position(z)%x
                dy = y1-position(z)%y
                dz = z1-position(z)%z
                
                !ncache = ncache + 1
                
                boxdx = dx*ibox; boxdy = dy*ibox; boxdz = dz*ibox
                boxdx = (boxdx+sign(1/(epsilon(boxdx)),boxdx)) -sign(1/epsilon(boxdx),dx)
                boxdy = (boxdy+sign(1/(epsilon(boxdy)),boxdy)) -sign(1/epsilon(boxdy),dy)
                boxdz = (boxdz+sign(1/(epsilon(boxdz)),boxdz)) -sign(1/epsilon(boxdz),dz)
                dx = dx-boxdx*box
                dy = dy-boxdy*box
                dz = dz-boxdz*box

                dr2 = dx*dx + dy*dy + dz*dz
                !dr2array(off3+ncache) = dr2

                dr2array(z+flag1) = dr2
             enddo
             ncache = ncache + (c2e-c2s)+1
          enddo
          

          nnpt   = 0
          ncache = 0
          
          do k = start_stencil,end_stencil
             cell_neigh = cnum(k)            
             c2s = start(cell_neigh); c2e = endposit(cell_neigh)
             
             do z= c2s,c2e                                  
                
                ncache = ncache + 1
                dr2    = dr2array(off3+ncache)
                
                if(z.ne.i.and.dr2.lt.rv2)then
                   !---check to see if atoms are bonded
                   neigh_flag = 1
                   do l = 1,num3tmp
                      if(bonds(l).eq.z)then
                         neigh_flag = 0
                      endif
                   enddo
                   
                   if(neigh_flag.eq.1)then
                      nnpt = nnpt+1
                      
                      if(n0+nnpt.le.nlistsize)then
                         nlist(n0+nnpt) = z
                      endif
                      
                   endif
                endif
             enddo
          enddo      
    
    
          numneigh(i) = nnpt

          if(nnpt.gt.neigh_alloc)then        
             if(nnpt.gt.threadmax)then
                threadmax = nnpt           
             endif
          endif

         
         
       enddo
    
    
       if(threadmax.gt.max)then
          max= threadmax
       endif

       call system_clock(T2,clock_rate,clock_max)
       !print*,'elapsed time building neighbor test',real(T2-T1)/real(clock_rate)
       time_neigh = time_neigh + real(T2-T1)/real(clock_rate)


       if(max.gt.neigh_alloc)then
          print*,'we are in here',max,step
          neigh_alloc = max*1.2
          do while(mod(neigh_alloc,16).ne.0)
             neigh_alloc = neigh_alloc + 1
          enddo
          print*,'we are in with neigh_alloc',neigh_alloc

          nlistsize = neigh_alloc*np
          !---deallocate on host

          if(intel_flag.eq.1)then
             !---dealocate on MIc
             !dir$ offload begin target(mic:ourmic),&
             !dir$ nocopy(nlist: alloc_if(.false.) free_if(.false.))
             
             deallocate(nlist)
             
             !dir$ end offload
          else
             deallocate(nlist)
          endif



          if(intel_flag.eq.1)then
             !---allocate on MIC
             !dir$ offload begin target(mic:ourmic) in(nlistsize) nocopy(nlist:length(0) alloc_if(.false.) free_if(.false.))
             allocate(nlist(nlistsize))
             !dir$ end offload
          else
             allocate(nlist(nlistsize))
          endif
       else
          redo = 0
       end if

    enddo
       
     
     
  
  end subroutine build_neighbor_nonewton_nomic




  subroutine build_neighbor_n2_nonewton(step)
    implicit none
    integer :: i,j,nnpt,neigh_flag,m
    integer :: off,step,num3tmp
    real*4  :: dx,dy,dz,dr2
    real*4 :: boxdx,boxdy,boxdz
    real*4 :: x1,y1,z1,x2,y2,z2

    numneigh(:) =0

    do i = 1,np
       x1 = position(i)%x; y1 =position(i)%y; z1 = position(i)%z
       off = (i-1)*neigh_alloc
       num3tmp = num3bond(i)
       nnpt = 0
       do j = 1,np
          if(i.eq.j)cycle

          dx = x1-position(j)%x
          dy = y1-position(j)%y
          dz = z1-position(j)%z
          boxdx = dx*ibox; boxdy = dy*ibox; boxdz = dz*ibox
          boxdx = (boxdx+sign(1/(epsilon(boxdx)),boxdx)) -sign(1/epsilon(boxdx),dx)
          boxdy = (boxdy+sign(1/(epsilon(boxdy)),boxdy)) -sign(1/epsilon(boxdy),dy)
          boxdz = (boxdz+sign(1/(epsilon(boxdz)),boxdz)) -sign(1/epsilon(boxdz),dz)
          dx = dx-box*boxdx; dy = dy-box*boxdy; dz = dz-box*boxdz

          dr2 = dx*dx + dy*dy + dz*dz

          neigh_flag = 0
          do m = 1,num3tmp
             if(specbond(m,i).eq.j)then
                neigh_flag = 1
             endif
          enddo
         
          if(neigh_flag.eq.0)then
             if(dr2.lt.rv2)then
                nnpt = nnpt+1
                nlist(off+nnpt) = j
             endif
          endif

       enddo
       numneigh(i)= nnpt
    enddo

  end subroutine build_neighbor_n2_nonewton

  
  subroutine check_neighbor_safe(nflag,step)
    implicit none
    double precision :: dx,dy,dz,dr2
    double precision :: max1,max2
    double precision :: dis
    integer :: nflag
    integer :: i,step
    
    
    max1 = 0.0d0; max2 = 0.0d0
    nflag = 0
    do i =1,np
       dr2  = drx(i)*drx(i) + dry(i)*dry(i) + drz(i)*drz(i)
       
       if(dr2.gt.max1)then
          max1 = dr2
       elseif(dr2.gt.max2)then
          max2 = dr2
       endif
    enddo

    dis = (max1 + max2)
    if(dis.ge.rvrc2)then
       nflag = 1 
    endif
  end subroutine check_neighbor_safe





  subroutine check_neighbor(nflag,step)
    implicit none
    double precision :: dx,dy,dz,dr2
    double precision :: max1,max2
    double precision :: dis
    integer :: nflag
    integer :: i,step
    
    
    max1 = 0.0d0; max2 = 0.0d0
    nflag = 0
    do i =1,np
       dr2  = drx(i)*drx(i) + dry(i)*dry(i) + drz(i)*drz(i)
       
       if(dr2.ge.rvrc2)then
          nflag = 1 
       endif
    enddo
  end subroutine check_neighbor






  !........................................!
  !     constructs LL and HOC lists        !
  !........................................!
  subroutine new_nlist(step)
    
    implicit none 
    integer :: i,icell,step
    
    !    -- initialize head of chain
    HOC(:) = 0
    
    
    !    -- make linked list
    do i=1,np
       
       !    -- determine cell number
       icell = xyz2bin(position(i)%x,position(i)%y,position(i)%z)

       !    -- update linked-list and head of chain
       LL(i) = HOC(icell)
       HOC(icell) = i
    enddo
  end subroutine new_nlist

  
  
  subroutine sort_posit_charge(step)
    integer :: i,j,k,count,step
    integer :: npz,opz
    integer :: tag1,tag2,tag3
    integer :: flag1,flag2


    count = 1
    do i = 0,ncellT-1
       
       j = HOC(i)
       start(i) = count
       
       do while (j.ne.0)
          
          tt(j) = count
          !ttb(count) = j

          ps(count)%x = position(j)%x; ps(count)%y = position(j)%y; ps(count)%z = position(j)%z
          ps(count)%type = position(j)%type
 
          qs(count) = q(j)
          vs(count)%x = v(j)%x; vs(count)%y = v(j)%y; vs(count)%z = v(j)%z
          
          molsort(:,count) = mol(:,j)
          atombin(count) = i
          
          global_id_ss(count) = global_id(j)
          
          integrate_flag_ss(count) = integrate_flag(j)

          
          j = LL(j)
          count = count + 1
       enddo
       endposit(i) = count - 1
    end do
    

    !---sort position
    do i = 1,np
       position(i)%x = ps(i)%x
       position(i)%y = ps(i)%y
       position(i)%z = ps(i)%z
       position(i)%type = ps(i)%type


       v(i)%x = vs(i)%x; v(i)%y = vs(i)%y; v(i)%z = vs(i)%z
    enddo



    
    !---sort charge
    global_id(:) = global_id_ss(:)
    q(1:np) = qs(1:np)
    mol(:,:) = molsort(:,:)
    integrate_flag(:) = integrate_flag_ss(:)
   
  end subroutine sort_posit_charge


  


  subroutine sort_specbond_numbond(step)
    implicit none
    integer :: step
    integer :: i,j,k
    integer :: np2,np1

    !---transfer data
    do i = 1,np

       np1 = tt(i)

       !---sort number of bond arrays
       num1bondss(np1) = num1bond(i)
       num2bondss(np1) = num2bond(i)
       num3bondss(np1) = num3bond(i)
       
       do j = 1,num3bond(i)

          np2 = tt( specbond(j,i) )
          specbond_ss(j,np1)  = np2

       enddo
    enddo

    specbond(:,:) = specbond_ss(:,:)
    num1bond(:) = num1bondss(:)
    num2bond(:) = num2bondss(:)
    num3bond(:) = num3bondss(:)
  end subroutine sort_specbond_numbond



  subroutine sort_improlist(step)
    implicit none
    integer :: step
    integer :: i,j,k,count
    integer :: npz,opz
    integer :: tag1,tag2,tag3
    integer :: flag1,flag2
    
    !---move data to sort array
    do i = 1,numimpro
       
       do j = 1,4
          npz = tt( improlist(j,i) )
          
          improlist(j,i) = npz
          
       end do
    enddo
  end subroutine sort_improlist



  subroutine sort_dihedlist(step)
    implicit none
    integer :: step
    integer :: i,j,k,count
    integer :: npz,opz
    integer :: tag1,tag2,tag3
    integer :: flag1,flag2
    
    !---move data to sort array
    do i = 1,numdihed
       
       do j = 1,4

          npz = tt( dihedlist(j,i) )
          
          dihedlist(j,i) = npz
          
       end do
    enddo
  end subroutine sort_dihedlist
  
  
  subroutine sort_anglelist(step)
    implicit none
    integer :: step
    integer :: i,j,k,count
    integer :: npz,opz
    integer :: tag1,tag2,tag3
    integer :: flag1,flag2
    
    !---move data to sort array
    do i = 1,numangle
       
       do j = 1,3
          npz = tt( anglelist(j,i) )
          
          anglelist(j,i) = npz
          
       end do
    enddo
  end subroutine sort_anglelist

  subroutine sort_bondlist(step)
    implicit none
    integer :: step
    integer :: i,j,npz
    
    !---move data to sort array
    do i = 1,numbond
       
       do j = 1,2
          npz = tt( bondlist(j,i) )
          
          bondlist(j,i) = npz
          
       end do
    enddo

  end subroutine sort_bondlist


  subroutine sort_shakeatom(step)
    implicit none
    integer :: step
    integer :: num
    integer :: i,j,npz,npz2,npz3,npz4
    
    !---move data to sort array
    do i = 1,nshake
       

       num = shakeatom(i)%num
       
       if(num.eq.1)then
          npz  = tt( shakeatom(i)%atom1 )
          npz2 = tt( shakeatom(i)%atom2 )
          
          shakeatom(i)%atom1 = npz
          shakeatom(i)%atom2 = npz2
       elseif(num.eq.-3.or.num.eq.3)then
          npz  = tt( shakeatom(i)%atom1 )
          npz2 = tt( shakeatom(i)%atom2 )
          npz3 = tt( shakeatom(i)%atom3 )

          
          shakeatom(i)%atom1 = npz
          shakeatom(i)%atom2 = npz2
          shakeatom(i)%atom3 = npz3
       elseif(num.eq.4)then
          npz  = tt( shakeatom(i)%atom1 )
          npz2 = tt( shakeatom(i)%atom2 )
          npz3 = tt( shakeatom(i)%atom3 )
          npz4 = tt( shakeatom(i)%atom4 ) 

          
          shakeatom(i)%atom1 = npz
          shakeatom(i)%atom2 = npz2
          shakeatom(i)%atom3 = npz3
          shakeatom(i)%atom4 = npz4
       endif
          
    enddo

  end subroutine sort_shakeatom

  subroutine build_bond_list(step)
    implicit none
    integer :: i,j,k,l
    integer :: num,tag
    integer :: count
    integer :: type1,type2
    integer :: bond
    integer :: step



    count = 0
    do i = 1 ,np
       num = num1bond(i)
       type1 = position(i)%type

       do j = 1,num
          tag   = specbond(j,i)
          
          if(tag.gt.i)then
             
             count = count + 1
             type2 = position(tag)%type
             bond  = atom2bond(type1,type2) 

             
             bondlist(1,count) = i
             bondlist(2,count) = tag
             bondlist(3,count) = bond
          endif
       enddo
    enddo



  end subroutine build_bond_list

  


  !------------------------------------------------------------------------------------------------------------
  !*** routine calculates the cell pairs where i>j
  !*** 
  subroutine  neighcell_newton
    
    implicit none
    integer :: offsets(3,13)
    integer :: in,i,count
    integer :: ic,itel,icz,icy,icx
    integer :: ix,iy,iz
    integer :: iccx,iccy,iccz
    integer :: x,y,z
    


    !--- build offsets array
    offsets(1,1)=1; offsets(2,1)=0;  offsets(3,1) = 0
    offsets(1,2)=1; offsets(2,2)=1;  offsets(3,2) = 0
    offsets(1,3)=1; offsets(2,3)=-1; offsets(3,3) = 0
    offsets(1,4)=0; offsets(2,4)=-1; offsets(3,4) = 0
    
    offsets(1,5)=1; offsets(2,5)=0;  offsets(3,5) = 1
    offsets(1,6)=1; offsets(2,6)=1;  offsets(3,6) = 1
    offsets(1,7)=1; offsets(2,7)=-1; offsets(3,7) = 1
    offsets(1,8)=0; offsets(2,8)=-1; offsets(3,8) = 1
    
    offsets(1,9)=1;  offsets(2,9)=0;  offsets(3,9) = -1
    offsets(1,10)=1; offsets(2,10)=1; offsets(3,10) = -1
    offsets(1,11)=1; offsets(2,11)=-1;offsets(3,11) = -1
    offsets(1,12)=0; offsets(2,12)=-1;offsets(3,12) = -1
    
    offsets(1,13) = 0; offsets(2,13)=0; offsets(3,13) = 1
    
    
    count = 0
    do ic = 0,ncellT-1
             
       
       !--- obtain x,y,z coordinates of cell of interest: ic
       !--- this is doing using integer arithmetic 
       icz = ic/(ncellD*ncellD)
       icy = (ic - icz*ncellD*ncellD)/(ncellD)
       icx = ic -icy*ncellD - icz*ncellD*ncellD
       
       
       
       do i = 1,13
          
          
          ix = offsets(1,i); iy = offsets(2,i); iz = offsets(3,i)
          
          iccz = icz+iz
          
          if(iccz.lt.0) then
             z = iccz + ncellD                              
          else if (iccz.ge.ncellD) then
             z = iccz - ncellD                              
          else
             z = iccz
          endif
          
          
          iccy = icy+iy
          
          if(iccy.lt.0) then
             y = iccy + ncellD                                            
          else if (iccy.ge.ncellD) then
             y = iccy - ncellD                                                                                 
          else
             y = iccy
          endif
          
          
          iccx = icx+ix
          if(iccx.lt.0) then
             x = iccx + ncellD             
          else if (iccx.ge.ncellD) then
             x = iccx - ncellD                      
          else
             x = iccx    
          endif
          
          in = (x+y*ncellD + z*ncellD*ncellD)   
          
          cnum(count) = in

          count = count + 1
       enddo
     
    enddo

    
  end subroutine neighcell_newton
 
  

  !------------------------------------------------------------------------------------------------------------
  !*** routine calculates the all cell pairs
  !*** 
  subroutine  neighcell_nonewton
    
    implicit none
    integer :: in,i,count
    integer :: ic,itel,icz,icy,icx
    integer :: ix,iy,iz
    integer :: iccx,iccy,iccz
    integer :: x,y,z
    integer :: xc,yc,zc
    
    count = 0
    do ic = 0,ncellT-1
             
       cnum(count) = ic

       !---obtain x y z coordinates of cell of interest
       !---this is done using integer arithmethic
       icz =  ic/(ncellD**2)
       icy = (ic - icz*ncellD*ncellD)/(ncellD)
       icx =  ic - icy*ncellD-icz*ncellD*ncellD
       
       itel = 0
       
       do iz = -1,1
          do iy = -1,1
             do ix = -1,1
                
                iccz = icz+iz
                if(iccz.lt.0)then
                   zc = iccz+ncellD
                elseif(iccz.ge.ncellD)then
                   zc = iccz-ncellD
                else
                   zc = iccz
                endif
                
                
                iccy = icy+iy             
                if(iccy.lt.0)then
                   yc = iccy+ncellD
                elseif(iccy.ge.ncellD)then
                   yc = iccy-ncellD
                else
                   yc = iccy
                endif
                
                iccx = icx+ix             
                if(iccx.lt.0)then
                   xc = iccx+ncellD
                elseif(iccx.ge.ncellD)then
                   xc = iccx-ncellD
                else
                   xc = iccx
                endif

                if(iz.ne.0.or.iy.ne.0.or.ix.ne.0)then

                   count = count + 1
                   in = (xc+yc*ncellD+zc*ncellD*ncellD)      
                   cnum(count) = in

                endif

             enddo
          enddo
       end do
       count = count + 1
    end do
    
  end subroutine neighcell_nonewton

  
  integer function xyz2bin(x,y,z)
    implicit none
    real*4  :: x,y,z
    integer :: ix,iy,iz,numlength
    double precision :: xdp,ydp,zdp,irn,boxdp

    !===lets convert to higher precision for binning of atoms
    xdp =x; ydp =y ; zdp = z; irn = 1.0d0/rn; boxdp = box

    if(xdp.lt.boxdp)then    
       ix = int(xdp*irn)
    else
       ix = ncellD-1
    endif
    
    if(ydp.lt.boxdp)then    
       iy = int(ydp*irn)
    else
       iy = ncellD-1
    endif
    
    if(zdp.lt.boxdp)then    
       iz = int(zdp*irn)
    else
       iz = ncellD-1
    endif
    
    if(ix.eq.ncellD)then
       print*,'initializing x xyz2bin error'
       ix = ncellD-1
    endif
    if(iy.eq.ncellD)then
       print*,'initializing y xyz2bin error'
       iy = ncellD-1
    endif
    if(iz.eq.ncellD)then
       print*,'initializing z xyz2bin error'
       iz = ncellD-1
    endif

    xyz2bin = ix+iy*ncellD+iz*ncellD*ncellD
    
    return
  end function xyz2bin
  
  

  subroutine dealloc_list
    implicit none
    deallocate(ll,HOC,start,endposit,cnum,vlist,numneigh)
  end subroutine dealloc_list

    
  subroutine check_allocation
    implicit none
    real*4 :: current_box
    integer :: ncell_cur
    logical :: exist
    


    ncell_cur = ncellT

    
    !--- check number of cells
    if(cut_coul.gt.rcut)then
       rn         = box/int(box/cut_coul)
    else
       rn = box/int(box/rcut)
    endif
    ncellT      = nint(vol/(rn**3))
    ncellD      = nint(ncellT**(1.0d0/3.0d0))



    if(ncellT.gt.ncell_alloc)then

      
       ncell_alloc = (ncellD +1)*(ncellD+1)*(ncellD+1)




       !=======================================================
       ! BEGIN FILE WRITE
       inquire(file=trim(data_files_path)//'/'//'neigh_realloc.dat',exist=exist)
       if(exist)then
          open(unit = 99999, file = trim(data_files_path)//'/'//'neigh_realloc.dat',status='old',position='append',action='write')
       else
          open(unit = 99999, file = trim(data_files_path)//'/'//'neigh_realloc.dat',status='new',action='write')
       endif
       write(99999,*)'we are reallocating neighbor list',rank
       write(99999,*)'old number of cells',ncell_cur
       write(99999,*)'new number of cells',ncellT
       write(99999,*)'new ncell alloc',ncell_alloc



       !========================================================================
       ! DEALLOCATED REALLOCATE
       if(intel_flag.eq.1)then
          !dir$ offload_transfer target(mic:ourmic)                   &
          !dir$ in(start: alloc_if(.false.) free_if(.true.)),         &
          !dir$ in(endposit: alloc_if(.false.) free_if(.true.)),      &
          !dir$ in(cnum: alloc_if(.false.) free_if(.true.))                                            
          
          deallocate(HOC,start,endposit,cnum)
          allocate(HOC(0:ncell_alloc),start(0:ncell_alloc),endposit(0:ncell_alloc),&
               cnum(0:27*ncell_alloc))
          call neighcell_nonewton


          !dir$ offload_transfer target(mic:ourmic)                   &
          !dir$ in(start: alloc_if(.true.) free_if(.false.)),         &
          !dir$ in(endposit: alloc_if(.true.) free_if(.false.)),      &
          !dir$ in(cnum: alloc_if(.true.) free_if(.false.)) 


       else
          
          deallocate(HOC,start,endposit,cnum)
          allocate(HOC(0:ncell_alloc),start(0:ncell_alloc),endposit(0:ncell_alloc),&
               cnum(0:27*ncell_alloc))
          call neighcell_nonewton
       endif

       !===========================================================================
       ! TRACK SUCCESFUL MEMORY REALLOCATION
       write(99999,*)'WE HAVE SUCCESFULLY REALLOCATED THE MEMORY'
       write(99999,*)
       write(99999,*)
       write(99999,*)
       close(99999)




    elseif(ncellT.ne.ncell_cur)then
       

       call neighcell_nonewton

    endif


  end subroutine check_allocation

  subroutine init_list
    implicit none
    double precision :: scale
    integer :: i

    rn         = box/int(box/cut_coul)
 !   rn = box/int(box/rcut)
    ncellT     = nint(vol/(rn**3))
    ncellD     = nint(ncellT**(1.0d0/3.0d0))
    scale      = 2.0d0
    ncell_alloc=ncellT*1.2

 

    call print_init_list
    
    
    
  end subroutine init_list


  subroutine print_init_list
    implicit none
    write(1234,*)
    write(1234,*)'--------------- INTIALIZING LISTS-----------------------------------'
    write(1234,*)'           box               :',box
    write(1234,*)'           rn                :',rn
    write(1234,*)'           rv                :',rv
    write(1234,*)'           ncellT            :',ncellT
    write(1234,*) '          ncellD            :',ncellD
    write(1234,*) '          scale            :', scale
    write(1234,*)'           nghbr list buff  :', neigh_alloc
    write(1234,*)'           MAKE SURE NCELLD CUBED IS NCELLT'
    write(1234,*)'--------------- END INTIALIZING LISTS-------------------------------'
    write(1234,*)
  end subroutine print_init_list
  
  
end module mod_neighbor
