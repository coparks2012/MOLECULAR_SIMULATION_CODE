module mod_neighbor_cluster
  
  use global
  use mod_pbc
contains
  
  
  !........................................!
  !     constructs LL and HOC lists        !
  !     step is step number in code
  !     flag is mol type
  !........................................!

  
  subroutine cluster_neighbor_wrapper(step,flag)
    implicit none
    integer :: step, flag,i
    integer :: part


    call check_allocation_cluster
    call mol_new_nlist(step,flag)
    call mol_build_neighbor_nosort(step)


    !open(unit = 99999, file = 'dens.dat')
    !do i= 1 ,nummolarray(1)
    !   write(99999,*)denstrack(i),com(i)%x,com(i)%y,com(i)%z
    !enddo
    !close(99999)
   ! call build_neighbor_n2_triclinic
   ! call histogram_cluster_neighbors


    !call force_n2_dens(0)
    !call force_vlist_dens
  end subroutine cluster_neighbor_wrapper

  
  subroutine cluster_neighbor_wrapper_debug(step,flag)
    implicit none
    integer :: step, flag
    integer :: part


    call check_allocation_cluster
    call mol_new_nlist(step,flag)
    call mol_build_neighbor_nosort(step)

    call force_n2_dens(0)                                                                                                                                         
    call force_vlist_dens    

  end subroutine cluster_neighbor_wrapper_debug



  subroutine mol_build_neighbor_nosort(step)
    implicit none
    double precision :: x1,y1,z1,x2,y2,z2
    double precision :: dx,dy,dz,dr2
    double precision :: dr2i,dr6i,dr12i
    double precision :: shiftx,shifty,shiftz
    integer :: i,j,k,l
    integer :: particle1,particle2
    integer :: cell2
    integer :: step
    integer :: cell_neigh
    integer :: num,densnum,offset


  
    mnum(:) = 0
    mvl(:) = 0
    denstrack(:) = 0
    densvl(:) = 0
    do i = 0,mncellT-1
      
       particle1 = mHOC(i)

       do while(particle1.ne.0)
          

          x1 = com(particle1)%x; y1 = com(particle1)%y; z1 = com(particle1)%z

          !--- loop over cell neighbors
          do k = i*26,26*i+25                                                                                                                 
             
             cell_neigh = mcnum(k)      
                     
             particle2 = mHOC(cell_neigh)
     
             
             do while(particle2.ne.0)
                x2 = com(particle2)%x; y2 = com(particle2)%y; z2 = com(particle2)%z
                
                dx = x1-x2
                dy = y1-y2
                dz = z1-z2
        
                !dx = x2-x1
                !dy = y2-y1
                !dz = z2-z1

                shiftx = -box*nint(dx*ibox)
                shifty = -box*nint(dy*ibox)
                shiftz = -box*nint(dz*ibox)
                
                dx = dx + shiftx
                dy = dy + shifty
                dz = dz + shiftz
                
                dr2 = dx*dx+dy*dy+dz*dz 

       
                if(dr2.lt.mrcut2)then  

                   num  = mnum(particle1) + 1
                   mnum(particle1) = num
                   offset  = mol_alloc*(particle1-1) + num
                   mvl(offset) = particle2
                   shift(offset)%x = -shiftx; shift(offset)%y = -shifty; shift(offset)%z = -shiftz

                   if(dr2.lt.denscut2)then
                      densnum = denstrack(particle1) + 1
                      denstrack(particle1) = densnum
                      offset  = densvlalloc*(particle1-1) + densnum
                      densvl(offset) = particle2
                   endif

                                         
                endif
                particle2 = mll(particle2)
             enddo
          enddo

          !---loop over same cell
          particle2 = mHOC(i)
          do while(particle2.ne.0)
             
             if(particle2.ne.particle1)then
                x2 = com(particle2)%x; y2 = com(particle2)%y; z2 = com(particle2)%z
                
                dx = x1-x2
                dy = y1-y2
                dz = z1-z2
                
                
                dr2 = dx*dx+dy*dy+dz*dz 
           
                if(dr2.lt.mrcut2)then  
                   
                   num  = mnum(particle1) + 1
                   mnum(particle1) = num
                   offset  = mol_alloc*(particle1-1) + num
                   mvl(offset) = particle2
                   shift(offset)%x = 0.0d0; shift(offset)%y = 0.0d0; shift(offset)%z = 0.0d0
                   
                   
                   if(dr2.lt.denscut2)then
                      densnum = denstrack(particle1) + 1
                      denstrack(particle1) = densnum
                      offset  = densvlalloc*(particle1-1) + densnum
                      densvl(offset) = particle2
                   endif
                endif
             endif

             particle2 = mll(particle2)
          enddo

          particle1 = mll(particle1)
       enddo
    end do


  end subroutine mol_build_neighbor_nosort

  subroutine mol_new_nlist(step,flag)
    
    implicit none 
    integer :: i,icell,step
    integer :: flag
    
    !    -- initialize head of chain
    mHOC(:) = 0
    
    
    !    -- make linked list
    do i=1,numMolArray(flag)
       
       !    -- determine cell number
       icell = mol_xyz2bin(com(i)%x,com(i)%y,com(i)%z)
       
       !    -- update linked-list and head of chain
       mll(i) = mHOC(icell)
       mHOC(icell) = i
    enddo
    
  end subroutine mol_new_nlist
  



  subroutine build_neighbor_n2_triclinic
    implicit none
    double precision :: dx,dy,dz,dr2
    double precision :: x1,y1,z1
    double precision :: x2,y2,z2
    double precision :: shiftx,shifty,shiftz
    integer :: i,j,offset,num,densnum
    
    mnum(:) = 0
    mvl(:) = 0
    denstrack(:) = 0
    densvl(:) = 0
    
    do i = 1,nummolarray(1)
       x1 = com(i)%x; y1 = com(i)%y; z1 = com(i)%z
     
       do j = 1,nummolarray(1)

          if(j.eq.i)cycle
          x2 = com(j)%x; y2 = com(j)%y; z2 = com(j)%z
          
          dx = x1-x2
          dy = y1-y2
          dz = z1-z2
          
          
          call minimum_image_triclinic_cluster(dx,dy,dz,shiftx,shifty,shiftz)

          dx = dx + shiftx
          dy = dy + shifty
          dz = dz + shiftz
          
          !call minimum_image_triclinic(dx,dy,dz)
          dr2 = dx*dx+dy*dy+dz*dz 
          
          
       !   print*,'what is distance'
       !   print*,x1,y1,z1
       !   print*,x2,y2,x2
       !   print*,dr2,shiftx,shifty,shiftz
       !   print*,sqrt(dr2),dr2,mrcut2
         
          if(dr2.lt.mrcut2)then  
             
             num  = mnum(i) + 1
             mnum(i) = num
             offset  = mol_alloc*(i-1) + num
             mvl(offset) = j
             shift(offset)%x = -shiftx; shift(offset)%y = -shifty; shift(offset)%z = -shiftz
             
             if(dr2.lt.denscut2)then
                densnum = denstrack(i) + 1
                denstrack(i) = densnum
                offset  = densvlalloc*(i-1) + densnum
                densvl(offset) = j
             endif
             
             
          endif
       enddo
       
    enddo
    
 
    
  end subroutine build_neighbor_n2_triclinic
  
  subroutine sort_mol(step,mol)
    integer :: i,j,k,count,step
    integer :: mol
    integer :: npz,opz
    integer :: tag1,tag2,tag3
    integer :: flag1,flag2
    
    
    count = 1
    do i = 0,mncellT-1
       
       j = mHOC(i)
       mstart(i) = count
       
       do while (j.ne.0)
          
          mttb(count) = j
          
          scom(count)%x = com(j)%x; scom(count)%y = com(j)%y; scom(count)%z = com(j)%z
          
          
          j = mll(j)
          count = count + 1
         enddo
         mendposit(i) = count - 1
      end do
      
      !---sort position
      do i = 1,numMolArray(mol)
         com(i)%x = scom(i)%x
         com(i)%y = scom(i)%y
         com(i)%z = scom(i)%z
      enddo
      
      
      
    end subroutine sort_mol
    
    
    subroutine mol_build_neighbor_newton(step)
      implicit none
      double precision :: x1,y1,z1,x2,y2,z2
      double precision :: dx,dy,dz,dr2
      double precision :: dr2i,dr6i,dr12i
      double precision :: shiftx,shifty,shiftz
      integer :: i,j,k,l,m
      integer :: c1s,c1e,c2s,c2e
      integer :: cell_neigh,step
      integer :: neigh_flag
      integer :: num
      integer :: offset
      
      mnum(:) = 0
      mvl(:) = 0

      !--- loop over all cells
      do i = 0,mncellT-1
         
         c1s   = mstart(i)
         c1e   = mendposit(i)
         
         !--- loop over cell neighbors
         do k = i*13,13*i+12                                                                                                                 
            
            cell_neigh = mcnum(k)       
            c2s = mstart(cell_neigh); c2e = mendposit(cell_neigh)       
            
            
            do j = c1s,c1e   
               x1  = com(j)%x; y1 = com(j)%y ; z1 = com(j)%z
           
               !dir$ simd
               do l= c2s,c2e
                  x2 = com(l)%x; y2 = com(l)%y; z2 = com(l)%z
                  
                  dx = x2-x1
                  dy = y2-y1
                  dz = z2-z1

                  shiftx = -box*nint(dx*ibox)
                  shifty = -box*nint(dy*ibox)
                  shiftz = -box*nint(dz*ibox)
                  
                  dx = dx + shiftx
                  dy = dy + shifty
                  dz = dz + shiftz
                  
                  dr2 = dx*dx+dy*dy+dz*dz 
                  
                  if(dr2.lt.mrcut2)then  
                     
                     num  = mnum(j) + 1
                     mnum(j) = num
                     offset  = mol_alloc*(j-1) + num
                     mvl(offset) = l
                     shift(offset)%x = shiftx; shift(offset)%y = shifty; shift(offset)%z = shiftz

                     num     = mnum(l) + 1
                     mnum(l) = num
                     offset  = mol_alloc*(l-1) + num
                     mvl(offset) = j

                     shift(offset)%x = -shiftx; shift(offset)%y = -shifty; shift(offset)%z = -shiftz

                     
                  endif
               enddo
            enddo
         enddo
         
         do j = c1s,c1e-1
            x1  = com(j)%x; y1 = com(j)%y ; z1 = com(j)%z
    
            !dir$ simd
            do l = j+1,c1e
               x2 = com(l)%x; y2 = com(l)%y; z2 = com(l)%z
               
               dx = x2-x1
               dy = y2-y1
               dz = z2-z1
               
               dr2 = dx*dx+dy*dy+dz*dz
               
               if(dr2.lt.mrcut2)then    
                  
                  num  = mnum(j) + 1
                  mnum(j) = num
                  offset  = mol_alloc*(j-1) + num
                  mvl(offset) = l
                  shift(offset)%x = 0.0d0; shift(offset)%y = 0.0d0; shift(offset)%z = 0.0d0


                  
                  num     = mnum(l) + 1
                  mnum(l) = num
                  offset  = mol_alloc*(l-1) + num
                  mvl(offset) = j
                  shift(offset)%x = 0.0d0; shift(offset)%y = 0.0d0; shift(offset)%z = 0.0d0
                                                      
               endif
               
            enddo
            
            if(mnum(j).gt.maxneigh)then
               maxneigh = mnum(j)
            endif
            
         end do
      enddo
      
      
      
      
      if(maxneigh.gt.0.80*mol_alloc)then
         print*,'we are getting close for cluster',maxneigh,mol_alloc
      endif
      
    end subroutine mol_build_neighbor_newton

    subroutine mol_build_neighbor_nonewton(step)
      implicit none
      double precision :: x1,y1,z1,x2,y2,z2
      double precision :: dx,dy,dz,dr2
      double precision :: dr2i,dr6i,dr12i
      double precision :: shiftx,shifty,shiftz
      integer :: i,j,k,l,m
      integer :: c1s,c1e,c2s,c2e
      integer :: cell_neigh,step
      integer :: neigh_flag
      integer :: num
      integer :: offset

      
      mnum(:) = 0
      mvl(:) = 0

      !--- loop over all cells
      do i = 0,mncellT-1
         
         c1s   = mstart(i)
         c1e   = mendposit(i)
         
         !--- loop over cell neighbors
         do k = i*26,26*i+25                                                                                                                 
            
            cell_neigh = mcnum(k)       

            c2s = mstart(cell_neigh); c2e = mendposit(cell_neigh)       
            
            
            do j = c1s,c1e   
               x1  = com(j)%x; y1 = com(j)%y ; z1 = com(j)%z
           
               !dir$ simd
               do l= c2s,c2e
                  x2 = com(l)%x; y2 = com(l)%y; z2 = com(l)%z
                  
                  dx = x2-x1
                  dy = y2-y1
                  dz = z2-z1

                  shiftx = -box*nint(dx*ibox)
                  shifty = -box*nint(dy*ibox)
                  shiftz = -box*nint(dz*ibox)
                  

                  dx = dx + shiftx
                  dy = dy + shifty
                  dz = dz + shiftz
                  
                  dr2 = dx*dx+dy*dy+dz*dz 
                  
                  if(dr2.lt.mrcut2)then  
                     
                     num  = mnum(j) + 1
                     mnum(j) = num
                     offset  = mol_alloc*(j-1) + num
                     mvl(offset) = l
                     shift(offset)%x = shiftx; shift(offset)%y = shifty; shift(offset)%z = shiftz


                  endif
               enddo
            enddo
         enddo
         
         do j = c1s,c1e-1
            x1  = com(j)%x; y1 = com(j)%y ; z1 = com(j)%z
    
            !dir$ simd
            do l = j+1,c1e
               x2 = com(l)%x; y2 = com(l)%y; z2 = com(l)%z
               
               dx = x2-x1
               dy = y2-y1
               dz = z2-z1
               
               dr2 = dx*dx+dy*dy+dz*dz
               
               if(dr2.lt.mrcut2)then    
                  
                  num  = mnum(j) + 1
                  mnum(j) = num
                  offset  = mol_alloc*(j-1) + num
                  mvl(offset) = l
                  shift(offset)%x = 0.0d0; shift(offset)%y = 0.0d0; shift(offset)%z = 0.0d0

                  
               endif
               
            enddo
            
            if(mnum(j).gt.maxneigh)then
               maxneigh = mnum(j)
               print*,maxneigh,j,mrcut,mrcut2
            endif
            
         end do
      enddo
      
      
      if(maxneigh.gt.0.80*mol_alloc)then
         print*,'we are getting close for cluster',maxneigh,mol_alloc
      endif
      
    end subroutine mol_build_neighbor_nonewton

    subroutine check_allocation_cluster
      implicit none
      real*4 :: current_box
      integer :: ncell_cur
      
      ncell_cur = mncellT
      
      !--- check number of cells
      mrn         = box/int(box/mrcut)
      mncellT     = nint(vol/(mrn**3))
      mncellD     = nint(mncellT**(1.0d0/3.0d0))
      
     
      if(mncellT.gt.mncell_alloc)then
         print*,'reallocating cell lists for clustering'
         mncell_alloc = (mncellD +1)*(mncellD+1)*(mncellD+1)
         print*,'old number of cells',ncell_cur
         print*,'new number of cells',mncellT
         print*,'new ncell alloc',mncell_alloc
         
            
         deallocate(mHOC,mstart,mendposit,mcnum)
         allocate(mHOC(0:mncell_alloc),mstart(0:mncell_alloc),mendposit(0:mncell_alloc),&
                 mcnum(0:27*mncell_alloc))
         call mol_neigh
            
            
     
         
      elseif(mncellT.ne.ncell_cur)then
         print*,'redoing neighbors for cell clustering'
         call mol_neigh
         
      endif
      
      
    end subroutine check_allocation_cluster
    

    !--- use this for when all pairs (i,j) need to be searched
    subroutine mol_neigh
      
      implicit none
      integer :: in
      integer :: icz,icy,icx,itel
      integer :: iccx,iccy,iccz
      integer :: ix,iy,iz
      integer :: xc,yc,zc,count,cell
      
      !--- initialize counter
      count = 0
      itel = 0


      do cell = 0,mncellT-1
         
         
         !---obtain x y z coordinates of cell of interest
         !---this is done using integer arithmethic
         icz =  cell/(mncellD**2)
         icy = (cell - icz*mncellD*mncellD)/(mncellD)
         icx =  cell - icy*mncellD-icz*mncellD*mncellD
         
         
         do iz = -1,1
            do iy = -1,1
               do ix = -1,1
                  
                  iccz = icz+iz
                  if(iccz.lt.0)then
                     zc = iccz+mncellD
                  elseif(iccz.ge.mncellD)then
                     zc = iccz-mncellD
                  else
                     zc = iccz
                  endif
                  
                  
                  iccy = icy+iy             
                  if(iccy.lt.0)then
                     yc = iccy+mncellD
                  elseif(iccy.ge.mncellD)then
                     yc = iccy-mncellD
                  else
                     yc = iccy
                  endif
                  
                  iccx = icx+ix             
                  if(iccx.lt.0)then
                     xc = iccx+mncellD
                  elseif(iccx.ge.mncellD)then
                     xc = iccx-mncellD
                  else
                     xc = iccx
                  endif
                  
                  !--- determine cell number with offset coordinates (ix,iy,iz)
                  !--- do not include cel self interactions here
                  
                  if(iz.ne.0.or.iy.ne.0.or.ix.ne.0)then
                     in = (xc+yc*mncellD+zc*mncellD*mncellD)
                     
                     mcnum(itel) = in
                     itel = itel +1

                  end if
    
               end do
            end do
         end do
      end do

    end subroutine mol_neigh
 


    subroutine histogram_cluster_neighbors
      implicit none
      integer :: i,j,k
      
      num_dens_ob = num_dens_ob+1
      do i = 1,nummolarray(1)
         histogram_dens( denstrack(i) ) = histogram_dens( denstrack(i) ) + 1
      enddo

  

    end subroutine histogram_cluster_neighbors
    
    
    integer function mol_xyz2bin(x,y,z)
      implicit none
      double precision :: imrn,boxdp
      double precision :: x,y,z
      integer :: ix,iy,iz,numlength
      
      boxdp = box; imrn = 1.0d0/mrn
      if(x.lt.boxdp)then    
         ix = int(x*imrn)
      else
         ix = mncellD-1
      endif
      
      if(y.lt.boxdp)then    
         iy = int(y*imrn)
      else
         iy = mncellD-1
      endif
      
      if(z.lt.boxdp)then    
         iz = int(z*imrn)
      else
         iz = mncellD-1
      endif

      if(ix.eq.mncellD)then
         ix = mncellD-1
         print*,'error in cluster with  x bin'
      endif
      if(iy.eq.mncellD)then
         iy = mncellD-1
         print*,'error in cluster with  y bin'
      endif
      if(iz.eq.mncellD)then
         iz = mncellD-1
         print*,'error in cluster with  z bin'
      endif
      
      mol_xyz2bin = ix+iy*mncellD+iz*mncellD*mncellD
      
      return
    end function mol_xyz2bin
    


    subroutine force_n2_COM(step)
    implicit none
    double precision :: force,forcelj
    double precision :: fudge_factor
    double precision :: x1,y1,z1,x2,y2,z2
    double precision :: ffx,ffy,ffz
    double precision :: dx,dy,dz,dr,dr2,dr2i,dr6i,dr12i
    double precision :: potential_N2
    integer :: i,j,step,l
    integer :: itype,jtype,icell
    integer :: offset
    integer :: num1tmp,num2tmp,num3tmp
    character(len=3) :: char,char2


    potential_N2 = 0.0d0
    do i = 1 ,numMolArray(1)-1
       x1 = com(i)%x; y1 = com(i)%y; z1 = com(i)%z;
       do j=i+1,numMolArray(1)
          
          x2 = com(j)%x; y2 = com(j)%y; z2 = com(j)%z; 
          
          dx = x1-x2; dy = y1-y2; dz = z1-z2
          dx = dx-box*nint(dx*ibox)
          dy = dy-box*nint(dy*ibox)
          dz = dz-box*nint(dz*ibox)

          dr2 = dx*dx + dy*dy + dz*dz

     
          if(dr2.lt.mrcut2)then
 
             dr2i = 1.0d0/dr2
             dr6i = dr2i*dr2i*dr2i
           
             potential_N2 = potential_N2 + dr6i*(dr6i-1.0d0)

          endif
       enddo
    enddo
    print*
    print*,'WHAT IS FINAL POTENTIAL N2',potential_N2,box,hbox,ibox
    print*
  end subroutine force_n2_COM



  subroutine force_vlist_com
    implicit none
    double precision :: force,forcelj
    double precision :: x1,y1,z1,x2,y2,z2
    double precision :: dx,dy,dz,dr,dr2,dr2i,dr6i,dr12i
    double precision :: ffx,ffy,ffz

    double precision :: fudge_factor
    double precision :: potential_v2
    integer :: i,j,l,step
    integer :: itype,jtype,neigh
    integer :: tid,num
    integer :: offset,icell
    integer :: num1tmp,num2tmp,num3tmp
    character(len=3) :: char,char2

    potential_v2 = 0.0d0
    do i = 1 ,numMolArray(1)
       x1 = com(i)%x; y1 = com(i)%y; z1 = com(i)%z

       print*,mnum(i)
       do j=1,mnum(i)
         

          neigh = mvl(j+mol_alloc*(i-1))
      
          if(neigh.gt.i)then
             
             x2 = com(neigh)%x; y2 = com(neigh)%y; z2 = com(neigh)%z
             
             dx = x1-x2; dy = y1-y2; dz = z1-z2
             dx = dx-box*nint(dx*ibox)
             dy = dy-box*nint(dy*ibox)
             dz = dz-box*nint(dz*ibox)
             
             dr2 = dx*dx + dy*dy + dz*dz
             
             
             if(dr2.lt.mrcut2)then
                
       
                dr2i = 1.0d0/dr2
                dr6i = dr2i*dr2i*dr2i
                potential_v2 = potential_v2 + dr6i*(dr6i-1.0d0)
         

             endif
          endif
       enddo
    enddo
    print*
    print*,'ENERGY VLIST',potential_v2,box,hbox,ibox
    print*
  end subroutine force_vlist_com
  
  


  subroutine force_n2_dens(step)
    implicit none
    double precision :: force,forcelj
    double precision :: fudge_factor
    double precision :: x1,y1,z1,x2,y2,z2
    double precision :: ffx,ffy,ffz
    double precision :: dx,dy,dz,dr,dr2,dr2i,dr6i,dr12i
    double precision :: potential_N2
    integer :: i,j,step,l
    integer :: itype,jtype,icell
    integer :: offset
    integer :: num1tmp,num2tmp,num3tmp
    character(len=3) :: char,char2

    potential_N2 = 0.0d0
    do i = 1 ,numMolArray(1)-1
       x1 = com(i)%x; y1 = com(i)%y; z1 = com(i)%z;
       do j=i+1,numMolArray(1)
          
          x2 = com(j)%x; y2 = com(j)%y; z2 = com(j)%z; 
          
          dx = x1-x2; dy = y1-y2; dz = z1-z2
          dx = dx-box*nint(dx*ibox)
          dy = dy-box*nint(dy*ibox)
          dz = dz-box*nint(dz*ibox)

          dr2 = dx*dx + dy*dy + dz*dz

     
          if(dr2.lt.denscut2)then
 

             dr2i = 1.0d0/dr2
             dr6i = dr2i*dr2i*dr2i
           
             potential_N2 = potential_N2 + dr6i*(dr6i-1.0d0)

          endif
       enddo
    enddo
    potential_N2 = 4.0d0*potential_N2
    print*
    print*,'WHAT IS FINAL POTENTIAL N2 DENS',potential_N2,box,hbox,ibox,denscut,denscut2
    print*

    write(1234,*)'WHAT IS FINAL POTENTIAL N2 DENS',potential_N2,mncellT,mrn
  end subroutine force_n2_dens



  subroutine force_vlist_dens
    implicit none
    double precision :: force,forcelj
    double precision :: x1,y1,z1,x2,y2,z2
    double precision :: dx,dy,dz,dr,dr2,dr2i,dr6i,dr12i
    double precision :: ffx,ffy,ffz
    double precision :: fudge_factor
    double precision :: potential_v2
    integer :: i,j,l,step
    integer :: itype,jtype,neigh
    integer :: tid,num
    integer :: offset,icell
    integer :: num1tmp,num2tmp,num3tmp
    character(len=3) :: char,char2


    potential_v2 = 0.0d0
    do i = 1 ,numMolArray(1)
       x1 = com(i)%x; y1 = com(i)%y; z1 = com(i)%z

       do j=1,denstrack(i)
         

          neigh = densvl(j+densvlalloc*(i-1))
          
            
          if(neigh.gt.i)then
             
             x2 = com(neigh)%x; y2 = com(neigh)%y; z2 = com(neigh)%z
             
             dx = x1-x2; dy = y1-y2; dz = z1-z2
             dx = dx-box*nint(dx*ibox)
             dy = dy-box*nint(dy*ibox)
             dz = dz-box*nint(dz*ibox)
             
             dr2 = dx*dx + dy*dy + dz*dz
             
             
             if(dr2.lt.denscut2)then
                
       
                dr2i = 1.0d0/dr2
                dr6i = dr2i*dr2i*dr2i
                potential_v2 = potential_v2 + dr6i*(dr6i-1.0d0)
         

             endif
          endif
       enddo
    enddo
    potential_v2 = potential_v2*4.0d0
    print*
    print*,'ENERGY VLIST DENS',potential_v2,box,ibox,hbox,denscut,denscut2
    print*

    write(1234,*)'ENERGY VLIST DENS',potential_v2,mncellT,mrn
  end subroutine force_vlist_dens
  
  
  subroutine init_list_cluster
    implicit none
    double precision :: scale
    integer :: i
    mrn         = box/int(box/mrcut)
    mncellT     = nint(vol/(mrn**3))
    mncellD     = nint(mncellT**(1.0d0/3.0d0))
    mncell_alloc = (mncellD+1)**3
    
    call print_init_list_cluster
    
  end subroutine init_list_cluster
  
  
  
  
  subroutine print_init_list_cluster
    implicit none
    write(1234,*)
    write(1234,*)'--------------- INTIALIZING CLUSTER LISTS-----------------------------------'
    write(1234,*)'           box               :',box
    write(1234,*)'           rn                :',mrn
    write(1234,*)'           ncellT            :',mncellT
    write(1234,*) '          ncellD            :',mncellD
    write(1234,*)'           nghbr list buff  :', mol_alloc
    write(1234,*)'           xyz        buff  :', xyz_alloc
    
    write(1234,*)'           MAKE SURE NCELLD CUBED IS NCELLT'
    write(1234,*)'--------------- END INTIALIZING CLUSTER LISTS-------------------------------'
    write(1234,*)
  end subroutine print_init_list_cluster
  



  end module mod_neighbor_cluster
