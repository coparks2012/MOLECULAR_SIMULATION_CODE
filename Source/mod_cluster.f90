module mod_cluster
  use global
  use mod_neighbor_cluster
  use mod_qcprot
  use mod_sort
  use mod_Grofile
  use mod_matrix_multiply
  use mod_prune
  use mod_build_ptta
   
  implicit none
  contains


    integer function get_greatestCluster(step,type,greatestBelongsto)
      implicit none
      integer :: step,type,ti
      integer :: numberofClusters
      integer :: greatestBelongsto
      integer :: T1,T2,clock_rate,clock_max
      integer :: i,numc,surface_molecules

      call system_clock(T1,clock_rate,clock_max)

      call build_ptta_wrap_cluster2
      !=============================================================
      ! use these when attempting RMSD
      !call prune
      !call prune_debug
      !call  xyzcenter_all(type)
      

      !=============================================================
      ! attempting this to resolve clustering issue with increasing time
      ! will see if problem hast to deal with integration or PBC


    

      call  xyzcenter_molecule(type)
      call  pbc_molcom(step,type)
      call  cluster_neighbor_wrapper(step,type)



      !call  cluster(step,type)

      numberofClusters = buildClusters()

      clusterhisto(:) = 0

      
      !! builds a histo of cluster sizes. THESE ONLY INCLUDE PARTICLES WITH HIGH DENSITY
      do ti = 1,numMolArray(1)
         if(belongsto(ti).ne.0) then
            clusterhisto(belongsto(ti)) = clusterhisto(belongsto(ti)) + 1
         endif
      enddo
      

      !! loop over particles to determine cluster which has the biggest size
      get_greatestCluster = 0
      do ti =1,numberofClusters
         if(clusterhisto(ti).gt.get_greatestCluster) then
            get_greatestCluster  = clusterhisto(ti)
            greatestBelongsto    = ti
         endif
      enddo


      !! APPEND SURFACE PARTICLES HERE
      call append_surface_molecules(greatestBelongsto,surface_molecules)
      get_greatestCluster = get_greatestCluster + surface_molecules


      call system_clock(T2,clock_rate,clock_max)
      time_clus = time_clus + real(T2-T1)/real(clock_rate)


 

      
    end function get_greatestCluster





    integer function clusterCriterium(ti)
      implicit none
      double precision :: dx,dy,dz,dr2
      double precision :: dcutoff
      integer :: ti,loop,count
      integer :: neighbor
             

      
      if(denstrack(ti).ge.densmin)then
         clusterCriterium = 1


      elseif(denstrack(ti).lt.densmin)then
         count = 0
         
      
         do loop=1,denstrack(ti)
            neighbor  = densvl( (ti-1)*densvlalloc + loop )
             
           if(denstrack(neighbor).ge.densmin)then
               dx = com(ti)%x-com(neighbor)%x
               dy = com(ti)%y-com(neighbor)%y
               dz = com(ti)%z-com(neighbor)%z
             
               dx = dx-box*nint(dx*ibox)
               dy = dy-box*nint(dy*ibox)
               dz = dz-box*nint(dz*ibox)
             
               dr2 = dx*dx + dy*dy +dz*dz
              
            

               if(dr2.lt.denscut2)then
                  count = count + 1
               endif
            endif
         enddo
         if(ti.eq.1094)then
            PRINT*,'ISOLATED CULPRIT',denstrack(ti)
            print*,'WHAT IS COUNT',count
         endif

         if(count.ge.1)then
            ClusterCriterium = 1
         else
            clusterCriterium = 0
         endif
      else
         clusterCriterium = 0
      endif
        


    end function clusterCriterium



    integer function clusterCriterium_surface(ti)
      implicit none
      double precision :: dx,dy,dz,dr2
      double precision :: dcutoff
      integer :: ti,loop,count
      integer :: neighbor
             
      clusterCriterium_surface = 0
      if(denstrack(ti).ge.densmin)then
         clusterCriterium_surface = 1
      endif
        


    end function clusterCriterium_surface

    subroutine append_surface_molecules(cluster,surface_number)
      implicit none
      double precision :: dx,dy,dz,dr2
      double precision :: dcutoff
      integer :: ti,loop,count,i,surface_number
      integer :: neighbor,cluster
             
      surfaceclus(:) = 0
      surface_number = 0
      do ti = 1,nummolarray(1)
      
         if(belongsto(ti).eq.cluster)then
            do loop=1,denstrack(ti)
               neighbor  = densvl( (ti-1)*densvlalloc + loop )
         
       
               if(denstrack(neighbor).lt.densmin)then
                  surfaceclus(neighbor) = surfaceclus(neighbor) + 1
               endif
            enddo
         endif
      enddo

      do i =1,nummolarray(1)
         if(surfaceclus(i).ge.1)then
            surface_number = surface_number + 1
            belongsto(i) = cluster
         endif
      enddo


    end subroutine append_surface_molecules

    



    recursive subroutine harvestCluster(ti ,numberofCluster)
      implicit none
      integer :: loop
      integer :: numberofCluster, ti
      integer :: neighbor
      
     
      do loop=1,denstrack(ti)
         
         neighbor  = densvl( (ti-1)*densvlalloc + loop )
         
         if( belongsto(neighbor).eq.0.and.clusterCriterium_surface(neighbor).eq.1) then
            
            belongsto(neighbor) = numberofCluster    
            call harvestCluster(neighbor,numberofCluster)
            
         endif
         
      enddo
    end subroutine harvestCluster



    integer function buildClusters()
      implicit none
      integer :: numberofCluster
      integer :: ti
      
      numberofCluster=0
      
      belongsto(:) = 0
      
      do ti =1,numMolArray(1)
         
         
         
         !! if a molecule is solid, but doesn't belong to a cluster yet
         if(belongsto(ti).eq.0.and.clusterCriterium_surface(ti).eq.1)then          
            numberofCluster = numberofCluster + 1
            belongsto(ti) = numberofCluster
            call  harvestCluster(ti ,numberofCluster)
         endif

         
      enddo
      
      buildClusters  = numberofCluster
    end function buildClusters
    
    

    subroutine cluster(step,type)
      implicit none
      double precision :: rotTEMP(3,numatomtemp)
      double precision :: local_template(3,numatomtemp)
      double precision :: rmsd,rot(9)
      double precision :: minalpha,minbeta,euc,dr2,dx,dy,dz
      double precision :: q0,q1,q2,q3
      integer :: i,j,k,imin,type,step,l,a,h,loc
      integer :: offset,length
      integer :: numalpha, numbeta


      print*,'======================================================='
      print*,'WE ARE NOW IN CLUSTER'
      
      open(unit = 1789, file = 'rmsd_debug.dat')
      open(unit = 1800, file = 'rmsd.dat')
      numalpha = 0
      numbeta  = 0
 
      do i = 1,numMolArray(type)

         minalpha = 1000.0d0
         minbeta  = 1000.0d0

         print*,denstrack(i)
         if(denstrack(i).ge.densmin)then
            call build_xyz(i,type)
         
            !---shift array xyz such that mol i has com {0,0,0}
            call center_array_xyz(type)
            
       
           
            !--- loop over number of templates for each polymorph
            do k = 1,numtemplate
               !---build local template
               local_template(:,:) = tempa(k)%temp(:,:)
               
               
               !---rotate template to minimize rmsd with first and center molecule
               call rotate_template(i,type,local_template,numatomtemp) 

           
               
               !---take xyz array which contains coordinates of all atoms
               !---find COM of all tagged molecules from those atoms
               !---comnumber determined here
               call build_xyz_com(i,type)
               
               
               !---take local templarte  array which contains coordinates of all atoms
               !---find COM of all molecules in template
               call build_localtemp_com(local_template,numatomtemp,type)
               
               
               !---for each coordinate in COM array, find closest neighboring coordinate in template COM
               call match(local_template,numatomtemp)

              
               
               length = uniquepairs*numatompertemp
               if(uniquepairs.gt.1)then
                  rmsd = calcrmsd(xyzf(:,1:length), tempf(:,1:length),length , rot,q0,q1,q2,q3)
               
                  !if(i.eq.109)then
                  !   PRINT*,'MADE IT IN'
                  !   if(k.eq.1)then
                  !      call xyzf_template_gro(xyzf,local_template,'cluster_debug1.gro')
                  !   elseif(k.eq.2)then
                  !      call xyzf_template_gro(xyzf,local_template,'cluster_debug2.gro')
                  !   elseif(k.eq.3)then
                  !      call xyzf_template_gro(xyzf,local_template,'cluster_debug3.gro')
                  !   elseif(k.eq.4)then
                  !      call xyzf_template_gro(xyzf,local_template,'cluster_debug4.gro')
                  !   endif
                  !endif
       
                  if(k.le.4)then
                     if(minalpha.gt.rmsd)then
                        minalpha = rmsd
                     endif
                  elseif(k.ge.5)then
                     if(minbeta.gt.rmsd)then
                        minbeta = rmsd
                     endif
                  endif                  
               endif
               
   
               
            enddo
       

            if(minalpha.lt.minbeta)then
               write(1789,*)'ALPHA',i,uniquepairs,denstrack(i),minalpha,minbeta
               numalpha = numalpha + 1
               minrmsd(i) = 0
            elseif(minbeta.lt.minalpha)then
               write(1789,*)'beta',i,uniquepairs,denstrack(i),minalpha,minbeta
               numbeta = numbeta + 1
               minrmsd(i) = 1
            else
               minrmsd(i) = 0
            endif
            write(1800,*)minalpha,minbeta
         endif 
      enddo
      print*,'num alpha and beta',numalpha,numbeta
      
    end subroutine cluster
    

    subroutine build_xyz_final
      implicit none
      double precision :: x(3,numatompertemp)
      integer :: neigh, i,j,ptr
      integer :: f1,f2
  
      
      do i = 1,comnumber
         neigh = sortcom(i)


         ptr   = int(xyzcom2match(2,neigh))

         f1 = (ptr-1)*numatompertemp+1
         f2 = f1 + numatompertemp-1

         call get_coords_xyz(x,neigh)

         xyzf(:,f1:f2) = x(:,:)
      end do
    end subroutine build_xyz_final


    subroutine get_coords_xyz(x,neigh)
      implicit none
      double precision :: x(3,numatompertemp)
      integer :: neigh
      integer :: i,j
      integer :: f1,f2
      
      f1 = (neigh-1)*numatompertemp+1
      f2 = f1 + numatompertemp-1
      
      x(:,:) = xyz(:,f1:f2)
    end subroutine get_coords_xyz

      
   

    !--------------------------------------------------------------!
    ! rotate template to match tagged molecule
    ! molecules and template must be centered
    subroutine rotate_template(molnum,moltype,local_template,len)
      implicit none
      integer :: molnum, moltype,len
      integer :: j,a
      double precision :: xsum,ysum,zsum
      double precision :: xx(3,numatompertemp), xtemp(3,numatompertemp)
      double precision :: rot(9)
      double precision :: local_template(3,len)
      double precision :: rmsd,q1,q2,q3,q4


      xx(:,:) = xyz(:,1:numatompertemp)
      

      !--- assumes first molecule is center molecule
      do j = 1,numatompertemp
         xtemp(:,j) = local_template(:,j)
      end do


      rmsd     = calcrmsd(xx, xtemp, numatompertemp, rot,q1,q2,q3,q4)
      call matrix_multiply(local_template,numatomtemp,rot)

     ! call xyz_template_centroid(local_template,numatomtemp,'both.gro')
     ! stop
      if(rmsd.gt.1.00.or.rmsd.ne.rmsd)then
         print*,'we have an error in rotate',rmsd,molnum
         call xyz_template_centroid(local_template,numatomtemp,'both.gro')
         call template_Gro(local_template,numatomtemp,'rotate.gro',0)
         call ptta_gro('ptta-fuckup.gro')
         stop
      endif


    end subroutine rotate_template





    !----------------------------------------------!
    ! xyz should be an array xyz(3,0:numneigh)
    ! xyz(:,0) should be reserved for tagged molecules
    subroutine build_xyz(molnumber,type)
      implicit none
      integer :: molnumber
      integer :: i,j,offset,type
      integer :: neighcom,neigh
      integer :: flag,count,count2


      !---check if array reallocation is needed
      if(mnum(molnumber)*numMolAtom(type).gt.xyz_alloc)then
         deallocate(xyz)
         xyz_alloc  = int(ceiling(1.2*mnum(molnumber)*numMolAtom(type)))
         allocate(xyz(3,xyz_alloc))
         print*,'WE HAVE REALLOCATED xyz',xyz_alloc
      endif

      

      !---make molnumber first element in array
      flag  = (molnumber-1)*numatompertemp+1
      count = 0
      do i = flag,flag + numatompertemp-1
         count = count +1
         xyz(1,count) = pruned(i)%x + pbcshift(1,molnumber)
         xyz(2,count) = pruned(i)%y + pbcshift(2,molnumber)
         xyz(3,count) = pruned(i)%z + pbcshift(3,molnumber)
      enddo

      !---track length of array
      xyz_number = numatompertemp + mnum(molnumber)*numatompertemp


      !---which part of neighbor list to loop over
      offset  = mol_alloc*(molnumber-1) 


      !---zero counter
      count2 = 0
      do i = 1,mnum(molnumber)
         
         !---number of COM neighbor
         neighcom = mvl(offset + i)
         !print*,'found these neighbors',neighcom
         !---starting atom for this COM neighbor
         neigh    = (neighcom-1)*numatompertemp+1
         
         do j = neigh , neigh+numatompertemp-1
            count2 = count2 + 1

            xyz(1,count+count2) = pruned(j)%x + pbcshift(1,neighcom) + shift(offset+i)%x
            xyz(2,count+count2) = pruned(j)%y + pbcshift(2,neighcom) + shift(offset+i)%y
            xyz(3,count+count2) = pruned(j)%z + pbcshift(3,neighcom) + shift(offset+i)%z
           

          enddo
      end do

    end subroutine build_xyz



    !----------------------------------------------!
    ! center array w/respect to centroid of tagged molecule (molecule 1)
    !***
    subroutine center_array_xyz(type)
      implicit none
      integer :: len,type
      double precision :: xx,yy,zz
      integer :: i

      xx = 0.0d0; yy = 0.0d0; zz =0.0d0
      
      do i = 1,numatompertemp
         xx = xx + xyz(1,i); yy = yy + xyz(2,i); zz = zz + xyz(3,i)
      end do
      
      xx = xx/real(numatompertemp)
      yy = yy/real(numatompertemp)
      zz = zz/real(numatompertemp)



      xyz(1,1:xyz_number) = xyz(1,1:xyz_number) - xx
      xyz(2,1:xyz_number) = xyz(2,1:xyz_number) - yy
      xyz(3,1:xyz_number) = xyz(3,1:xyz_number) - zz
    end subroutine center_array_xyz
      

    !----------------------------------------------!
    ! center array w/respect to centroid of tagged molecule (molecule 1)
    !***
    subroutine center_array_template(type)
      implicit none
      integer :: len,type
      integer :: i,j
      double precision :: xx,yy,zz

      do j = 1,numtemplate
         
         xx = 0.0d0; yy = 0.0d0; zz =0.0d0
         
         do i = 1,numatompertemp
            xx = xx + tempa(j)%temp(1,i)
            yy = yy + tempa(j)%temp(2,i)
            zz = zz + tempa(j)%temp(3,i)
         end do
         
         xx = xx/real(numatompertemp)
         yy = yy/real(numatompertemp)
         zz = zz/real(numatompertemp)
         
         tempa(j)%temp(1,:) = tempa(j)%temp(1,:) - xx
         tempa(j)%temp(2,:) = tempa(j)%temp(2,:) - yy
         tempa(j)%temp(3,:) = tempa(j)%temp(3,:) - zz


      end do
    end subroutine center_array_template





    subroutine  xyzcenter_one(type,mol,xcom,ycom,zcom)
      implicit none
      integer :: type,mol,num
      integer :: j
      integer :: f1,f2
      double precision :: xcom,ycom,zcom
      
      f1 = numatompertemp*(mol-1)+1
      f2 = f1 + (numatompertemp-1)
      
      !allocate the center of mass of each molecule
      xcom = 0.0d0
      ycom = 0.0d0
      zcom = 0.0d0
      do j = f1,f2
         xcom = xcom + pruned(j)%x
         ycom = ycom + pruned(j)%y
         zcom = zcom + pruned(j)%z
      end do
      xcom = xcom / real(numatompertemp)
      ycom = ycom / real(numatompertemp)
      zcom = zcom / real(numatompertemp)
      
      
    end subroutine xyzcenter_one
    
    subroutine build_xyz_com(molnumber,type)
      implicit none
      double precision :: xx,yy,zz
      integer :: molnumber
      integer :: i,j,type
      integer :: neigh
      integer :: offset
      integer :: atom
      
      !---check if array reallocation is needed
      if(mnum(molnumber).gt.xyzcom_alloc)then
         !deallocate(xyz,xyzcom2match,sortcom)

         deallocate(xyzcom2match,xyzcom)
         xyzcom_alloc  = int(ceiling(1.2*mnum(molnumber)))
         allocate(xyzcom2match(2,xyzcom_alloc),xyzcom(3,xyzcom_alloc))
         
         !allocate(xyzcom(3,xyzcom_alloc),xyzcom2match(2,xyzcom_alloc),sortcom(xyzcom_alloc))
         print*,'WE HAVE REALLOCATED',xyzcom_alloc
      endif
      
      comnumber = xyz_number/numatompertemp

      if(mod(xyz_number,numatompertemp).ne.0)then
         print*,'xyz number not divisible by numatompertemp'
         print*,'xyz number:',xyz_number
         print*,'numatompertemp:',numatompertemp
         stop
      endif

      do i =1,comnumber
         xx = 0.0d0; yy = 0.0d0; zz = 0.0d0

         do j = (i-1)*numatompertemp+1,(i-1)*numatompertemp+numatompertemp
            xx = xx + xyz(1,j)
            yy = yy + xyz(2,j)
            zz = zz + xyz(3,j)
         enddo

         xx = xx/real(numatompertemp)
         yy = yy/real(numatompertemp)
         zz = zz/real(numatompertemp)

         xyzcom(1,i) = xx
         xyzcom(2,i) = yy
         xyzcom(3,i) = zz
      enddo
         
            
    end subroutine build_xyz_com
      



   subroutine build_localtemp_com(local_template,len,type)
     implicit none
     double precision :: local_template(3,len)
     double precision :: xx,yy,zz
     integer :: len
     integer :: i,j,k
     integer :: type,localnummol

     do i =1,nummoltemp
        xx = 0.0d0; yy = 0.0d0; zz = 0.0d0
        
        do j = (i-1)*numatompertemp+1,(i-1)*numatompertemp+numatompertemp
           xx = xx + local_template(1,j)
           yy = yy + local_template(2,j)
           zz = zz + local_template(3,j)

        enddo
        
        xx = xx/real(numatompertemp)
        yy = yy/real(numatompertemp)
        zz = zz/real(numatompertemp)

        
        comtemp(1,i) = xx
        comtemp(2,i) = yy
        comtemp(3,i) = zz
   
     enddo
   
    end subroutine build_localtemp_com


    subroutine match(local_template,len)
      implicit none
      integer :: i,j,flag
      integer :: f1,f2,f3,f4,len
      double precision :: dx,dy,dz
      double precision :: x(3,numatompertemp)
      double precision :: local_template(3,len)
      double precision :: x1,y1,z1
      double precision :: x2,y2,z2
      double precision :: dr2,min
      double precision :: neigh
      double precision :: dis1,dis2

      uniquepairs = 0

      !---------------------------------------------------------
      ! loop over all molecules in template
      ! find molecule in xyz array that is closest to each molecule in template
      do i =1,nummoltemp
         x1 = comtemp(1,i); y1 = comtemp(2,i); z1 = comtemp(3,i)

         min = 1000.0d0
         do j = 1,comnumber
            x2 = xyzcom(1,j); y2 = xyzcom(2,j); z2 = xyzcom(3,j)
            dx = x1-x2; dy = y1-y2; dz = z1-z2
            dr2 = dx*dx + dy*dy + dz*dz
                        

            if(dr2.lt.min)then
               min = dr2
               xyzcom2match(1,i) = dr2
               xyzcom2match(2,i) = real(j)

            endif
         enddo
      enddo

      !---only pair molecule and template if no other molecule in template is closer to xyz molecule
      uniquepairs = 0
      do i =1,nummoltemp
         dis1  = xyzcom2match(1,i)
         neigh = xyzcom2match(2,i)
         
         flag = 0
         do j = 1,nummoltemp
            if(i.eq.j)cycle

            if(neigh.eq.xyzcom2match(2,j))then
               if(dis1.gt.xyzcom2match(1,j))then
                  flag = flag + 1
               endif
            endif
         enddo
         
         if(flag.eq.0)then


            uniquepairs = uniquepairs + 1
            !---get solution pointers
            f1 = (uniquepairs-1)*numatompertemp +1 
            f2 = (uniquepairs-1)*numatompertemp+numatompertemp



            call get_coords_xyz(x,int(neigh))

            xyzf(:,f1:f2)  = x(:,:)
 

            f3 = (i-1)*numatompertemp + 1
            f4 = (i-1)*numatompertemp + numatompertemp
            tempf(:,f1:f2) = local_template(:,f3:f4)

     
         endif
      enddo
  
     
    end subroutine match


    subroutine  xyzcenter_all(type)
      implicit none
      integer :: type,mol,num
      integer :: j
      integer :: f1,f2
      double precision :: xcom,ycom,zcom
      

      do mol = 1,numMolArray(type)
         f1 = numatompertemp*(mol-1)+1
         f2 = f1 + (numatompertemp-1)
         
         !determinr the center of mass of each molecule
         xcom = 0.0d0
         ycom = 0.0d0
         zcom = 0.0d0
         do j = f1,f2
            xcom = xcom + pruned(j)%x
            ycom = ycom + pruned(j)%y
            zcom = zcom + pruned(j)%z
         end do
        

         com(mol)%x = xcom / real(numatompertemp)
         com(mol)%y = ycom / real(numatompertemp)
         com(mol)%z = zcom / real(numatompertemp)
      end do

      open(unit = 99999,file = 'com_atom.dat')
      do j = 1,nummolarray(type)
         write(99999,*)com(j)%x,com(j)%y,com(j)%z
      enddo
      close(99999)
      !stop


    end subroutine xyzcenter_all

    subroutine  xyzcenter_molecule(type)
      implicit none
      integer :: type,mol,num
      integer :: j
      integer :: f1,f2,atomnum
      double precision :: xcom,ycom,zcom
     
      atomnum = nummolatom(1)
      do mol = 1,numMolArray(type)
         f1 = atomnum*(mol-1)+1
         f2 = f1 + (atomnum-1)
         
         !determinr the center of mass of each molecule
         xcom = 0.0d0
         ycom = 0.0d0
         zcom = 0.0d0
         do j = f1,f2
            xcom = xcom + ptta_full(j)%x
            ycom = ycom + ptta_full(j)%y
            zcom = zcom + ptta_full(j)%z

         end do

  
         com(mol)%x = xcom / real(atomnum)
         com(mol)%y = ycom / real(atomnum)
         com(mol)%z = zcom / real(atomnum)
         
      end do




    end subroutine xyzcenter_molecule


    

    subroutine pbc_molcom(step,type)
      implicit none
      integer :: type
      integer :: i,step

      pbcshift(:,:) = 0.0d0
      
      do i = 1,numMolArray(type)

         do while (com(i)%x.gt.box)
            com(i)%x = com(i)%x - box 
            pbcshift(1,i) = pbcshift(1,i) - box
         end do
         do while  (com(i)%x .lt.0.d0)               
            com(i)%x = com(i)%x + box   
            pbcshift(1,i) = pbcshift(1,i) + box
         end do
         
         
         do while (com(i)%y.gt.box) 
            com(i)%y = com(i)%y - box
            pbcshift(2,i) = pbcshift(2,i) - box
         enddo
         do while  (com(i)%y.lt.0.d0) 
            com(i)%y = com(i)%y + box    
            pbcshift(2,i) = pbcshift(2,i) + box 
         enddo
         

         do while (com(i)%z.gt.box) 
            com(i)%z = com(i)%z - box 
            pbcshift(3,i) = pbcshift(2,i) - box
         enddo
         do while  (com(i)%z.lt.0.d0) 
            com(i)%z = com(i)%z + box    
            pbcshift(3,i) = pbcshift(3,i) + box
         enddo
         
      enddo

   !   call GroFile_COM('com.gro',1)
   !   stop

    end subroutine pbc_molcom

    subroutine pbc_molcom_triclinic
      implicit none
      double precision :: dx,dy,dz
      double precision :: box_x,box_y,box_z
      double precision :: x1,y1,z1
      double precision :: xbound_lo,xbound_hi
      integer :: i
      
      
      print*,'boxes',xprd,yprd,zprd
      print*,'xy xz yz',xy,xz,yz
      print*,'x',xlo,xhi
      print*,'y',ylo,yhi
      print*,'z',zlo,zhi
      
      do i = 1,nummolarray(1)
         x1 =com(i)%x
         y1 =com(i)%y
         z1 =com(i)%z
         
         
         
         if (z1.lt.zlo)then 
            z1 = z1 + zprd
            y1 = y1 + yz
            x1 = x1 + xz
            
            pbcshift(3,i) = pbcshift(1,i) + zprd
            pbcshift(2,i) = pbcshift(1,i) + yz
            pbcshift(1,i) = pbcshift(1,i) + xz
            
            
         elseif(z1.gt.zhi)then
            z1 = z1 - zprd
            y1 = y1 - yz
            x1 = x1 - xz
            
            pbcshift(3,i) = pbcshift(1,i) - zprd
            pbcshift(2,i) = pbcshift(1,i) - yz
            pbcshift(1,i) = pbcshift(1,i) - xz
            
         endif
         
         if (y1.lt.ylo)then
            y1 = y1 + yprd
            x1 = x1 + xy
            
            pbcshift(2,i) = pbcshift(1,i) + yprd
            pbcshift(1,i) = pbcshift(1,i) + xy
            
         elseif(y1.gt.yhi)then
            y1 = y1 - yprd
            x1 = x1 - xy
            
            pbcshift(2,i) = pbcshift(1,i) - yprd
            pbcshift(1,i) = pbcshift(1,i) - xy
            
         endif
         
         xbound_hi = (xy/yprd)*y1+xprd
         xbound_lo = (xy/yprd)*y1
         
       if(x1.lt.xbound_lo)then
          x1 = x1 + xprd
          
          pbcshift(1,i) = pbcshift(1,i) + xprd
          
       elseif(x1.gt.xbound_hi)then
          x1 = x1 -xprd
          
          pbcshift(1,i) = pbcshift(1,i) - xprd
          
       endif
       
       
       com(i)%x=x1
       com(i)%y=y1
       com(i)%z=z1
       call GroFile_COM('com_tri.gro',1)

    enddo
  end subroutine pbc_molcom_triclinic

    subroutine print_cluster_size(w1,greatest)
      implicit none
      integer :: greatest,w1
      w1 = get_greatestCluster(0,1,greatest)
      print*,'NUMBER IN CLUSTER',w1
    end subroutine print_cluster_size



  
  end module mod_cluster
