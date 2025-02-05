module mod_force_nonbond
  use global
  use omp_lib
  

  contains
    !dir$ attributes offload:mic::lj_cut_coul_dsf_nonewton_MIC
    subroutine lj_cut_coul_dsf_nonewton_MIC(step)
       implicit none
      double precision :: force,forcelj,forcecoul
      double precision :: x1,y1,z1,x2,y2,z2
      double precision :: dx,dy,dz,dr,dr2,dr2i,dr6i,dr12i,dri
      double precision :: ffx,ffy,ffz
      double precision :: qtmp,r,prefactor,erfcc,erfcd,t
      double precision :: boxdx,boxdy,boxdz
      double precision :: bondscalef,bondscale
      integer :: i,j,l,step
      integer :: itype,jtype,neigh
      integer :: tid,num
      integer :: offset,ioffset,neigh_off
      integer :: T1,T2,clock_rate,clock_max
      
      !call system_clock(T1,clock_rate,clock_max)
      
     
      !$omp parallel do schedule(dynamic) reduction(+:potential,e_coul,ffx,ffy,ffz) default(firstprivate),&
      !$omp& shared(position,ff,nlist,numneigh,q,lj1,lj2,lj3,lj4)
      do i = 1 ,np
         x1 = position(i)%x; y1 = position(i)%y; z1 = position(i)%z; itype = position(i)%type
         qtmp = q(i)
         ioffset = (itype-1)*numAtomType
         neigh_off = neigh_alloc*(i-1)
         num = numneigh(i)
         ffx = 0.0d0; ffy = 0.0d0; ffz = 0.0d0   
         
         !dir$ vector aligned
         !dir$ simd reduction(+:potential,e_coul,ffx,ffy,ffz)
         do j= 1,num
            
            neigh      =  nlist(neigh_off+j)
            
            dx = x1-position(neigh)%x
            dy = y1-position(neigh)%y
            dz = z1-position(neigh)%z
            jtype = position(neigh)%type
            
            boxdx = dx*ibox; boxdy = dy*ibox; boxdz = dz*ibox
            boxdx = (boxdx+sign(1/(epsilon(boxdx)),boxdx)) -sign(1/epsilon(boxdx),dx)
            boxdy = (boxdy+sign(1/(epsilon(boxdy)),boxdy)) -sign(1/epsilon(boxdy),dy)
            boxdz = (boxdz+sign(1/(epsilon(boxdz)),boxdz)) -sign(1/epsilon(boxdz),dz)
            
            dx = dx-box*boxdx; dy = dy-box*boxdy; dz = dz-box*boxdz
            dr2 = dx*dx + dy*dy + dz*dz
            
            !---lennard jones interactions
            dr2i = 1.0d0/dr2
            dr6i = dr2i*dr2i*dr2i
            if(dr2.gt.rcut2)dr6i=0.0d0
            offset = ioffset + jtype
            forcelj = dr6i*(lj1(offset)*dr6i-lj2(offset))
            potential = potential + dr6i*(dr6i*lj3(offset)-lj4(offset))      
            
            
            !---electrostatic calculations     
            r = sqrt(dr2)
            dri = 1.0d0/r
            prefactor =  coulpre*qtmp*q(neigh)*dri
            !prefactor = qtmp*q(neigh)*dri
            if(dr2.gt.cut_coulsq)prefactor =0.0d0
            erfcd = exp(-alpha*alpha*r*r)
            t = 1.0 / (1.0 + EWALD_P*alpha*r)
            erfcc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * erfcd
            forcecoul = prefactor * (erfcc*dri + 2.0*alpha*MY_PIS_INV*erfcd +&
                 r*f_shift) * r
            
            e_coul = e_coul + prefactor*(erfcc-r*e_shift-dr2*f_shift)             
            
            force   = (forcecoul+forcelj)*dr2i
            
    
          


            ffx = ffx + dx*force
            ffy = ffy + dy*force
            ffz = ffz + dz*force
            
         enddo
         ff(i)%x = ffx; ff(i)%y = ffy; ff(i)%z = ffz
      enddo
      !$omp end parallel do 

      !call system_clock(T2,clock_rate,clock_max)
      
      !time_nonbond = time_nonbond + real(T2-T1)/real(clock_rate)
      potential = 0.50d0*potential
      e_coul = 0.50d0*e_coul


      
    end subroutine lj_cut_coul_dsf_nonewton_MIC
    
    
    
    subroutine lj_cut_coul_dsf_nonewton_NOMIC(step)
      implicit none
      double precision :: force,forcelj,forcecoul
      double precision :: x1,y1,z1,x2,y2,z2
      double precision :: dx,dy,dz,dr,dr2,dr2i,dr6i,dr12i,dri
      double precision :: ffx,ffy,ffz
      double precision :: qtmp,r,prefactor,erfcc,erfcd,t
      double precision :: boxdx,boxdy,boxdz
      double precision :: bondscalef,bondscale
      integer :: i,j,l,step
      integer :: itype,jtype,neigh
      integer :: tid,num
      integer :: offset,ioffset,neigh_off
      integer :: T1,T2,clock_rate,clock_max
      
      !call system_clock(T1,clock_rate,clock_max)
      
     
      !$omp parallel do schedule(dynamic) reduction(+:potential,e_coul,ffx,ffy,ffz) default(firstprivate),&
      !$omp& shared(position,ff,nlist,numneigh,q,lj1,lj2,lj3,lj4)
      do i = 1 ,np
         x1 = position(i)%x; y1 = position(i)%y; z1 = position(i)%z; itype = position(i)%type
         qtmp = q(i)
         ioffset = (itype-1)*numAtomType
         neigh_off = neigh_alloc*(i-1)
         num = numneigh(i)
         ffx = 0.0d0; ffy = 0.0d0; ffz = 0.0d0   
         
         !dir$ vector aligned
         !dir$ simd reduction(+:potential,e_coul,ffx,ffy,ffz)
         do j= 1,num
            
            neigh      =  nlist(neigh_off+j)
            
            dx = x1-position(neigh)%x
            dy = y1-position(neigh)%y
            dz = z1-position(neigh)%z
            jtype = position(neigh)%type
            
            boxdx = dx*ibox; boxdy = dy*ibox; boxdz = dz*ibox
            boxdx = (boxdx+sign(1/(epsilon(boxdx)),boxdx)) -sign(1/epsilon(boxdx),dx)
            boxdy = (boxdy+sign(1/(epsilon(boxdy)),boxdy)) -sign(1/epsilon(boxdy),dy)
            boxdz = (boxdz+sign(1/(epsilon(boxdz)),boxdz)) -sign(1/epsilon(boxdz),dz)
            
            dx = dx-box*boxdx; dy = dy-box*boxdy; dz = dz-box*boxdz
            dr2 = dx*dx + dy*dy + dz*dz
            
            !---lennard jones interactions
            dr2i = 1.0d0/dr2
            dr6i = dr2i*dr2i*dr2i
            if(dr2.gt.rcut2)dr6i=0.0d0
            offset = ioffset + jtype
            forcelj = dr6i*(lj1(offset)*dr6i-lj2(offset))
            potential = potential + dr6i*(dr6i*lj3(offset)-lj4(offset))      
            
            
            !---electrostatic calculations     
            r = sqrt(dr2)
            dri = 1.0d0/r
            prefactor =  coulpre*qtmp*q(neigh)*dri
            !prefactor = qtmp*q(neigh)*dri
            if(dr2.gt.cut_coulsq)prefactor =0.0d0
            erfcd = exp(-alpha*alpha*r*r)
            t = 1.0 / (1.0 + EWALD_P*alpha*r)
            erfcc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * erfcd
            forcecoul = prefactor * (erfcc*dri + 2.0*alpha*MY_PIS_INV*erfcd +&
                 r*f_shift) * r
            
            e_coul = e_coul + prefactor*(erfcc-r*e_shift-dr2*f_shift)             
            
            force   = (forcecoul+forcelj)*dr2i
            
    
          


            ffx = ffx + dx*force
            ffy = ffy + dy*force
            ffz = ffz + dz*force
            
         enddo
         ff(i)%x = ffx; ff(i)%y = ffy; ff(i)%z = ffz
      enddo
      !$omp end parallel do 

      !call system_clock(T2,clock_rate,clock_max)
      
      !time_nonbond = time_nonbond + real(T2-T1)/real(clock_rate)
      potential = 0.50d0*potential
      e_coul = 0.50d0*e_coul
      
    end subroutine lj_cut_coul_dsf_nonewton_NOMIC
    
    
    
    !dir$ attributes offload:mic::lj_cut_coul_long_nonewton_MIC
    subroutine lj_cut_coul_long_nonewton_MIC
      implicit none
      real*4  :: force,forcelj,forcecoul
      real*4  :: x1,y1,z1,x2,y2,z2
      real*4  :: dx,dy,dz,dr,dr2,dr2i,dr6i,dr12i,dri
      double precision :: ffx,ffy,ffz,vvx,vvy,vvz
      real*4  :: qtmp,r,prefactor,erfc,t
      real*4  :: boxdx,boxdy,boxdz
      real*4  :: grij,expm2
      integer :: i,j,l,step
      integer :: itype,jtype,neigh
      integer :: tid,num
      integer :: offset,ioffset,neigh_off
      integer :: T1,T2,clock_rate,clock_max

                                         
      vvx = 0.0d0
      vvy = 0.0d0
      vvz = 0.0d0
      !$omp parallel do schedule(dynamic) reduction(+:potential,e_coul,ffx,ffy,ffz,vvx,vvy,vvz) default(firstprivate),&
      !$omp& shared(position,ff,nlist,numneigh,q,lj1,lj2,lj3,lj4)
         do i = 1 ,np
         x1 = position(i)%x; y1 = position(i)%y; z1 = position(i)%z; itype = position(i)%type
         qtmp = q(i)
         ioffset = (itype-1)*numAtomType
         neigh_off = neigh_alloc*(i-1)
         num = numneigh(i)
         
         ffx = 0.0d0; ffy = 0.0d0; ffz = 0.0d0   
         
         !dir$ vector aligned
         !dir$ simd reduction(+:potential,e_coul,ffx,ffy,ffz,vvx,vvy,vvz)
         do j= 1,num
            
            neigh           = nlist(neigh_off+j)
            
            dx = x1-position(neigh)%x
            dy = y1-position(neigh)%y
            dz = z1-position(neigh)%z
            jtype = position(neigh)%type
            
            boxdx = dx*ibox; boxdy = dy*ibox; boxdz = dz*ibox
            boxdx = (boxdx+sign(1/(epsilon(boxdx)),boxdx)) -sign(1/epsilon(boxdx),dx)
            boxdy = (boxdy+sign(1/(epsilon(boxdy)),boxdy)) -sign(1/epsilon(boxdy),dy)
            boxdz = (boxdz+sign(1/(epsilon(boxdz)),boxdz)) -sign(1/epsilon(boxdz),dz)
            
            dx = dx-box*boxdx; dy = dy-box*boxdy; dz = dz-box*boxdz
            dr2 = dx*dx + dy*dy + dz*dz
            
            !---lennard jones interactions
            dr2i = 1.0d0/dr2
            dr6i = dr2i*dr2i*dr2i
            if(dr2.gt.rcut2)dr6i = 0.0d0
            offset = ioffset + jtype
            forcelj = dr6i*(lj1(offset)*dr6i-lj2(offset))
            potential = potential + dr6i*(dr6i*lj3(offset)-lj4(offset))  


            !===electrostatic calculation
            r = sqrt(dr2)
            prefactor =  coulpre*qtmp*q(neigh)/r
            grij = gewald*r
            expm2 = exp(-grij*grij)
            t = 1.0/(1.0 + EWALD_P*grij)
            erfc = t*(A1+t*(A2+t*(A3+t*(A4+t*A5))))*expm2
            if(dr2.gt.cut_coulsq)prefactor = 0.0d0
            forcecoul = prefactor * (erfc + EWALD_F*grij*expm2)
            e_coul = e_coul + prefactor*erfc    


            
            force   = (forcecoul+forcelj)*dr2i
            ffx = ffx + dx*force
            ffy = ffy + dy*force
            ffz = ffz + dz*force
            
            vvx = vvx + dx*dx*force
            vvy = vvy + dy*dy*force
            vvz = vvz + dz*dz*force
            
 
         enddo
         ff(i)%x = ffx; ff(i)%y = ffy; ff(i)%z = ffz

      enddo
      !$omp end parallel do

      virialx = 0.50*vvx
      virialy = 0.50*vvy
      virialz = 0.50*vvz
      potential = 0.50d0*potential
      e_coul = 0.50d0*e_coul

     

    end subroutine lj_cut_coul_long_nonewton_MIC
    
    
    
    subroutine lj_cut_coul_long_nonewton_NOMIC
      implicit none
      real*4 :: force,forcelj,forcecoul
      real*4  :: x1,y1,z1,x2,y2,z2
      real*4  :: dx,dy,dz,dr,dr2,dr2i,dr6i,dr12i,dri
      real*4 :: grij,expm2
      double precision :: ffx,ffy,ffz,vvx,vvy,vvz
      real*4 :: qtmp,r,prefactor,erfc,t
      real*4 :: boxdx,boxdy,boxdz
      integer :: i,j,l,step
      integer :: itype,jtype,neigh
      integer :: tid,num
      integer :: offset,ioffset,neigh_off
      integer :: T1,T2,clock_rate,clock_max
      integer :: numt,threadid

      vvx  =0.0d0
      vvy  =0.0d0
      vvz  =0.0d0

      !$omp parallel do schedule(dynamic) reduction(+:potential,e_coul,ffx,ffy,ffz,vvx,vvy,vvz) default(firstprivate),&
      !$omp& shared(position,ff,nlist,numneigh,q,lj1,lj2,lj3,lj4)
      do i = 1 ,np
         x1 = position(i)%x; y1 = position(i)%y; z1 = position(i)%z; itype = position(i)%type
         qtmp = q(i)
         ioffset = (itype-1)*numAtomType
         neigh_off = neigh_alloc*(i-1)
         num = numneigh(i)
         
         ffx = 0.0d0; ffy = 0.0d0; ffz = 0.0d0   
         
         !dir$ vector aligned
         !dir$ simd reduction(+:potential,e_coul,ffx,ffy,ffz,vvx,vvy,vvz)
         do j= 1,num
            
            neigh           = nlist(neigh_off+j)
            
            dx = x1-position(neigh)%x
            dy = y1-position(neigh)%y
            dz = z1-position(neigh)%z
            jtype = position(neigh)%type
            
            boxdx = dx*ibox; boxdy = dy*ibox; boxdz = dz*ibox
            boxdx = (boxdx+sign(1/(epsilon(boxdx)),boxdx)) -sign(1/epsilon(boxdx),dx)
            boxdy = (boxdy+sign(1/(epsilon(boxdy)),boxdy)) -sign(1/epsilon(boxdy),dy)
            boxdz = (boxdz+sign(1/(epsilon(boxdz)),boxdz)) -sign(1/epsilon(boxdz),dz)
            
            dx = dx-box*boxdx; dy = dy-box*boxdy; dz = dz-box*boxdz
            dr2 = dx*dx + dy*dy + dz*dz
            
            !---lennard jones interactions
            dr2i = 1.0d0/dr2
            dr6i = dr2i*dr2i*dr2i
            if(dr2.gt.rcut2)dr6i = 0.0d0
            offset = ioffset + jtype
            forcelj = dr6i*(lj1(offset)*dr6i-lj2(offset))
            potential = potential + dr6i*(dr6i*lj3(offset)-lj4(offset))  


            !===electrostatic calculation
            r = sqrt(dr2)
            prefactor =  coulpre*qtmp*q(neigh)/r
            grij = gewald*r
            expm2 = exp(-grij*grij)
            t = 1.0/(1.0 + EWALD_P*grij)
            erfc = t*(A1+t*(A2+t*(A3+t*(A4+t*A5))))*expm2
            if(dr2.gt.cut_coulsq)prefactor = 0.0d0
            forcecoul = prefactor * (erfc + EWALD_F*grij*expm2)
            e_coul = e_coul + prefactor*erfc    

            
            force   = (forcecoul+forcelj)*dr2i
            ffx = ffx + dx*force
            ffy = ffy + dy*force
            ffz = ffz + dz*force
            
            vvx = vvx + dx*dx*force
            vvy = vvy + dy*dy*force
            vvz = vvz + dz*dz*force
            
 
         enddo
         ff(i)%x = ffx; ff(i)%y = ffy; ff(i)%z = ffz

      enddo
      !$omp end parallel do


      virialx = 0.50*vvx
      virialy = 0.50*vvy
      virialz = 0.50*vvz
      


      
      !time_nonbond = time_nonbond + real(T2-T1)/real(clock_rate)
      potential = 0.50d0*potential
      e_coul = 0.50d0*e_coul
      

  end subroutine lj_cut_coul_long_nonewton_NOMIC



  



  
  subroutine lj_cut_coul_dsf_N2(step)
    implicit none
    double precision  :: force,forcelj,forcecoul
    double precision  :: x1,y1,z1,x2,y2,z2
    double precision  :: dx,dy,dz,dr,dr2,dr2i,dr6i,dr12i,dri
    double precision  :: bondscale,bondscalef
    double precision  :: grij,expm2
    double precision  :: fudge_factor
    double precision :: ffx,ffy,ffz
    double precision :: qtmp,r,prefactor,erfc,t
    double precision :: boxdx,boxdy,boxdz,erfcd,erfcc
    integer :: i,j,l,step
    integer :: itype,jtype,neigh
    integer :: tid,num
    integer :: offset,ioffset,neigh_off
    integer :: T1,T2,clock_rate,clock_max
    integer :: numt,threadid
    integer :: num1tmp,num2tmp,num3tmp
      

    do i = 1 ,np-1
       x1 = position(i)%x; y1 = position(i)%y; z1 = position(i)%z; itype = position(i)%type
       num1tmp = num1bond(i)
       num2tmp = num2bond(i)
       num3tmp = num3bond(i)
       qtmp = q(i)
       
       do j=i+1,np
          
          x2 = position(j)%x; y2 = position(j)%y; z2 = position(j)%z; jtype = position(j)%type
          
          dx = x1-x2; dy = y1-y2; dz = z1-z2
          dx = dx-box*nint(dx*ibox)
          dy = dy-box*nint(dy*ibox)
          dz = dz-box*nint(dz*ibox)

          dr2 = dx*dx + dy*dy + dz*dz

     
          if(dr2.lt.rcut2)then

             fudge_factor = 1.0d0
             do l = 1,num3tmp
                if(specbond(l,i).eq.j.and.l.le.num2tmp)then
                   fudge_factor = 0.0d0
                elseif(specbond(l,i).eq.j.and.l.gt.num2tmp)then
                   fudge_factor = 0.50d0
                endif
             enddo
             

             dr2i = 1.0d0/dr2
             dr6i = dr2i*dr2i*dr2i
           
             offset  = numAtomType*(itype-1)+jtype
             forcelj = dr6i*(lj1(offset)*dr6i-lj2(offset))
             potential = potential + fudge_factor*dr6i*(lj3(offset)*dr6i-lj4(offset))
             

             !---electrostatic calculations     
             r = sqrt(dr2)
             dri = 1.0d0/r
             prefactor =  coulpre*qtmp*q(j)*dri
             if(dr2.gt.cut_coulsq)prefactor =0.0d0
             erfcd = exp(-alpha*alpha*r*r)
             t = 1.0 / (1.0 + EWALD_P*alpha*r)
             erfcc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * erfcd
             forcecoul = prefactor * (erfcc*dri + 2.0*alpha*MY_PIS_INV*erfcd +&
                  r*f_shift) * r
             
             e_coul = e_coul + fudge_factor*prefactor*(erfcc-r*e_shift-dr2*f_shift)             
             
             
             force   = fudge_factor*(forcecoul+forcelj)*dr2i
             

             ff(i)%x = ff(i)%x + dx*force
             ff(i)%y = ff(i)%y + dy*force
             ff(i)%z = ff(i)%z + dz*force


             ff(j)%x = ff(j)%x - dx*force
             ff(j)%y = ff(j)%y - dy*force
             ff(j)%z = ff(j)%z - dz*force




          endif
       enddo
    enddo


  end subroutine lj_cut_coul_dsf_N2

  subroutine force_N2
    implicit none
    double precision  :: force,forcelj,forcecoul
    double precision  :: x1,y1,z1,x2,y2,z2
    double precision  :: dx,dy,dz,dr,dr2,dr2i,dr6i,dr12i,dri
    double precision  :: bondscale,bondscalef
    double precision  :: grij,expm2
    double precision  :: fudge_factor
    double precision :: ffx,ffy,ffz
    double precision :: qtmp,r,prefactor,erfc,t
    double precision :: boxdx,boxdy,boxdz
    integer :: i,j,l,step
    integer :: itype,jtype,neigh
    integer :: tid,num
    integer :: offset,ioffset,neigh_off
    integer :: T1,T2,clock_rate,clock_max
    integer :: numt,threadid,bondflag
    integer :: num1tmp,num2tmp,num3tmp
      

    do i = 1 ,np-1
       x1 = position(i)%x; y1 = position(i)%y; z1 = position(i)%z; itype = position(i)%type
       num1tmp = num1bond(i)
       num2tmp = num2bond(i)
       num3tmp = num3bond(i)
       qtmp = q(i)
       
       do j=i+1,np
          
          x2 = position(j)%x; y2 = position(j)%y; z2 = position(j)%z; jtype = position(j)%type
          
          dx = x1-x2; dy = y1-y2; dz = z1-z2
          dx = dx-box*nint(dx*ibox)
          dy = dy-box*nint(dy*ibox)
          dz = dz-box*nint(dz*ibox)

          dr2 = dx*dx + dy*dy + dz*dz

     
          if(dr2.lt.rcut2)then

             fudge_factor = 1.0d0
             bondflag = 1
             do l = 1,num3tmp
                if(specbond(l,i).eq.j.and.l.le.num2tmp)then
                   fudge_factor = 0.0d0
                   bondflag = 0
                elseif(specbond(l,i).eq.j.and.l.gt.num2tmp)then
                   fudge_factor = 0.50d0
                   bondflag = 0
                endif
             enddo
             
             forcelj = 0.0d0
             dr2i = 1.0d0/dr2
             dr6i = dr2i*dr2i*dr2i
             if(dr2.gt.rcut2)dr6i =0.0d0
             offset  = numAtomType*(itype-1)+jtype
             forcelj = dr6i*(lj1(offset)*dr6i-lj2(offset))
             potential = potential + fudge_factor*dr6i*(lj3(offset)*dr6i-lj4(offset))
  
           
             !***check to see if force correction is correct
      
             forcecoul =0.0d0
             r = sqrt(dr2)
             prefactor =  coulpre*qtmp*q(j)/r
             grij = gewald*r
             expm2 = exp(-grij*grij)
             t = 1.0/(1.0 + EWALD_P*grij)
             erfc = t*(A1+t*(A2+t*(A3+t*(A4+t*A5))))*expm2
             if(dr2.gt.cut_coulsq)prefactor = 0.0d0
             forcecoul = prefactor * (erfc + EWALD_F*grij*expm2)
             forcecoul = forcecoul - (1.0-fudge_factor)*prefactor
             
             e_coul = e_coul + prefactor*erfc    
             e_coul = e_coul - (1.0-fudge_factor)*prefactor
             

             force   = (forcecoul+fudge_factor*forcelj)*dr2i
             

             ff(i)%x = ff(i)%x + dx*force
             ff(i)%y = ff(i)%y + dy*force
             ff(i)%z = ff(i)%z + dz*force


             ff(j)%x = ff(j)%x - dx*force
             ff(j)%y = ff(j)%y - dy*force
             ff(j)%z = ff(j)%z - dz*force




          endif
       enddo
    enddo

  end subroutine force_N2



   subroutine force_N2_14
    implicit none
    double precision  :: force,forcelj,forcecoul
    double precision  :: x1,y1,z1,x2,y2,z2
    double precision  :: dx,dy,dz,dr,dr2,dr2i,dr6i,dr12i,dri
    double precision  :: bondscale,bondscalef
    double precision  :: grij,expm2
    double precision  :: fudge_factor
    double precision :: ffx,ffy,ffz
    double precision :: qtmp,r,prefactor,erfc,t
    double precision :: boxdx,boxdy,boxdz
    integer :: i,j,l,step
    integer :: itype,jtype,neigh
    integer :: tid,num
    integer :: offset,ioffset,neigh_off
    integer :: T1,T2,clock_rate,clock_max
    integer :: numt,threadid,bondflag
    integer :: num1tmp,num2tmp,num3tmp
      

    do i = 1 ,np-1
       x1 = position(i)%x; y1 = position(i)%y; z1 = position(i)%z; itype = position(i)%type
       num1tmp = num1bond(i)
       num2tmp = num2bond(i)
       num3tmp = num3bond(i)
       qtmp = q(i)
       
       do j=i+1,np
          
          x2 = position(j)%x; y2 = position(j)%y; z2 = position(j)%z; jtype = position(j)%type
          
          dx = x1-x2; dy = y1-y2; dz = z1-z2
          dx = dx-box*nint(dx*ibox)
          dy = dy-box*nint(dy*ibox)
          dz = dz-box*nint(dz*ibox)

          dr2 = dx*dx + dy*dy + dz*dz

    
          fudge_factor = 1.0d0
          bondflag = 1
          do l = 1,num3tmp
             if(specbond(l,i).eq.j.and.l.le.num2tmp)then
                fudge_factor = 0.0d0
                bondflag = 0
             elseif(specbond(l,i).eq.j.and.l.gt.num2tmp)then
                fudge_factor = 0.50d0
                bondflag = 0
             endif
          enddo
          
          
          forcecoul =0.0d0
          forcelj = 0.0d0
          dr2i = 0.0d0
          if(bondflag.eq.1)then
             dr2i = 1.0d0/dr2
             dr6i = dr2i*dr2i*dr2i
             if(dr2.gt.rcut2)dr6i =0.0d0
             offset  = numAtomType*(itype-1)+jtype
             forcelj = dr6i*(lj1(offset)*dr6i-lj2(offset))
             potential = potential + dr6i*(lj3(offset)*dr6i-lj4(offset))
             
             
             !***check to see if force correction is correct
             r = sqrt(dr2)
             prefactor =  coulpre*qtmp*q(j)/r
             grij = gewald*r
             expm2 = exp(-grij*grij)
             t = 1.0/(1.0 + EWALD_P*grij)
             erfc = t*(A1+t*(A2+t*(A3+t*(A4+t*A5))))*expm2
             if(dr2.gt.cut_coulsq)prefactor = 0.0d0
             forcecoul = prefactor * (erfc + EWALD_F*grij*expm2)             
             e_coul = e_coul + prefactor*erfc    
             
             
          endif
          force   = (forcecoul+forcelj)*dr2i
          
          ff(i)%x = ff(i)%x + dx*force
          ff(i)%y = ff(i)%y + dy*force
          ff(i)%z = ff(i)%z + dz*force
          
          
          ff(j)%x = ff(j)%x - dx*force
          ff(j)%y = ff(j)%y - dy*force
          ff(j)%z = ff(j)%z - dz*force
          
          


          
       enddo
    enddo
  end subroutine force_N2_14





  

  subroutine RDF_calc
      implicit none
      double precision  :: force,forcelj,forcecoul
      double precision  :: x1,y1,z1,x2,y2,z2
      double precision  :: dx,dy,dz,dr,dr2,dr2i,dr6i,dr12i,dri
      double precision  :: bondscale,bondscalef
      double precision  :: grij,expm2
      double precision  :: ffx,ffy,ffz
      double precision  :: qtmp,r,prefactor,erfc,t
      double precision  :: boxdx,boxdy,boxdz
      integer :: i,j,l,step
      integer :: itype,jtype,neigh
      integer :: tid,num
      integer :: offset,ioffset,neigh_off
      integer :: T1,T2,clock_rate,clock_max
      integer :: numt,threadid
      integer :: node,neighcheck

      !dir$ offload begin target(mic:ourmic),                                       &
      !dir$ in(position: alloc_if(.false.),free_if(.false.)),                       &
      !dir$ inout(rdf: alloc_if(.true.) free_if(.true.)),                           &
      !dir$ nocopy(nlist: alloc_if(.false.) free_if(.false.)),                      &
      !dir$ nocopy(numneigh: alloc_if(.false.) free_if(.false.)),                   &
      !dir$ in(rcut2),in(cut_coulsq),in(box),in(ibox),in(np),                       &
      !dir$ in(numAtomType),in(neigh_alloc),                                        &
      !dir$ inout(numob_rdf)

      numob_rdf = numob_rdf + 1
      do i = 1 ,np
         x1 = position(i)%x; y1 = position(i)%y; z1 = position(i)%z; itype = position(i)%type

         if(itype.eq.3)neighcheck=3
         !if(itype.eq.8)neighcheck=3


         if(itype.eq.3)then
            neigh_off = neigh_alloc*(i-1)
            num = numneigh(i)
            
            do j= 1,num
               
               neigh           = nlist(neigh_off+j)
               if(position(neigh)%type.eq.neighcheck)then
                  
                  dx = x1-position(neigh)%x
                  dy = y1-position(neigh)%y
                  dz = z1-position(neigh)%z
                  jtype = position(neigh)%type
                  
                  boxdx = dx*ibox; boxdy = dy*ibox; boxdz = dz*ibox
                  boxdx = (boxdx+sign(1/(epsilon(boxdx)),boxdx)) -sign(1/epsilon(boxdx),dx)
                  boxdy = (boxdy+sign(1/(epsilon(boxdy)),boxdy)) -sign(1/epsilon(boxdy),dy)
                  boxdz = (boxdz+sign(1/(epsilon(boxdz)),boxdz)) -sign(1/epsilon(boxdz),dz)
                  
                  dx = dx-box*boxdx; dy = dy-box*boxdy; dz = dz-box*boxdz
                  dr2 = dx*dx + dy*dy + dz*dz

                  if(dr2.lt.rcut2)then
                     node = sqrt(dr2)/rdf_delta
                     rdf(node) = rdf(node) + 1
                  endif
               endif
               
            enddo
         endif
      enddo
      !dir$ end offload
    end subroutine RDF_CALC


end module mod_force_nonbond
