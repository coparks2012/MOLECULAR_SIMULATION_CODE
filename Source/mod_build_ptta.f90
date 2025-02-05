module mod_build_ptta
  use global

  implicit none
  contains

    subroutine build_ptta
      implicit none
      integer :: i

      do i = 1,np
         if(mol(1,i).eq.1)then
            ptta( mol(2,i) )%x = position(i)%x
            ptta( mol(2,i) )%y = position(i)%y
            ptta( mol(2,i) )%z = position(i)%z
          

         endif
      enddo

 
    end subroutine build_ptta

       
    subroutine build_ptta_wrap
      implicit none
      double precision :: x1,y1,z1
      double precision :: x2,y2,z2
      double precision :: dx,dy,dz
      integer :: i,j,count,countnp

      count = 0
      countnp = 0
      print*,'building ptta wrap',nummolarray(1),nummolatom(1)
      do i = 1,nummolarray(1)
         do j = 1,nummolatom(1)
            
            countnp = countnp + 1
            if(j.eq.1.and.mol(2,countnp).ne.0)then
               count = count + 1
               x1 = position(countnp)%x
               y1 = position(countnp)%y
               z1 = position(countnp)%z
               
               ptta(count)%x = x1
               ptta(count)%y = y1
               ptta(count)%z = z1

            elseif(mol(2,countnp).ne.0)then
               count = count + 1
               x2 = position(countnp)%x
               y2 = position(countnp)%y
               z2 = position(countnp)%z
            
               dx = x2-x1
               dy = y2-y1
               dz = z2-z1
               
               ptta(count)%x = x2 -box*nint(dx*ibox)
               ptta(count)%y = y2 -box*nint(dy*ibox)
               ptta(count)%z = z2 -box*nint(dz*ibox)
            end if
         enddo
      enddo
    end subroutine build_ptta_wrap


    subroutine build_ptta_wrap_cluster
      implicit none
      double precision :: x2,y2,z2,x1,y1,z1
      double precision :: dx,dy,dz
      integer :: i,j,count,countnp,id
            

      do i = 1 ,np
         id = global_id(i)
         if(mol(1,id).eq.1)then
            ptta(id)%x = position(id)%x
            ptta(id)%y = position(id)%y
            ptta(id)%z = position(id)%z
         endif
      end do

      count = 0
      countnp = 0
      do i = 1,nummolarray(1)
         do j = 1,nummolatom(1)
            
            countnp = countnp + 1 
            if(j.eq.1.and.mol(2,countnp).ne.0)then
               count = count + 1
               x1 = ptta(count)%x
               y1 = ptta(count)%y
               z1 = ptta(count)%z
               

            elseif(mol(2,countnp).ne.0)then
               count = count + 1
               x2 = ptta(count)%x
               y2 = ptta(count)%y
               z2 = ptta(count)%z
            
               dx = x2-x1
               dy = y2-y1
               dz = z2-z1
               
               ptta(count)%x = x2 -box*nint(dx*ibox)
               ptta(count)%y = y2 -box*nint(dy*ibox)
               ptta(count)%z = z2 -box*nint(dz*ibox)
            end if
         enddo
      enddo
    end subroutine build_ptta_wrap_cluster


    subroutine build_ptta_wrap_cluster2
      implicit none
      double precision :: x2,y2,z2,x1,y1,z1
      double precision :: dx,dy,dz
      integer :: i,j,count,countnp,id
      integer :: flag
      
      flag = nummolarray(1)*nummolatom(1)

      do i = 1 ,np
         id = global_id(i)
         if(id.lt.flag)then
            ptta_full(id)%x = position(i)%x
            ptta_full(id)%y = position(i)%y
            ptta_full(id)%z = position(i)%z
         endif
      end do

   
      count = 0
      countnp = 0
      do i = 1,nummolarray(1)
         do j = 1,nummolatom(1)
            
            countnp = countnp + 1 
            if(j.eq.1)then
               count = count + 1
               x1 = ptta_full(count)%x
               y1 = ptta_full(count)%y
               z1 = ptta_full(count)%z
               

            else
               count = count + 1
               x2 = ptta_full(count)%x
               y2 = ptta_full(count)%y
               z2 = ptta_full(count)%z
            
               dx = x2-x1
               dy = y2-y1
               dz = z2-z1
               
               ptta_full(count)%x = x2 -box*nint(dx*ibox)
               ptta_full(count)%y = y2 -box*nint(dy*ibox)
               ptta_full(count)%z = z2 -box*nint(dz*ibox)
            end if
         enddo
      enddo
    end subroutine build_ptta_wrap_cluster2


         
  end module mod_build_ptta
