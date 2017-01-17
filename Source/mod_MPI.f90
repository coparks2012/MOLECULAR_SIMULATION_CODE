module mod_MPI
  !use ifport
  use global
  use mpi
  
  contains
    
    subroutine MPI_ALLOC()
      implicit none
      
      allocate(p2n(0:size-1),p2w(0:size-1))

    end subroutine MPI_ALLOC
    


    subroutine parallel_temper(w0,n0)
      implicit none
      logical :: swap
      double precision :: E0,E1
      double precision :: delta
      double precision :: rand1
      double precision :: k0temp,k1temp
      integer :: histo1temp(0:nucMax)
      integer :: histo2temp(0:nucMax)
      integer :: len
      integer :: msgtag,ierr
      integer :: status(MPI_STATUS_SIZE)
      integer :: loop
      integer :: w0,n0
      integer :: n0temp,n1temp
      integer :: p0temp,p1temp
      integer :: w0temp,w1temp

      !call mpi_barrier(MPI_COMM_WORLD,ierr)

   
      if(rank.eq.0)then
         
         !---initialize values for rank zero
         
         p2n(0) = n0; p2w(0) = w0; p2kn(0) = kn

         do loop =1,size-1
            call MPI_RECV(p2n(loop),1,MPI_INT,loop,loop,MPI_COMM_WORLD,status,ierr)
            
         enddo
         do loop =1,size-1
            call MPI_RECV(p2w(loop),1,MPI_INT,loop,loop,MPI_COMM_WORLD,status,ierr)
         enddo
         
         do loop = 1,size-1
            call MPI_RECV(replicaEx(:,loop),300,MPI_INT,loop,loop,MPI_COMM_WORLD,status,ierr)
         enddo

      
         do loop = 1,size-1
            call MPI_RECV(p2kn(loop),1,MPI_DOUBLE,loop,loop,MPI_COMM_WORLD,status,ierr)
         enddo



         replicaEx(:,0) = nuchisto(:)
         
         do loop = 0,size-2
            
            n0temp = p2n(loop); n1temp = p2n(loop+1)
            w0temp = p2w(loop); w1temp = p2w(loop+1) 
            k0temp = p2kn(loop); k1temp  = p2kn(loop+1)



            E0 = 0.50d0*k0temp*(w0temp-n0temp)**2 + 0.50d0*k1temp*(w1temp-n1temp)**2 
            E1 = 0.50d0*k1temp*(w0temp-n1temp)**2 + 0.50d0*k0temp*(w1temp-n0temp)**2
            
            delta = (E1-E0)/(kboltz*temp)
            delta = exp(-delta)
            rand1 = rand()


            
            if(rand1.lt.delta)then
               

               histo1temp(:) = replicaEx(:,loop); histo2temp(:) = replicaEx(:,loop+1)
               p2n(loop) = n1temp; p2n(loop+1) = n0temp
               replicaEx(:,loop) = histo2temp(:); replicaEx(:,loop+1) = histo1temp(:)
               p2kn(loop)= k1temp; p2kn(loop+1) = k0temp 
                              

               print*,'succesful swap:',p2n(loop),p2n(loop+1)
                              
            endif
            
         enddo
         
         
         !---now send new values back to processors
         do loop = 1, size -1
            call MPI_SEND(p2n(loop),1,MPI_INT,loop,loop,MPI_COMM_WORLD,ierr)         
         enddo

         do loop = 1, size-1
            call MPI_SEND(replicaEx(:,loop),300,MPI_INT,loop,loop,MPI_COMM_WORLD,ierr)
         enddo

         
         do loop = 1, size-1
            call MPI_SEND(p2kn(loop),1,MPI_DOUBLE,loop,loop,MPI_COMM_WORLD,ierr)
         enddo

         n0 = p2n(0)
         nuchisto(:) = replicaEx(:,0)
         kn = p2kn(0)
         
      endif
      
      if(rank.ne.0)then
         call MPI_SEND(n0,1,MPI_INT,0,rank,MPI_COMM_WORLD,ierr)
         call MPI_SEND(w0,1,MPI_INT,0,rank,MPI_COMM_WORLD,ierr)
         call MPI_SEND(nuchisto(:),300,MPI_INT,0,rank,MPI_COMM_WORLD,ierr)
         call MPI_SEND(kn,1,MPI_DOUBLE,0,rank,MPI_COMM_WORLD,ierr)

      endif
      


      if(rank.ne.0)then
         call MPI_RECV(n0,1,MPI_INT,0,rank,MPI_COMM_WORLD,status,ierr)
         call MPI_RECV(nuchisto(:),300,MPI_INT,0,rank,MPI_COMM_WORLD,status,ierr)
         call MPI_RECV(kn,1,MPI_DOUBLE,0,rank,MPI_COMM_WORLD,status,ierr)
      endif



    end subroutine parallel_temper

  end module mod_MPI
