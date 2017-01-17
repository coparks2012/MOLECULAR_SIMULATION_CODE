program check
  implicit none
  double precision :: sum
  double precision :: x1,y1,z1
  double precision :: x2,y2,z2
  double precision :: dx,dy,dz
  integer :: i,j,k

  open(unit = 1000, file = 'check.dat')

  do j = 1,4
     sum = 0.0d0
     do i = 1,6
        read(1000,*)x1,y1,z1
        read(1000,*)x2,y2,z2
        read(1000,*)
     
        dx = x1-x2; dy = y1-y2; dz = z2-z1
        dx = dx*3.50d0
        dy = dy*3.50d0
        dz = dz*3.50d0
        sum = sum+dx*dx + dy*dy + dz*dz
     end do
     print*,'RMSD:',j,sqrt(sum/10.0d0)
     read(1000,*)
  end do

end program check
