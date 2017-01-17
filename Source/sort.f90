program sort
  use mod_sort
  implicit none
  double precision :: r(2,5)
  !integer :: b(5)
  integer :: i


  r(1,:) = (/4.1d0,2.1d0,2.05d0,-1.5d0,4.2d0/)
  r(2,:) = (/1.0d0, 3.0d0, 2.0d0, 80.0d0, 6.0d0/)
  !b = rargsort(r,5)
  r(1,:) = r(1,rargsort(r(1,:),5))
 
  do i = 1,5
     print*,r(1,i),r(2,i)
  enddo


end program sort
