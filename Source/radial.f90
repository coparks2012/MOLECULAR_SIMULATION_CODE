program radial
  implicit none
  integer :: i
  integer :: one,two
  double precision :: three, four, five
  open(unit = 1, file = 'radial.log')
  open(unit = 2, file = 'radial.plot')

  do i = 1 ,99
     read(1,*)one,two,three,four,five
     write(2,*)three,five
  enddo

end program radial
