module mod_matrix_multiply
  implicit none
  
  contains
    
    subroutine matrix_multiply(A,len,rot)
      integer :: i,len
      double precision :: A(3,len), rot(9)
      double precision :: x,y,z
     
      do i=1, len 
         
         x = rot(1)*A(1,i) + rot(2)*A(2,i) + rot(3)*A(3,i)
         y = rot(4)*A(1,i) + rot(5)*A(2,i) + rot(6)*A(3,i)
         z = rot(7)*A(1,i) + rot(8)*A(2,i) + rot(9)*A(3,i)
         
         A(1,i) = x
         A(2,i) = y
         A(3,i) = z
      enddo
    end subroutine matrix_multiply

  end module mod_matrix_multiply
