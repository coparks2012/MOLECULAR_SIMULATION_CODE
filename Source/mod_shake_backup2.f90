module mod_shake
  use global
  use mod_minimg
  
  implicit none
  contains
    
    subroutine shake(iflag)
      implicit none
      integer :: i
      integer :: i0, i1,i2,i3
      integer :: iflag
      integer :: num
      double precision :: bond1,bond2,bond3
      double precision :: dtsq_factor
      integer :: T1,T2,clock_rate,clock_max
      
      call system_clock(T1,clock_rate,clock_max)

      if(iflag.eq.0)dtsq_factor = 0.50d0*dt*dt*kcal2amu
      if(iflag.eq.1)dtsq_factor = dt*dt*kcal2amu

      call shake_update(iflag)


!      !$omp parallel do default(firstprivate),&
!      !$omp& shared(shakeatom,mass,position,ff)
      do i = 1,nshake

         num = shakeatom(i)%num


         if(num.eq.1)then
            bond1 = shakeatom(i)%bond
            i0 = shakeatom(i)%atom1
            i1 = shakeatom(i)%atom2

            call shake2(i0,i1,bond1,dtsq_factor)
         elseif(num.eq.-3)then
            bond1 = shakeatom(i)%bond
            bond2 = shakeatom(i)%bond2
            bond3 = shakeatom(i)%bond3

            i0 = shakeatom(i)%atom1
            i1 = shakeatom(i)%atom2
            i2 = shakeatom(i)%atom3

            call shake3angle(i0,i1,i2,bond1,bond2,bond3,dtsq_factor)


         elseif(num.eq.3)then
            bond1 = shakeatom(i)%bond
            bond2 = shakeatom(i)%bond2

            i0 = shakeatom(i)%atom1
            i1 = shakeatom(i)%atom2
            i2 = shakeatom(i)%atom3

            call shake3(i0,i1,i2,bond1,bond2,dtsq_factor)

            

         elseif(num.eq.4)then
            bond1 = shakeatom(i)%bond
            bond2 = shakeatom(i)%bond2
            bond3 = shakeatom(i)%bond3

            i0 = shakeatom(i)%atom1
            i1 = shakeatom(i)%atom2
            i2 = shakeatom(i)%atom3
            i3 = shakeatom(i)%atom4

            call shake4(i0,i1,i2,i3,bond1,bond2,bond3,dtsq_factor)
            
         endif
      enddo
 !     !$omp end parallel do
      
      call system_clock(T2,clock_rate,clock_max)

      shake_time = shake_time + real(T2-T1)/real(clock_rate)



    end subroutine shake



    subroutine shake2(i0,i1,bond1,dtsq_factor)
      implicit none

      ! argument variables

      integer i0,i1,ifactor,vflag
      double precision :: bond1

      ! local variables
      double precision :: ff0x,ff0y,ff0z
      double precision :: ff1x,ff1y,ff1z
      double precision :: x0,y0,z0
      double precision :: x1,y1,z1
      double precision :: r01sq,s01sq
      double precision :: a,b,c,determ
      double precision :: lamda,lamda1,lamda2
      double precision :: r01(3),s01(3)
      double precision :: invmass0,invmass1
      double precision :: dtsq_factor

      invmass0 = 1.0/mass(position(i0)%type)
      invmass1 = 1.0/mass(position(i1)%type)

      ff0x = ff(i0)%x; ff0y = ff(i0)%y; ff0z = ff(i0)%z
      ff1x = ff(i1)%x; ff1y = ff(i1)%y; ff1z = ff(i1)%z

      ! r01 = distance vec between atoms, with PBC
      r01(1) = position(i0)%x - position(i1)%x
      r01(2) = position(i0)%y - position(i1)%y
      r01(3) = position(i0)%z - position(i1)%z
      call minimg(r01(1),r01(2),r01(3))



      ! s01 = distance vec after unconstrained update, with PBC   
      s01(1) = xshake(i0)%x-xshake(i1)%x
      s01(2) = xshake(i0)%y-xshake(i1)%y
      s01(3) = xshake(i0)%z-xshake(i1)%z
      call minimg(s01(1),s01(2),s01(3))

      ! scalar distances between atoms

      r01sq = r01(1)*r01(1) + r01(2)*r01(2) + r01(3)*r01(3)
      s01sq = s01(1)*s01(1) + s01(2)*s01(2) + s01(3)*s01(3)

      ! a,b,c = coeffs in quadratic equation for lamda

      a = (invmass0+invmass1)*(invmass0+invmass1) * r01sq
      b = 2.0 * (invmass0+invmass1) * &
         (s01(1)*r01(1) + s01(2)*r01(2) + s01(3)*r01(3))
      c = s01sq - bond1*bond1

      ! error check

      determ = b*b - 4.0*a*c
      if (determ < 0.0) then
        print*, 'WARNING: SHAKE determinant < 0.0',i0,i1
        determ = 0.0
        print*,b,a,c
        stop
      endif

      ! exact quadratic solution for lamda

      lamda1 = (-b+sqrt(determ)) / (2.0*a)
      lamda2 = (-b-sqrt(determ)) / (2.0*a)

      if (abs(lamda1) <= abs(lamda2)) then
        lamda = lamda1
      else
        lamda = lamda2
      endif

      ! update forces if atom is owned by this processor

      lamda = lamda/dtsq_factor

 

      ff0x = ff0x + lamda*r01(1)
      ff0y = ff0y + lamda*r01(2)
      ff0z = ff0z + lamda*r01(3)
      
      ff1x = ff1x - lamda*r01(1)
      ff1y = ff1y - lamda*r01(2)
      ff1z = ff1z - lamda*r01(3)

      ff(i0)%x = ff0x
      ff(i0)%y = ff0y
      ff(i0)%z = ff0z

      
      ff(i1)%x = ff1x
      ff(i1)%y = ff1y
      ff(i1)%z = ff1z


      virialx = virialx + r01(1)*r01(1)*lamda
      virialy = virialy + r01(2)*r01(2)*lamda
      virialz = virialz + r01(3)*r01(3)*lamda



      return
    end subroutine shake2



    ! -------------------------------------------------------------------------
    ! 3-atom SHAKE

    subroutine shake3(i0,i1,i2,bond1,bond2,dtsq_factor)
      implicit none

      ! argument variables
      integer i0,i1,i2,ifactor,vflag
      real*8 bond1,bond2,dtsq_factor

      ! local variables
      double precision :: ff0x,ff0y,ff0z
      double precision :: ff1x,ff1y,ff1z
      double precision :: ff2x,ff2y,ff2z
      real*8 s01sq,s02sq
      real*8 a11,a12,a21,a22,b1,b2
      real*8 a11inv,a12inv,a21inv,a22inv
      real*8 determ,determinv
      real*8 lamda01,lamda02
      real*8 r01(3),r02(3),s01(3),s02(3)
      real*8 invmass0,invmass1,invmass2
      real*8 quad1,quad2
      real*8 quad1_0101,quad1_0202,quad1_0102
      real*8 quad2_0101,quad2_0202,quad2_0102
      real*8 lamda01_new,lamda02_new
      real*8 r01sq,r02sq,r0102
      integer niter
      logical done

      invmass0 = 1.0/mass(position(i0)%type)
      invmass1 = 1.0/mass(position(i1)%type)
      invmass2 = 1.0/mass(position(i2)%type)
      ff0x = ff(i0)%x; ff0y = ff(i0)%y; ff0z = ff(i0)%z
      ff1x = ff(i1)%x; ff1y = ff(i1)%y; ff1z = ff(i1)%z
      ff2x = ff(i2)%x; ff2y = ff(i2)%y; ff2z = ff(i2)%z

      ! r01,r02 = distance vecs between atoms, with PBC

      r01(1) = position(i0)%x - position(i1)%x
      r01(2) = position(i0)%y - position(i1)%y
      r01(3) = position(i0)%z - position(i1)%z
      call minimg(r01(1),r01(2),r01(3))

      r02(1) = position(i0)%x - position(i2)%x
      r02(2) = position(i0)%y - position(i2)%y
      r02(3) = position(i0)%z - position(i2)%z
      call minimg(r02(1),r02(2),r02(3))

      ! s01,s02 = distance vecs after unconstrained update, with PBC

      s01(1) = xshake(i0)%x - xshake(i1)%x
      s01(2) = xshake(i0)%y - xshake(i1)%y
      s01(3) = xshake(i0)%z - xshake(i1)%z
      call minimg(s01(1),s01(2),s01(3))

      s02(1) = xshake(i0)%x - xshake(i2)%x
      s02(2) = xshake(i0)%y - xshake(i2)%y
      s02(3) = xshake(i0)%z - xshake(i2)%z
      call minimg(s02(1),s02(2),s02(3))

      ! scalar distances between atoms

      s01sq = s01(1)*s01(1) + s01(2)*s01(2) + s01(3)*s01(3)
      s02sq = s02(1)*s02(1) + s02(2)*s02(2) + s02(3)*s02(3)

      ! matrix coeffs and rhs for lamda equations

      a11 = 2.0 * (invmass0+invmass1) *&
          (s01(1)*r01(1) + s01(2)*r01(2) + s01(3)*r01(3))
      a12 = 2.0 * invmass0 *&
          (s01(1)*r02(1) + s01(2)*r02(2) + s01(3)*r02(3))
      a21 = 2.0 * invmass0 *&
          (s02(1)*r01(1) + s02(2)*r01(2) + s02(3)*r01(3))
      a22 = 2.0 * (invmass0+invmass2) *&
          (s02(1)*r02(1) + s02(2)*r02(2) + s02(3)*r02(3))

      b1 = bond1*bond1 - s01sq
      b2 = bond2*bond2 - s02sq

      ! inverse of matrix

      determ = a11*a22 - a12*a21
      if (determ == 0.0) call abort('SHAKE determinant = 0.0')
      determinv = 1.0/determ

      a11inv = a22*determinv
      a12inv = -a12*determinv
      a21inv = -a21*determinv
      a22inv = a11*determinv

      !quadratic correction coeffs

      r01sq = r01(1)*r01(1) + r01(2)*r01(2) + r01(3)*r01(3)
      r02sq = r02(1)*r02(1) + r02(2)*r02(2) + r02(3)*r02(3)
      r0102 = (r01(1)*r02(1) + r01(2)*r02(2) + r01(3)*r02(3))

      quad1_0101 = (invmass0+invmass1)**2 * r01sq
      quad1_0202 = invmass0*invmass0 * r02sq
      quad1_0102 = 2.0 * (invmass0+invmass1)*invmass0 * r0102

      quad2_0202 = (invmass0+invmass2)**2 * r02sq
      quad2_0101 = invmass0*invmass0 * r01sq
      quad2_0102 = 2.0 * (invmass0+invmass2)*invmass0 * r0102

      ! iterate until converged

      lamda01 = 0.0
      lamda02 = 0.0
      niter = 0
      done = .FALSE.

      do while (.not. done .and. niter < maxiter)

        quad1 = quad1_0101 * lamda01*lamda01 + &
            quad1_0202 * lamda02*lamda02 + &
            quad1_0102 * lamda01*lamda02

        quad2 = quad2_0101 * lamda01*lamda01 + &
            quad2_0202 * lamda02*lamda02 + &
            quad2_0102 * lamda01*lamda02
        
        b1 = bond1*bond1 - s01sq - quad1
        b2 = bond2*bond2 - s02sq - quad2
        
        lamda01_new = a11inv*b1 + a12inv*b2
        lamda02_new = a21inv*b1 + a22inv*b2

        done = .TRUE.
        if (abs(lamda01_new-lamda01) .gt. tolerance) done = .FALSE.
        if (abs(lamda02_new-lamda02) .gt. tolerance) done = .FALSE.

        lamda01 = lamda01_new
        lamda02 = lamda02_new
        niter = niter + 1

      enddo

      if(niter.eq.maxiter)then
         print*,'warning:,niter is maxiter in shake_3'
         print*,'warning generated from rank',rank
         print*
      endif

      ! update forces if atom is owned by this processor

      lamda01 = lamda01/dtsq_factor
      lamda02 = lamda02/dtsq_factor
 
      !---- update forces in local accumulators
      ff0x = ff0x + lamda01*r01(1) + lamda02*r02(1)
      ff0y = ff0y + lamda01*r01(2) + lamda02*r02(2)
      ff0z = ff0z + lamda01*r01(3) + lamda02*r02(3)
      
      ff1x = ff1x - lamda01*r01(1)
      ff1y = ff1y - lamda01*r01(2)
      ff1z = ff1z - lamda01*r01(3)
      
      ff2x = ff2x - lamda02*r02(1)
      ff2y = ff2y - lamda02*r02(2)
      ff2z = ff2z - lamda02*r02(3)
      

      ff(i0)%x = ff0x
      ff(i0)%y = ff0y
      ff(i0)%z = ff0z
      
      ff(i1)%x = ff1x
      ff(i1)%y = ff1y
      ff(i1)%z = ff1z

      ff(i2)%x = ff2x
      ff(i2)%y = ff2y
      ff(i2)%z = ff2z
      

      virialx = virialx + (r01(1)*r01(1)*lamda01)
      virialy = virialy + (r01(2)*r01(2)*lamda01)
      virialz = virialz + (r01(3)*r01(3)*lamda01)
    
      virialx = virialx + (r02(1)*r02(1)*lamda02)
      virialy = virialy + (r02(2)*r02(2)*lamda02)
      virialz = virialz + (r02(3)*r02(3)*lamda02)
   

      
    end subroutine shake3


    ! -------------------------------------------------------------------------
    ! 4-atom SHAKE

    subroutine shake4(i0,i1,i2,i3,bond1,bond2,bond3,dtsq_factor)
      
      implicit none

      ! argument variables

      integer i0,i1,i2,i3
      real*8 bond1,bond2,bond3,dtsq_factor

      ! local variables
      double precision :: ff0x,ff0y,ff0z
      double precision :: ff1x,ff1y,ff1z
      double precision :: ff2x,ff2y,ff2z
      double precision :: ff3x,ff3y,ff3z
      real*8 s01sq,s02sq,s03sq
      real*8 a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3
      real*8 a11inv,a12inv,a13inv,a21inv,a22inv,a23inv
      real*8 a31inv,a32inv,a33inv
      real*8 ainv,binv,cinv,dinv,einv,finv,ginv,hinv,iinv
      real*8 determ,determinv
      real*8 lamda01,lamda02,lamda03
      real*8 r01(3),r02(3),r03(3),s01(3),s02(3),s03(3)
      real*8 invmass0,invmass1,invmass2,invmass3
      real*8 quad1,quad2,quad3
      real*8 quad1_0101,quad1_0202,quad1_0303
      real*8 quad1_0102,quad1_0103,quad1_0203
      real*8 quad2_0101,quad2_0202,quad2_0303
      real*8 quad2_0102,quad2_0103,quad2_0203
      real*8 quad3_0101,quad3_0202,quad3_0303
      real*8 quad3_0102,quad3_0103,quad3_0203
      real*8 lamda01_new,lamda02_new,lamda03_new
      real*8 r01sq,r02sq,r03sq,r0102,r0103,r0203
      integer niter
      logical done



      invmass0 = 1.0/mass(position(i0)%type)
      invmass1 = 1.0/mass(position(i1)%type)
      invmass2 = 1.0/mass(position(i2)%type)
      invmass3 = 1.0/mass(position(i3)%type)


      ff0x = ff(i0)%x; ff0y = ff(i0)%y; ff0z = ff(i0)%z
      ff1x = ff(i1)%x; ff1y = ff(i1)%y; ff1z = ff(i1)%z
      ff2x = ff(i2)%x; ff2y = ff(i2)%y; ff2z = ff(i2)%z
      ff3x = ff(i3)%x; ff3y = ff(i3)%y; ff3z = ff(i3)%z


      ! r01,r02,r03 = distance vecs between atoms, with PBC

      r01(1) = position(i0)%x - position(i1)%x
      r01(2) = position(i0)%y - position(i1)%y
      r01(3) = position(i0)%z - position(i1)%z
      call minimg(r01(1),r01(2),r01(3))

      r02(1) = position(i0)%x - position(i2)%x
      r02(2) = position(i0)%y - position(i2)%y
      r02(3) = position(i0)%z - position(i2)%z
      call minimg(r02(1),r02(2),r02(3))

      r03(1) = position(i0)%x - position(i3)%x
      r03(2) = position(i0)%y - position(i3)%y
      r03(3) = position(i0)%z - position(i3)%z
      call minimg(r03(1),r03(2),r03(3))

      ! s01,s02,s03 = distance vecs after unconstrained update, with PBC

      s01(1) = xshake(i0)%x - xshake(i1)%x
      s01(2) = xshake(i0)%y - xshake(i1)%y
      s01(3) = xshake(i0)%z - xshake(i1)%z
      call minimg(s01(1),s01(2),s01(3))

      s02(1) = xshake(i0)%x - xshake(i2)%x
      s02(2) = xshake(i0)%y - xshake(i2)%y
      s02(3) = xshake(i0)%z - xshake(i2)%z
      call minimg(s02(1),s02(2),s02(3))

      s03(1) = xshake(i0)%x - xshake(i3)%x
      s03(2) = xshake(i0)%y - xshake(i3)%y
      s03(3) = xshake(i0)%z - xshake(i3)%z
      call minimg(s03(1),s03(2),s03(3))

      ! scalar distances between atoms

      r01sq = r01(1)*r01(1) + r01(2)*r01(2) + r01(3)*r01(3)
      r02sq = r02(1)*r02(1) + r02(2)*r02(2) + r02(3)*r02(3)
      r03sq = r03(1)*r03(1) + r03(2)*r03(2) + r03(3)*r03(3)
      s01sq = s01(1)*s01(1) + s01(2)*s01(2) + s01(3)*s01(3)
      s02sq = s02(1)*s02(1) + s02(2)*s02(2) + s02(3)*s02(3)
      s03sq = s03(1)*s03(1) + s03(2)*s03(2) + s03(3)*s03(3)

      ! matrix coeffs and rhs for lamda equations

      a11 = 2.0 * (invmass0+invmass1) *&
          (s01(1)*r01(1) + s01(2)*r01(2) + s01(3)*r01(3))
      a12 = 2.0 * invmass0 *&
          (s01(1)*r02(1) + s01(2)*r02(2) + s01(3)*r02(3))
      a13 = 2.0 * invmass0 *&
          (s01(1)*r03(1) + s01(2)*r03(2) + s01(3)*r03(3))
      a21 = 2.0 * invmass0 *&
          (s02(1)*r01(1) + s02(2)*r01(2) + s02(3)*r01(3))
      a22 = 2.0 * (invmass0+invmass2) *&
          (s02(1)*r02(1) + s02(2)*r02(2) + s02(3)*r02(3))
      a23 = 2.0 * invmass0 *&
          (s02(1)*r03(1) + s02(2)*r03(2) + s02(3)*r03(3))
      a31 = 2.0 * invmass0 *&
          (s03(1)*r01(1) + s03(2)*r01(2) + s03(3)*r01(3))
      a32 = 2.0 * invmass0 *&
          (s03(1)*r02(1) + s03(2)*r02(2) + s03(3)*r02(3))
      a33 = 2.0 * (invmass0+invmass3) *&
          (s03(1)*r03(1) + s03(2)*r03(2) + s03(3)*r03(3))

      b1 = bond1*bond1 - s01sq
      b2 = bond2*bond2 - s02sq
      b3 = bond3*bond3 - s03sq

      !` inverse of matrix

      determ = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 -&
          a11*a23*a32 - a12*a21*a33 - a13*a22*a31
      if (determ == 0.0) call abort('SHAKE determinant = 0.0')
      determinv = 1.0/determ

      a11inv = determinv * (a22*a33 - a23*a32)
      a12inv = -determinv * (a12*a33 - a13*a32)
      a13inv = determinv * (a12*a23 - a13*a22)
      a21inv = -determinv * (a21*a33 - a23*a31)
      a22inv = determinv * (a11*a33 - a13*a31)
      a23inv = -determinv * (a11*a23 - a13*a21)
      a31inv = determinv * (a21*a32 - a22*a31)
      a32inv = -determinv * (a11*a32 - a12*a31)
      a33inv = determinv * (a11*a22 - a12*a21)

      !quadratic correction coeffs

      r01sq = r01(1)*r01(1) + r01(2)*r01(2) + r01(3)*r01(3)
      r02sq = r02(1)*r02(1) + r02(2)*r02(2) + r02(3)*r02(3)
      r03sq = r03(1)*r03(1) + r03(2)*r03(2) + r03(3)*r03(3)
      r0102 = (r01(1)*r02(1) + r01(2)*r02(2) + r01(3)*r02(3))
      r0103 = (r01(1)*r03(1) + r01(2)*r03(2) + r01(3)*r03(3))
      r0203 = (r02(1)*r03(1) + r02(2)*r03(2) + r02(3)*r03(3))

      quad1_0101 = (invmass0+invmass1)**2 * r01sq
      quad1_0202 = invmass0*invmass0 * r02sq
      quad1_0303 = invmass0*invmass0 * r03sq
      quad1_0102 = 2.0 * (invmass0+invmass1)*invmass0 * r0102
      quad1_0103 = 2.0 * (invmass0+invmass1)*invmass0 * r0103
      quad1_0203 = 2.0 * invmass0*invmass0 * r0203

      quad2_0101 = invmass0*invmass0 * r01sq
      quad2_0202 = (invmass0+invmass2)**2 * r02sq
      quad2_0303 = invmass0*invmass0 * r03sq
      quad2_0102 = 2.0 * (invmass0+invmass2)*invmass0 * r0102
      quad2_0103 = 2.0 * invmass0*invmass0 * r0103
      quad2_0203 = 2.0 * (invmass0+invmass2)*invmass0 * r0203

      quad3_0101 = invmass0*invmass0 * r01sq
      quad3_0202 = invmass0*invmass0 * r02sq
      quad3_0303 = (invmass0+invmass3)**2 * r03sq
      quad3_0102 = 2.0 * invmass0*invmass0 * r0102
      quad3_0103 = 2.0 * (invmass0+invmass3)*invmass0 * r0103
      quad3_0203 = 2.0 * (invmass0+invmass3)*invmass0 * r0203

      ! iterate until converged

      lamda01 = 0.0
      lamda02 = 0.0
      lamda03 = 0.0
      niter = 0
      done = .FALSE.

      do while (.not. done .and. niter < maxiter)

        quad1 = quad1_0101 * lamda01*lamda01 + &
            quad1_0202 * lamda02*lamda02 +&
            quad1_0303 * lamda03*lamda03 + &
            quad1_0102 * lamda01*lamda02 +&
            quad1_0103 * lamda01*lamda03 +&
            quad1_0203 * lamda02*lamda03

        quad2 = quad2_0101 * lamda01*lamda01 +& 
            quad2_0202 * lamda02*lamda02 +&
            quad2_0303 * lamda03*lamda03 + &
            quad2_0102 * lamda01*lamda02 +&
            quad2_0103 * lamda01*lamda03 +&
            quad2_0203 * lamda02*lamda03

        quad3 = quad3_0101 * lamda01*lamda01 +& 
            quad3_0202 * lamda02*lamda02 +&
            quad3_0303 * lamda03*lamda03 +&
            quad3_0102 * lamda01*lamda02 +&
            quad3_0103 * lamda01*lamda03 +&
            quad3_0203 * lamda02*lamda03

        b1 = bond1*bond1 - s01sq - quad1
        b2 = bond2*bond2 - s02sq - quad2
        b3 = bond3*bond3 - s03sq - quad3
        
        lamda01_new = a11inv*b1 + a12inv*b2 + a13inv*b3
        lamda02_new = a21inv*b1 + a22inv*b2 + a23inv*b3
        lamda03_new = a31inv*b1 + a32inv*b2 + a33inv*b3
        
        done = .TRUE.
        if (abs(lamda01_new-lamda01) .gt. tolerance) done = .FALSE.
        if (abs(lamda02_new-lamda02) .gt. tolerance) done = .FALSE.
        if (abs(lamda03_new-lamda03) .gt. tolerance) done = .FALSE.



        lamda01 = lamda01_new
        lamda02 = lamda02_new
        lamda03 = lamda03_new
        niter = niter + 1

      enddo


      if(niter.eq.maxiter)then
         print*,'warning:,niter is maxiter in shake_4',niter,maxiter
         print*,'what the hell is tolerance',tolerance
         print*,'warning generated from rank',rank
         print*,abs(lamda01_new-lamda01)
         print*,abs(lamda02_new-lamda02)
         print*,abs(lamda03_new-lamda03)
         print*
         
      endif
      
      ! update forces if atom is owned by this processor

      lamda01 = lamda01/dtsq_factor
      lamda02 = lamda02/dtsq_factor
      lamda03 = lamda03/dtsq_factor

      ff0x = ff0x +&
           lamda01*r01(1) + lamda02*r02(1) + lamda03*r03(1)
      ff0y = ff0y +&
           lamda01*r01(2) + lamda02*r02(2) + lamda03*r03(2)
      ff0z = ff0z +&
           lamda01*r01(3) + lamda02*r02(3) + lamda03*r03(3)
      
      ff1x = ff1x - lamda01*r01(1)
      ff1y = ff1y - lamda01*r01(2)
      ff1z = ff1z - lamda01*r01(3)
      
      ff2x = ff2x - lamda02*r02(1)
      ff2y = ff2y - lamda02*r02(2)
      ff2z = ff2z - lamda02*r02(3)
      
      ff3x = ff3x - lamda03*r03(1)
      ff3y = ff3y - lamda03*r03(2)
      ff3z = ff3z - lamda03*r03(3)


      ff(i0)%x = ff0x
      ff(i0)%y = ff0y
      ff(i0)%z = ff0z
     
      ff(i1)%x = ff1x
      ff(i1)%y = ff1y
      ff(i1)%z = ff1z
      
      ff(i2)%x = ff2x
      ff(i2)%y = ff2y
      ff(i2)%z = ff2z
      
      ff(i3)%x = ff3x
      ff(i3)%y = ff3y
      ff(i3)%z = ff3z
      
      
      virialx = virialx + r01(1)*r01(1)*lamda01
      virialy = virialy + r01(2)*r01(2)*lamda01
      virialz = virialz + r01(3)*r01(3)*lamda01
      
      virialx = virialx + r02(1)*r02(1)*lamda02
      virialy = virialy + r02(2)*r02(2)*lamda02
      virialz = virialz + r02(3)*r02(3)*lamda02
      
      virialx = virialx + r03(1)*r03(1)*lamda03
      virialy = virialy + r03(2)*r03(2)*lamda03
      virialz = virialz + r03(3)*r03(3)*lamda03
      
      

    end subroutine shake4





    ! -------------------------------------------------------------------------
    ! 3-atom SHAKE with fixed angle

    subroutine shake3angle(i0,i1,i2,bond1,bond2,bond12,dtsq_factor)
      
      implicit none
     
      !argument variables
      
      integer i0,i1,i2,ifactor,vflag
      real*8 bond1,bond2,bond12,dtsq_factor
      
      ! local variables
      double precision :: ff0x,ff0y,ff0z
      double precision :: ff1x,ff1y,ff1z
      double precision :: ff2x,ff2y,ff2z
      real*8 dtsq,s01sq,s02sq,s12sq
      real*8 a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3
      real*8 a11inv,a12inv,a13inv,a21inv,a22inv,a23inv
      real*8 a31inv,a32inv,a33inv
      real*8 ainv,binv,cinv,dinv,einv,finv,ginv,hinv,iinv
      real*8 determ,determinv
      real*8 lamda01,lamda02,lamda12
      real*8 s0(3),s1(3),s2(3),s3(3)
      real*8 r01(3),r02(3),r12(3),s01(3),s02(3),s12(3)
      real*8 invmass0,invmass1,invmass2
      real*8 quad1,quad2,quad3
      real*8 quad1_0101,quad1_0202,quad1_1212
      real*8 quad1_0102,quad1_0112,quad1_0212
      real*8 quad2_0101,quad2_0202,quad2_1212
      real*8 quad2_0102,quad2_0112,quad2_0212
      real*8 quad3_0101,quad3_0202,quad3_1212
      real*8 quad3_0102,quad3_0112,quad3_0212
      real*8 lamda01_new,lamda02_new,lamda12_new
      real*8 r01sq,r02sq,r12sq,r0102,r0112,r0212
      double precision :: xi0,yi0,zi0
      double precision :: xi1,yi1,z11
      double precision :: x12,y12,z12
      integer niter
      logical done
    
      invmass0 = 1.0/mass(position(i0)%type)
      invmass1 = 1.0/mass(position(i1)%type)
      invmass2 = 1.0/mass(position(i2)%type)

      ff0x = ff(i0)%x; ff0y = ff(i0)%y; ff0z = ff(i0)%z
      ff1x = ff(i1)%x; ff1y = ff(i1)%y; ff1z = ff(i1)%z
      ff2x = ff(i2)%x; ff2y = ff(i2)%y; ff2z = ff(i2)%z

      
      ! r01,r02,r12 = distance vecs between atoms, with PBC
      
      r01(1) = position(i0)%x - position(i1)%x
      r01(2) = position(i0)%y - position(i1)%y
      r01(3) = position(i0)%z - position(i1)%z
      call minimg(r01(1),r01(2),r01(3))
      
      r02(1) = position(i0)%x - position(i2)%x
      r02(2) = position(i0)%y - position(i2)%y
      r02(3) = position(i0)%z - position(i2)%z
      call minimg(r02(1),r02(2),r02(3))
      
      r12(1) = position(i1)%x - position(i2)%x
      r12(2) = position(i1)%y - position(i2)%y
      r12(3) = position(i1)%z - position(i2)%z
      call minimg(r12(1),r12(2),r12(3))
      
      ! s01,s02,s12 = distance vecs after unconstrained update, with PBC
     
      s01(1) = xshake(i0)%x - xshake(i1)%x
      s01(2) = xshake(i0)%y - xshake(i1)%y
      s01(3) = xshake(i0)%z - xshake(i1)%z
      call minimg(s01(1),s01(2),s01(3))
     
      s02(1) = xshake(i0)%x - xshake(i2)%x
      s02(2) = xshake(i0)%y - xshake(i2)%y
      s02(3) = xshake(i0)%z - xshake(i2)%z
      call minimg(s02(1),s02(2),s02(3))
     
      s12(1) = xshake(i1)%x - xshake(i2)%x
      s12(2) = xshake(i1)%y - xshake(i2)%y
      s12(3) = xshake(i1)%z - xshake(i2)%z
      call minimg(s12(1),s12(2),s12(3))
     
     !scalar distances between atoms
     
      r01sq = r01(1)*r01(1) + r01(2)*r01(2) + r01(3)*r01(3)
      r02sq = r02(1)*r02(1) + r02(2)*r02(2) + r02(3)*r02(3)
      r12sq = r12(1)*r12(1) + r12(2)*r12(2) + r12(3)*r12(3)
      s01sq = s01(1)*s01(1) + s01(2)*s01(2) + s01(3)*s01(3)
      s02sq = s02(1)*s02(1) + s02(2)*s02(2) + s02(3)*s02(3)
      s12sq = s12(1)*s12(1) + s12(2)*s12(2) + s12(3)*s12(3)
     
     ! matrix coeffs and rhs for lamda equations
     
      a11 = 2.0 * (invmass0+invmass1) *&
           (s01(1)*r01(1) + s01(2)*r01(2) + s01(3)*r01(3))
      a12 = 2.0 * invmass0 *&
           (s01(1)*r02(1) + s01(2)*r02(2) + s01(3)*r02(3))
      a13 = - 2.0 * invmass1 *&
           (s01(1)*r12(1) + s01(2)*r12(2) + s01(3)*r12(3))
      a21 = 2.0 * invmass0 *&
           (s02(1)*r01(1) + s02(2)*r01(2) + s02(3)*r01(3))
      a22 = 2.0 * (invmass0+invmass2) *&
           (s02(1)*r02(1) + s02(2)*r02(2) + s02(3)*r02(3))
      a23 = 2.0 * invmass2 *&
           (s02(1)*r12(1) + s02(2)*r12(2) + s02(3)*r12(3))
      a31 = - 2.0 * invmass1 *&
           (s12(1)*r01(1) + s12(2)*r01(2) + s12(3)*r01(3))
      a32 = 2.0 * invmass2 *&
           (s12(1)*r02(1) + s12(2)*r02(2) + s12(3)*r02(3))
      a33 = 2.0 * (invmass1+invmass2) *&
           (s12(1)*r12(1) + s12(2)*r12(2) + s12(3)*r12(3))
     

      determ = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 -&
           a11*a23*a32 - a12*a21*a33 - a13*a22*a31
      if (determ == 0.0) call abort('SHAKE determinant = 0.0')
      determinv = 1.0d0/determ


      a11inv = determinv * (a22*a33 - a23*a32)
      a12inv = -determinv * (a12*a33 - a13*a32)
      a13inv = determinv * (a12*a23 - a13*a22)
      a21inv = -determinv * (a21*a33 - a23*a31)
      a22inv = determinv * (a11*a33 - a13*a31)
      a23inv = -determinv * (a11*a23 - a13*a21)
      a31inv = determinv * (a21*a32 - a22*a31)
      a32inv = -determinv * (a11*a32 - a12*a31)
      a33inv = determinv * (a11*a22 - a12*a21)

      !----quadratic correction coeffs

      quad1_0101 = (invmass0+invmass1)*(invmass0+invmass1) * r01sq
      quad1_0202 = invmass0*invmass0 * r02sq
      quad1_1212 = invmass1*invmass1 * r12sq
      quad1_0102 = 2.0 * (invmass0+invmass1)*invmass0 * r0102
      quad1_0112 = - 2.0 * (invmass0+invmass1)*invmass1 * r0112
      quad1_0212 = - 2.0 * invmass0*invmass1 * r0212
     
      quad2_0101 = invmass0*invmass0 * r01sq
      quad2_0202 = (invmass0+invmass2)**2 * r02sq
      quad2_1212 = invmass2*invmass2 * r12sq
      quad2_0102 = 2.0 * (invmass0+invmass2)*invmass0 * r0102
      quad2_0112 = 2.0 * invmass0*invmass2 * r0112
      quad2_0212 = 2.0 * (invmass0+invmass2)*invmass2 * r0212
     
      quad3_0101 = invmass1*invmass1 * r01sq
      quad3_0202 = invmass2*invmass2 * r02sq
      quad3_1212 = (invmass1+invmass2)**2 * r12sq
      quad3_0102 = - 2.0 * invmass1*invmass2 * r0102
      quad3_0112 = - 2.0 * (invmass1+invmass2)*invmass1 * r0112
      quad3_0212 = 2.0 * (invmass1+invmass2)*invmass2 * r0212
     
       
     
     ! quadratic correction coeffs
     ! r01sq = r01(1)*r01(1) + r01(2)*r01(2) + r01(3)*r01(3)
     ! r02sq = r02(1)*r02(1) + r02(2)*r02(2) + r02(3)*r02(3)
     ! r12sq = r12(1)*r12(1) + r12(2)*r12(2) + r12(3)*r12(3)
     ! r0102 = (r01(1)*r02(1) + r01(2)*r02(2) + r01(3)*r02(3))
     ! r0112 = (r01(1)*r12(1) + r01(2)*r12(2) + r01(3)*r12(3))
     ! r0212 = (r02(1)*r12(1) + r02(2)*r12(2) + r02(3)*r12(3))
     
    
     ! iterate until converged
      lamda01 = 0.0d0
      lamda02 = 0.0d0
      lamda12 = 0.0d0
      niter = 0
      done = .FALSE.
     
      do while (.not. done .and. niter < maxiter)
        
         quad1 = quad1_0101 * lamda01*lamda01 + &
              quad1_0202 * lamda02*lamda02 + &
              quad1_1212 * lamda12*lamda12 + &
              quad1_0102 * lamda01*lamda02 + &
              quad1_0112 * lamda01*lamda12 + &
              quad1_0212 * lamda02*lamda12
        
         quad2 = quad2_0101 * lamda01*lamda01 + & 
              quad2_0202 * lamda02*lamda02 + &
              quad2_1212 * lamda12*lamda12 + &
              quad2_0102 * lamda01*lamda02 + &
              quad2_0112 * lamda01*lamda12 + &
              quad2_0212 * lamda02*lamda12
         
         quad3 = quad3_0101 * lamda01*lamda01 + &
              quad3_0202 * lamda02*lamda02 + &
              quad3_1212 * lamda12*lamda12 + & 
              quad3_0102 * lamda01*lamda02 + &
              quad3_0112 * lamda01*lamda12 + &
              quad3_0212 * lamda02*lamda12
        
         b1 = bond1*bond1 - s01sq - quad1
         b2 = bond2*bond2 - s02sq - quad2
         b3 = bond12*bond12 - s12sq - quad3
        
         lamda01_new = a11inv*b1 + a12inv*b2 + a13inv*b3
         lamda02_new = a21inv*b1 + a22inv*b2 + a23inv*b3
         lamda12_new = a31inv*b1 + a32inv*b2 + a33inv*b3
        
         done = .TRUE.
         if (abs(lamda01_new-lamda01) .gt. tolerance) done = .FALSE.
         if (abs(lamda02_new-lamda02) .gt. tolerance) done = .FALSE.
         if (abs(lamda12_new-lamda12) .gt. tolerance) done = .FALSE.
         
         lamda01 = lamda01_new
         lamda02 = lamda02_new
         lamda12 = lamda12_new
         niter = niter + 1
        
      enddo

       if(niter.eq.maxiter)then
         print*,'warning:,niter is maxiter in shake_3_angle'
      endif

     !update forces if atom is owned by this processor
     
      lamda01 = lamda01/dtsq_factor
      lamda02 = lamda02/dtsq_factor
      lamda12 = lamda12/dtsq_factor
     
      ff0x = ff0x + lamda01*r01(1) + lamda02*r02(1)
      ff0y = ff0y + lamda01*r01(2) + lamda02*r02(2)
      ff0z = ff0z + lamda01*r01(3) + lamda02*r02(3)
     
      ff1x = ff1x - lamda01*r01(1) + lamda12*r12(1)
      ff1y = ff1y - lamda01*r01(2) + lamda12*r12(2)
      ff1z = ff1z - lamda01*r01(3) + lamda12*r12(3)
     
      ff2x = ff2x - lamda02*r02(1) - lamda12*r12(1)
      ff2y = ff2y - lamda02*r02(2) - lamda12*r12(2)
      ff2z = ff2z - lamda02*r02(3) - lamda12*r12(3)


      ff(i0)%x = ff0x 
      ff(i0)%y = ff0y 
      ff(i0)%z = ff0z 
     
      ff(i1)%x = ff1x 
      ff(i1)%y = ff1y
      ff(i1)%z = ff1z
     
      ff(i2)%x = ff2x 
      ff(i2)%y = ff2y 
      ff(i2)%z = ff2z 

      virialx = virialx + (r01(1)*r01(1)*lamda01)
      virialy = virialy + (r01(2)*r01(2)*lamda01)
      virialz = virialz + (r01(3)*r01(3)*lamda01)

      virialx = virialx + (r02(1)*r02(1)*lamda02)
      virialy = virialy + (r02(2)*r02(2)*lamda02)
      virialz = virialz + (r02(3)*r02(3)*lamda02)

      virialx = virialx + (r12(1)*r12(1)*lamda12)
      virialy = virialy + (r12(2)*r12(2)*lamda12)
      virialz = virialz + (r12(3)*r12(3)*lamda12)


     
    end subroutine shake3angle


    ! -----------------------------------------------------------------------
    ! SHAKE update to unconstrained coords
    ! xshake = predicted new coords using current velocity/force
    ! assumes NVE update, seems to be accurate enough for NVT,NPH,NPT as well

    subroutine shake_update(iflag)
      implicit none
      integer iflag
      integer i
      real*8 invmass,dtsq_factor

      if (iflag == 0) dtsq_factor = 0.5*dt*dt*kcal2amu
      if (iflag == 1) dtsq_factor = dt*dt*kcal2amu


      do i = 1,np
        invmass = 1.0/mass(position(i)%type)
        xshake(i)%x = position(i)%x + dt*v(i)%x + dtsq_factor*invmass*ff(i)%x
        xshake(i)%y = position(i)%y + dt*v(i)%y + dtsq_factor*invmass*ff(i)%y
        xshake(i)%z = position(i)%z + dt*v(i)%z + dtsq_factor*invmass*ff(i)%z

      enddo
      
      
    end subroutine shake_update




     
  end module mod_shake
  
