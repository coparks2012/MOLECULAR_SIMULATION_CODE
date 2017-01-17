module mod_qcprot
  use mod_matrix_multiply
  implicit none
  
contains
  
  !! Superposition coords2 onto coords1 -- in other words, coords2 is rotated, coords1 is held fixed */
  double precision function calcrmsd(coords1, coords2, len,rot,q1,q2,q3,q4)
    implicit none
    integer :: len
    double precision :: x,y,z,euc_dist
    double precision :: coords1(3,len),coords2(3,len) 
    double precision :: A(9),rot(9)
    double precision :: E0
    double precision :: rmsd
    double precision :: q1,q2,q3,q4
    integer :: i
    



    !---calculate the (weighted) inner product of two structures */
    E0 = IP(A, coords1, coords2,len)
    
    !---calculate the RMSD & rotational matrix */
    call FastCalcRMSDAndRotation(rot, A, rmsd, E0, len, -1.0d0,q1,q2,q3,q4)
        


    !---rotate second set of coordinates
    call matrix_multiply(coords2,len,rot)


    !---calculate euclidean distance */
     euc_dist = 0.0d0
     do  i = 1, len
        x = coords1(1,i)-coords2(1,i)
        y = coords1(2,i)-coords2(2,i)
        z = coords1(3,i)-coords2(3,i)
        
        euc_dist = euc_dist + x*x + y*y + z*z
     end do

    calcrmsd = sqrt(euc_dist/(1.0d0*len))

  end function calcrmsd


  
 
  

  subroutine centercoords(coords,len)
    implicit none
    integer :: len
    double precision :: coords(3,len)
    integer :: i
    double precision :: x,y,z

    x = 0.0d0; y = 0.0d0; z = 0.0d0
    do i = 1,len
       x = x + coords(1,i) 
       y = y + coords(2,i) 
       z = z + coords(3,i) 
    enddo

    x = x/real(len)
    y = y/real(len)
    z = z/real(len)

    coords(1,:) = coords(1,:) - x
    coords(2,:) = coords(2,:) - y
    coords(3,:) = coords(3,:) - z
  end subroutine centercoords
       


  double precision function IP(A, coords1,  coords2, len)
    implicit none
    integer :: len
    double precision :: coords1(3,len),coords2(3,len)
    double precision :: x1, x2, y1, y2, z1, z2, G1,G2
    double precision :: A(9)
    double precision :: fx1(len), fy1(len), fz1(len)
    double precision :: fx2(len), fy2(len), fz2(len)
    integer :: i
    integer :: weight_flag

    fx1(:) = coords1(1,:)
    fy1(:) = coords1(2,:)
    fz1(:) = coords1(3,:)


    fx2(:) = coords2(1,:)
    fy2(:) = coords2(2,:)
    fz2(:) = coords2(3,:)


    A(:) = 0.0d0

        
    do i = 1,len
       
       x1 = fx1(i)
       y1 = fy1(i)
       z1 = fz1(i)
       
       G1  = G1 +x1 * x1 + y1 * y1 + z1 * z1
       
       x2 = fx2(i)
       y2 = fy2(i)
       z2 = fz2(i) 
       
       G2 = G2 + (x2 * x2 + y2 * y2 + z2 * z2)
       
       A(1) = A(1) + (x1 * x2)
       A(2) = A(2) + (x1 * y2)
       A(3) = A(3) + (x1 * z2)
       
       A(4) = A(4) + (y1 * x2)
       A(5) = A(5) + (y1 * y2)
       A(6) = A(6) + (y1 * z2)
       
       A(7) = A(7) + (z1 * x2)
       A(8) = A(8) + (z1 * y2)
       A(9) = A(9) + (z1 * z2)  
    end do
    IP = (G1 + G2) * 0.50d0
  end function IP
  


  subroutine FastCalcRMSDAndRotation(rot,  A,  rmsd,  E0,  len,  minScore,q1,q2,q3,q4)
    implicit none
    integer :: len
    double precision :: rot(9), A(9)
    double precision :: Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz
    double precision :: Szz2, Syy2, Sxx2, Sxy2, Syz2, Sxz2, Syx2, Szy2, Szx2
    double precision :: SyzSzymSyySzz2, Sxx2Syy2Szz2Syz2Szy2, Sxy2Sxz2Syx2Szx2
    double precision :: SxzpSzx, SyzpSzy, SxypSyx, SyzmSzy
    double precision :: SxzmSzx, SxymSyx, SxxpSyy, SxxmSyy
    double precision :: C(0:3)
    double precision :: rmsd,E0,minScore
    integer:: i
    double precision :: mxEigenV
    double precision :: oldg = 0.0
    double precision :: b, aa, delta, rms, qsqr
    double precision :: q1, q2, q3, q4, normq
    double precision :: a11, a12, a13, a14, a21, a22, a23, a24
    double precision :: a31, a32, a33, a34, a41, a42, a43, a44
    double precision :: a2, x2, y2, z2
    double precision :: xy, az, zx, ay, yz, ax
    double precision :: a3344_4334, a3244_4234, a3243_4233, a3143_4133,a3144_4134, a3142_4132
    double precision :: evecprec
    double precision :: evalprec 
    double precision :: a1324_1423, a1223_1322 ,a1123_1321,a1122_1221
    double precision :: a1224_1422,a1124_1421
    evecprec = 1e-6; evalprec = 1e-11
    oldg = 0.0d0
    

    Sxx = A(1); Sxy = A(2); Sxz = A(3)
    Syx = A(4); Syy = A(5); Syz = A(6)
    Szx = A(7); Szy = A(8); Szz = A(9)

    Sxx2 = Sxx * Sxx
    Syy2 = Syy * Syy
    Szz2 = Szz * Szz

    Sxy2 = Sxy * Sxy
    Syz2 = Syz * Syz
    Sxz2 = Sxz * Sxz

    Syx2 = Syx * Syx
    Szy2 = Szy * Szy
    Szx2 = Szx * Szx

    SyzSzymSyySzz2 = 2.0d0*(Syz*Szy - Syy*Szz);
    Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2

    C(2) = -2.0d0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2)
    C(1) = 8.0d0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz)

    SxzpSzx = Sxz + Szx
    SyzpSzy = Syz + Szy
    SxypSyx = Sxy + Syx
    SyzmSzy = Syz - Szy
    SxzmSzx = Sxz - Szx
    SxymSyx = Sxy - Syx
    SxxpSyy = Sxx + Syy
    SxxmSyy = Sxx - Syy
    Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2

    C(0) = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2 &
         + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2) &
         + (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz)) &
         + (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz)) &
         + (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz)) &
         + (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz))

    !* Newton-Raphson */
    mxEigenV = E0;
     
    do i = 1,50
       oldg = mxEigenV
       x2 = mxEigenV*mxEigenV
       b = (x2 + C(2))*mxEigenV
       aa = b + C(1)
       delta = ((aa*mxEigenV + C(0))/(2.0*x2*mxEigenV + b + aa))
       mxEigenV = mxEigenv- delta
       if (abs(mxEigenV - oldg).lt.abs(evalprec*mxEigenV))then         
          exit
       endif
    enddo
       

    if (i.eq.50) then
       print*,'nMore than 50 iterations needed'
    end if

    !/* the fabs() is to guard against extremely small, but *negative* numbers due to floating point error */
    rms = sqrt(abs(2.0 * (E0 - mxEigenV)/len));
    rmsd = rms

    if(minScore.gt.0) then
       if (rms.lt.minScore)then
          return
       endif
    endif
 
    a11 = SxxpSyy + Szz-mxEigenV; a12 = SyzmSzy; a13 = - SxzmSzx; a14 = SxymSyx
    a21 = SyzmSzy; a22 = SxxmSyy - Szz-mxEigenV; a23 = SxypSyx; a24= SxzpSzx
    a31 = a13; a32 = a23; a33 = Syy-Sxx-Szz - mxEigenV; a34 = SyzpSzy
    a41 = a14; a42 = a24; a43 = a34; a44 = Szz - SxxpSyy - mxEigenV
    a3344_4334 = a33 * a44 - a43 * a34; a3244_4234 = a32 * a44-a42*a34
    a3243_4233 = a32 * a43 - a42 * a33; a3143_4133 = a31 * a43-a41*a33
    a3144_4134 = a31 * a44 - a41 * a34; a3142_4132 = a31 * a42-a41*a32
    q1 =  a22*a3344_4334-a23*a3244_4234+a24*a3243_4233
    q2 = -a21*a3344_4334+a23*a3144_4134-a24*a3143_4133
    q3 =  a21*a3244_4234-a22*a3144_4134+a24*a3142_4132
    q4 = -a21*a3243_4233+a22*a3143_4133-a23*a3142_4132

    qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4

    !The following code tries to calculate another column in the adjoint matrix when the norm of the 
    !current column is too small.
    !Usually this block will never be activated.  To be absolutely safe this should be
    !uncommented, but it is most likely unnecessary.

    if (qsqr.lt.evecprec)then
    
        q1 =  a12*a3344_4334 - a13*a3244_4234 + a14*a3243_4233
        q2 = -a11*a3344_4334 + a13*a3144_4134 - a14*a3143_4133
        q3 =  a11*a3244_4234 - a12*a3144_4134 + a14*a3142_4132
        q4 = -a11*a3243_4233 + a12*a3143_4133 - a13*a3142_4132
        qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4

        if (qsqr.lt.evecprec)then
        
           a1324_1423 = a13 * a24 - a14 * a23; a1224_1422 = a12 * a24 - a14 * a22
           a1223_1322 = a12 * a23 - a13 * a22; a1124_1421 = a11 * a24 - a14 * a21
           a1123_1321 = a11 * a23 - a13 * a21; a1122_1221 = a11 * a22 - a12 * a21

            q1 =  a42 * a1324_1423 - a43 * a1224_1422 + a44 * a1223_1322
            q2 = -a41 * a1324_1423 + a43 * a1124_1421 - a44 * a1123_1321
            q3 =  a41 * a1224_1422 - a42 * a1124_1421 + a44 * a1122_1221
            q4 = -a41 * a1223_1322 + a42 * a1123_1321 - a43 * a1122_1221
            qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4

            if (qsqr.lt.evecprec)then
            
                q1 =  a32 * a1324_1423 - a33 * a1224_1422 + a34 * a1223_1322
                q2 = -a31 * a1324_1423 + a33 * a1124_1421 - a34 * a1123_1321
                q3 =  a31 * a1224_1422 - a32 * a1124_1421 + a34 * a1122_1221
                q4 = -a31 * a1223_1322 + a32 * a1123_1321 - a33 * a1122_1221
                qsqr = q1*q1 + q2 *q2 + q3*q3 + q4*q4;
                
                if (qsqr.lt.evecprec)then
                
                   !if qsqr is still too small, return the identity matrix. */
                   rot(1) = 1.0
                   rot(2) = 0.0
                   rot(3) = 0.0
                   rot(4) = 0.0
                   rot(5) = 1.0
                   rot(6) = 0.0
                   rot(7) = 0.0
                   rot(8) = 0.0
                   rot(9) = 1.0

                    return
                 endif
              endif
           endif
        endif
                
            
     
    normq = sqrt(qsqr)
    q1 = q1/normq
    q2 = q2/normq
    q3 = q3/normq
    q4 = q4/normq

    a2 = q1 * q1
    x2 = q2 * q2
    y2 = q3 * q3
    z2 = q4 * q4

    xy = q2 * q3
    az = q1 * q4
    zx = q4 * q2
    ay = q1 * q3
    yz = q3 * q4
    ax = q1 * q2

    rot(1) = a2 + x2 - y2 - z2
    rot(2) = 2 * (xy + az)
    rot(3) = 2 * (zx - ay)
    rot(4) = 2 * (xy - az)
    rot(5) = a2 - x2 + y2 - z2
    rot(6) = 2 * (yz + ax)
    rot(7) = 2 * (zx + ay)
    rot(8) = 2 * (yz - ax)
    rot(9) = a2 - x2 - y2 + z2
  end subroutine FastCalcRMSDAndRotation

   
   
end module mod_qcprot

