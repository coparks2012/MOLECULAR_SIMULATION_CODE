module mod_average
use global
contains
  subroutine average

    implicit none
    logical :: e
    integer :: dummy, i, n,io
    double precision :: dummyvar,avg,aux,var,aux1
    
    
    
    n=0
    avgE = 0.0d0
    varE = 0.0d0
    open(unit=1,file="PErun.dat")    
    do 
       read(1,*,IOSTAT=io)dummyvar,aux

       if(io.gt.0)then
          print*,'error in energy data,io gt 0'
          
       elseif(io.lt.0)then
          varE = sqrt(varE/(n*1.0d0)-(avgE/(n*1.0d0))**2)
          avgE = avgE/(n*1.d0)
          exit
          
       else
          avgE = aux+avgE
          varE = varE+aux**2
          n=n+1
       endif
    enddo
    close(1)


    avgT =0.0d0
    varT = 0.0d0
    n    =0
    open(unit=1,file="Temprun.dat")    
    do 
       read(1,*,IOSTAT=io)dummyvar,aux
       if(io.gt.0)then
          print*,'error in energy data,io gt 0'
          
       elseif(io.lt.0)then
          varT = sqrt(varT/(n*1.0d0)-(avgT/(n*1.0d0))**2)
          avgT = avgT/(n*1.d0)
          exit
          
       else
          avgT=aux+avgT
          varT = varT+aux**2
          n=n+1
       endif
    enddo
    close(1)

    avgP =0.0d0
    varP = 0.0d0
    n    =0
    open(unit=1,file="pressure_run.dat")    
    do 
       read(1,*,IOSTAT=io)dummyvar,aux
       if(io.gt.0)then
          print*,'error in energy data,io gt 0'
          
       elseif(io.lt.0)then
          varP = sqrt(varP/(n*1.0d0)-(avgP/(n*1.0d0))**2)
          avgP = avgP/(n*1.d0)
          exit
          
       else
          avgP=aux+avgP
          varP = varP+aux**2
          n=n+1
       endif
    enddo
    close(1)
    

    
      
  end subroutine average

end module mod_average
