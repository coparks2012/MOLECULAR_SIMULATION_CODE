module mod_print_results
  use global
  use mod_average
  implicit none

  contains 
    subroutine print_results

      print*,'--------RUN RESULTS---------------'
      print*,'                final potential run',potential*reps
      print*,'                final temperaturte run:',temp*rtemp
      
      call average
      

      print*,'                Average  Energy    :', avgE

      print*,'                Average  Energy per particle   :', avgE/np
      print*,'                Variance Energy per particle  :',  varE/np
      
      print*,'                Average temp :', avgT
      print*,'                Variance temp:', varT

      print*,'                Average Pres:',avgP
      print*,'                Var  P      :',varP
      
      
      print*,'--------END RUN RESULTS---------------'

    end subroutine print_results

  end module mod_print_results
