      program color_test
      
      use color_def
      use color_lib

      implicit none

      integer :: nct
      integer, parameter :: i110=65, i160=75, i390=12, i606=25, i814=39
      double precision :: feh, afe, logL, logT, logg, filter(maxnc)
      logical,parameter :: do_something_else = .false.
      character(len=128) :: test_data_dir

      test_data_dir = '/home/dotter/science/Spectra'

      call interactive
      
      contains

      subroutine interactive
            
      logT=99d0
      do while(logT > 0d0)
         write(*,*)  'enter nct, FeH, afe, logT, logg, logL:'
         read(*,*) nct, FeH,afe, logT, logg, logL

         call color_init(test_data_dir,afe,nct)
         call get_mags(nct, feh, afe, logL, logg, logT, filter)

         write(*,'(99a12)') 'logT', 'logg', 'logL', filter_name(nct,1:ncol(nct))
         write(*,'(99f12.4)') logT, logg, logL, filter(1:ncol(nct))
      enddo

      end subroutine interactive

      end program color_test
