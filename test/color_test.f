      program color_test
      
      use colors

      implicit none

      integer :: nct
      integer, parameter :: i110=65, i160=75, i390=12, i606=25, i814=39
      integer, parameter :: col=90
      double precision :: feh, afe, logL, logT, logg, filter(col)
      logical,parameter :: do_something_else = .false.
      character(len=128) :: test_data_dir

      test_data_dir = '/home/dotter/science/Spectra'
      nct=SDSS
      afe=0.2d0

      call color_init(test_data_dir,nct,afe)

      call interactive
      
      contains

      subroutine interactive
            
      logT=99d0
      do while(logT > 0d0)
         write(*,*)  'enter FeH, logT, logg, logL:'
         read(*,*) FeH, logT, logg, logL
         call get_mags(nct, feh, afe, logL, logg, logt, filter)
         write(*,'(99a8)') 'logT', 'logg', 'logL', phx_filters(1:ncol)
         write(*,'(99f8.4)') logT, logg, logL, filter(1:ncol)
      enddo

      end subroutine interactive

      end program color_test
