      module color_def

      !constants
      double precision, parameter :: solbol = 4.75d0
      double precision, parameter :: pi=3.14159265358979d0
      
      !definitions
      integer, parameter :: UBVRIJHKs=1
      integer, parameter :: WashDDO=2
      integer, parameter :: HST_WFPC2=3
      integer, parameter :: HST_ACS_WFC=4
      integer, parameter :: HST_ACS_HRC=5
      integer, parameter :: HST_WFC3=6
      integer, parameter :: CFHT=7
      integer, parameter :: SDSS=8
      integer, parameter :: PanSTARRS=9
      integer, parameter :: SPITZER=10
      integer, parameter :: UKIDSS=11
      integer, parameter :: WISE=12
      integer, parameter :: SkyMapper=13
      integer, parameter :: LSST=14
      integer, parameter :: maxct = LSST

      integer, parameter :: maxnc = 90
      integer :: ncol(maxct)
      character(len=12) :: filter_name(maxct,maxnc)
      logical :: colors_initialized(maxct) = .false.

      contains

      subroutine interp(a,b,x,n)
c {a} are the tabulated values for use in interpolation
c {b} are coefficients of the interpolating polynomial
c  x  is the abscissa to be interpolated
c  n  is the number of points to be used, interpolating polynomial
c     has order n-1 
      integer, intent(in) :: n
      double precision, intent(in) :: a(n), x
      double precision, intent(out) :: b(n)
      integer :: i,j
      do i=1,n
         b(i)=1.0d0
         do j=1,n
            if(j.ne.i) b(i)=b(i)*(x-a(j))/(a(i)-a(j))
         enddo
      enddo

      end subroutine interp

      end module color_def
