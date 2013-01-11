cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                 ALL_COLOR.F                              c
c Serves as a holding place for code related to Teff-Color transformations c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module colors

      use phoenix
      use castelli_kurucz

      implicit none

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
      integer :: ncol, icol(maxct)
      double precision :: vegamag(10)
      character(len=12) :: filter_name(maxnc)
      character(len=128) :: data_dir
      logical :: colors_initialized = .false.

      contains

      subroutine color_init(my_data_dir,nct,afe)
      character(len=128), intent(in) :: my_data_dir
      double precision, intent(in) :: afe
      integer, intent(in) :: nct
      
      data_dir = my_data_dir

      call phx_color_init(data_dir, afe, nct)
      call kur_color_init(data_dir, afe, nct)

      icol(:) = (/10,6,13,12,17,77,5,5,7,4,5,4,6,6/)
      ncol=icol(nct)

      filter_name(1:ncol) = phx_filters(1:ncol)

!     from Bessell, Castelli, & Plez 1998 (BCP98) for Vega:
!     V     U-B    B-V    V-R    V-I
!     0.03 -0.004 -0.002 -0.007 -0.003
!                  U       B       V     R       I       
      vegamag = (/ 2.6d-2, 2.8d-2, 3d-2, 3.7d-2, 3.3d-2, 0d0, 0d0, 0d0, 0d0, 0d0 /)

      colors_initialized = .true.
      end subroutine color_init
     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_mags(nct,feh,afe,logL,grav,teff,filter)

      integer, intent(in) :: nct
      double precision, intent(in) :: feh, afe, logL, grav, teff
      double precision, intent(out) :: filter(maxnc)
      double precision :: color(maxnc),clr1(maxnc),Fhi,Flo,mbol
      double precision :: colorp(maxnc),colork(maxnc),clr2(maxnc)
                     ! alternatively tlo=3.95d0, thi=3.97d0
      double precision, parameter :: tlo=3.87d0, thi=3.93d0
      integer :: j

      if(.not.colors_initialized) stop 'get_mags: need to call color_init'
      !'
      mbol = solbol - 2.5d0*logL

!***************************
!for PHOENIX + Castelli & Kurucz synthetic
!***************************
c  smooth out bumps due to PHOENIX - Kurucz transition by applying two
c  adjustments: first, ramp between PHX and C&K between TLO and THI
c  and second, apply an offset to all mags calculated at THI for all
c  mags above THI
      if(teff>=tlo)then
         call phx_color_interp(feh,tlo,grav,colorp)
         call kur_color_interp(feh,tlo,grav,colork)
         do j=1,ncol            ! for the low-T end
            clr1(j)=colork(j)-colorp(j)
         enddo
         call phx_color_interp(feh,thi,grav,colorp)
         call kur_color_interp(feh,thi,grav,colork)
         do j=1,ncol            !for the high-T end
            clr2(j)=colork(j)-colorp(j)
         enddo
      endif
c actual Teff and log g for each point are found here
      if(teff<=thi) call phx_color_interp(feh,teff,grav,colorp)
      if(teff>=tlo) call kur_color_interp(feh,teff,grav,colork)

c phoenix colors only go up to 10,000 K, so ramp between PHX and KUR
c for high temperatures, the ramp is achieved with Flo and Fhi:
      Flo=(thi-teff)/(thi-tlo)  !goes from 1 at tlo to 0 at thi
      Flo = 0.5d0*(1d0-cos(pi*Flo)) !smoothe the edges, a la Bill Paxton
      Fhi = 1.0d0-Flo           !goes from 0 at tlo to 1 at thi

      do j=1,ncol               !loop through colors
         if(teff<tlo) then
c use only PHX colors below TLO
            color(j)=colorp(j)
         elseif(teff>=tlo.and.teff<=thi) then
c in this transitional region, the color is made up of 3 pieces:
c     #1-PHX ramping down from 1 at TLO to 0 at THI
c     #2-C&K ramping  up  from 0 at TLO to 1 at THI
c     #3-offset between PHX and C&K
            color(j)=Flo*colorp(j) + Fhi*(colork(j)-Flo*clr1(j)-Fhi*clr2(j))
         elseif(teff>thi)then
c use only C&K colors above THI but include offset as found at THI for continuity
            color(j) = colork(j) - clr2(j) 
         endif
      enddo

      forall(j=1:ncol) filter(j)=mbol-color(j)

      !add here any non-zero mags from definition of photometric system
      if(nct==UBVRIJHKs)then    ! for UBVRI we take Vega mags from BCP98
         forall(j=1:ncol) filter(j)=filter(j)+vegamag(j)
      endif

      end subroutine get_mags

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                       END OF ALL_COLOR FILE                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end module colors
