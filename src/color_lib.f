!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                 ALL_COLOR.F                              c
! Serves as a holding place for code related to Teff-Color transformations c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module color_lib

      use color_def
      use phoenix
      use castelli_kurucz

      implicit none

      contains

      subroutine color_init(data_dir, afe, nct)
      character(len=128), intent(in) :: data_dir
      double precision, intent(in) :: afe
      integer, intent(in) :: nct

      if(colors_initialized(nct)) return

      ncol = [10,10,13,12,17,77,5,5,7,4,5,4,6,6]

      suffix = [ 'UBVRIJHKsKp', 'WashDDOuvby', 'HST_WFPC2  ', 'HST_ACSWF  ',
     >            'HST_ACSHR  ', 'HST_WFC3   ', 'CFHTugriz  ', 'SDSSugriz  ',
     >            'PanSTARRS  ', 'SPITZER    ', 'UKIDSS     ', 'WISE       ',
     >            'SkyMapper  ', 'LSST       '                              ]
      
      !initialize sub-modules'
      call phx_color_init(data_dir, afe, nct)
      call kur_color_init(data_dir, afe, nct)

!     from Bessell, Castelli, & Plez 1998 (BCP98) for Vega:
!     V     U-B    B-V    V-R    V-I
!     0.03 -0.004 -0.002 -0.007 -0.003
!                  U       B       V     R       I      J     H    Ks   Kp  D51 
      vegamag = [ 2.6d-2, 2.8d-2, 3d-2, 3.7d-2, 3.3d-2, 0d0, 0d0, 0d0, 0d0, 0d0 ]

      colors_initialized(nct) = .true.
      end subroutine color_init
     
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_mags(nct,feh,afe,logL,grav,teff,filter)
      integer, intent(in) :: nct
      double precision, intent(in) :: feh, afe, logL, grav, teff
      double precision, intent(out) :: filter(maxnc)
      double precision :: color(maxnc),clr1(maxnc),Fhi,Flo,mbol
      double precision :: colorp(maxnc),colork(maxnc),clr2(maxnc)
                     ! alternatively tlo=3.95d0, thi=3.97d0
      double precision, parameter :: tlo=3.90d0, thi=3.93d0
      integer :: j
      double precision :: my_alpha_Fe

      my_alpha_Fe = afe

      if(.not.colors_initialized(nct)) stop 'get_mags: need to call color_init'
      !'
      mbol = solbol - 2.5d0*logL

!***************************
!for PHOENIX + Castelli & Kurucz synthetic
!***************************
!  smooth out bumps due to PHOENIX - Kurucz transition by applying two
!  adjustments: first, ramp between PHX and C&K between TLO and THI
!  and second, apply an offset to all mags calculated at THI for all
!  mags above THI
      if(teff>=tlo)then
         call phx_color_interp(nct,feh,tlo,grav,colorp)
         call kur_color_interp(nct,feh,tlo,grav,colork)
         do j=1,ncol(nct)            ! for the low-T end
            clr1(j)=colork(j)-colorp(j)
         enddo
         call phx_color_interp(nct,feh,thi,grav,colorp)
         call kur_color_interp(nct,feh,thi,grav,colork)
         do j=1,ncol(nct)            !for the high-T end
            clr2(j)=colork(j)-colorp(j)
         enddo
      endif
! actual Teff and log g for each point are found here
      if(teff<=thi) call phx_color_interp(nct,feh,teff,grav,colorp)
      if(teff>=tlo) call kur_color_interp(nct,feh,teff,grav,colork)

! phoenix colors only go up to 10,000 K, so ramp between PHX and KUR
! for high temperatures, the ramp is achieved with Flo and Fhi:
      Flo=(thi-teff)/(thi-tlo)  !goes from 1 at tlo to 0 at thi
      Flo = 0.5d0*(1d0-cos(pi*Flo)) !smoothe the edges, a la Bill Paxton
      Fhi = 1.0d0-Flo           !goes from 0 at tlo to 1 at thi

      do j=1,ncol(nct)          !loop through colors
         if(teff<tlo) then
! use only PHX colors below TLO
            color(j)=colorp(j)
         elseif(teff>=tlo.and.teff<=thi) then
! in this transitional region, the color is made up of 3 pieces:
!     #1-PHX ramping down from 1 at TLO to 0 at THI
!     #2-C&K ramping  up  from 0 at TLO to 1 at THI
!     #3-offset between PHX and C&K
            color(j)=Flo*colorp(j) + Fhi*(colork(j)-Flo*clr1(j)-Fhi*clr2(j))
         elseif(teff>thi)then
! use only C&K colors above THI but include offset as found at THI for continuity
            color(j) = colork(j) - clr2(j) 
         endif
      enddo

      forall(j=1:ncol(nct)) filter(j)=mbol-color(j)

      !add here any non-zero mags from definition of photometric system
      if(nct==UBVRIJHKs)then    ! for UBVRI we take Vega mags from BCP98
         forall(j=1:ncol(nct)) filter(j)=filter(j)+vegamag(j)
      endif

      end subroutine get_mags

      end module color_lib
