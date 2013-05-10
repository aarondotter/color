CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       PHOENIX COLOR INTERPOLATION IN TEFF, LOG G, [Fe/H],& [a/Fe]    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

ccccccccccccccccccccccc subroutine phx_color_init cccccccccccccccccccccc
c this subroutine determines which filter set will be used, reads in   c
c the appropriate tables, and stores them in coltbl                    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module phoenix
      
      use color_def

      implicit none

      private
      integer, parameter :: nz_max=10
      integer, parameter :: nt=66, nfl=14, nmax=1000, phx_file=21
      double precision, parameter :: cln = dlog(1.0d1)
      double precision :: z(nz_max),coltbl(maxct,maxnc,nmax,nz_max)
      double precision :: teff(nz_max,nmax),ggl(nz_max,nmax)
      integer :: itemp(nz_max,nt),itnum(nz_max),isize(nz_max),pcol(nfl),nz
      logical :: phx_debug = .false.
      character(len=128) :: phx_data_dir

      public :: phx_color_init, phx_color_interp

      contains

      subroutine phx_color_init(data_dir,afe,nct)
      character(len=128) :: data_dir
      double precision, intent(in) :: afe
      integer, intent(in) :: nct
      integer :: iafe

      phx_data_dir = trim(data_dir) // '/phx/'

      z = [-4.0d0,-3.5d0,-3.0d0,-2.5d0,-2.0d0,-1.5d0,-1.0d0,-0.5d0,0.0d0,0.5d0]

      if(afe>0.2d0) then
         nz=9
      else
         nz=10
      endif

! set [a/Fe] index parameter
      iafe=0
      if(afe==-0.2d0) iafe=1
      if(afe== 0.0d0) iafe=2
      if(afe== 0.2d0) iafe=3
      if(afe== 0.4d0) iafe=4
      if(afe== 0.6d0) iafe=5
      if(afe== 0.8d0) iafe=6
      !if(feh<=-2.5d0) iafe=2
      if(iafe==0) stop "PROBLEM WITH [a/Fe] IN ISOCHRONE FILE"

c read in the required tables for a range of [Fe/H] at fixed [a/Fe]    
      call phx_color_read(iafe,nct)

c all set for interpolation in T_eff and log(g)      
      end subroutine phx_color_init


ccccccccccccccccccccccccc subroutine phx_color_read ccccccccccccccccccccccc
c reads in color tables based on specified Z, [a/Fe], filter set          c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine phx_color_read(iafe,nct)

      integer :: iafe,nct,iz,ia,j,ierr
      double precision :: feh0
      character(len=128) :: filename
      character(len=5) :: zfile(nz_max)
      character(len=3) :: afile(6)

      afile = ['am2','ap0','ap2','ap4','ap6','ap8']
      zfile = ['Zm4d0','Zm3d5','Zm3d0','Zm2d5','Zm2d0','Zm1d5','Zm1d0','Zm0d5','Zp0d0','Zp0d5']

      do iz=1,nz
c set up filenames based on the input variables
         if(iz<4)then
            ia = 2
         else
            ia=iafe
         endif
         filename = trim(phx_data_dir) // zfile(iz) // afile(ia) // "." // suffix(nct)
         if(phx_debug) write(0,*) trim(filename)
c open table for reading
         open(phx_file,file=trim(filename),status='old')
         read(phx_file,*) !skip header
         read(phx_file,'(17x,99a12)',iostat=ierr) filter_name(nct,1:ncol(nct))
         read(phx_file,*,iostat=ierr) !skip this line
c read data in table
         do j=1,nmax
            read(phx_file,*,iostat=ierr) teff(iz,j),ggl(iz,j),feh0,coltbl(nct,1:ncol(nct),j,iz)
            if(ierr/=0) exit
            if(feh0/=z(iz)) write(0,*) 'WARNING: table [Fe/H] incorrect'
         enddo
c close data file
         close(phx_file)
c isize counts number of data points, itnum counts number of distinct T's,
c itemp gives location of each T value in the teff-array
         isize(iz)=j-1
         itemp(iz,1)=1
         itnum(iz)=1
         do j=2,isize(iz)
            if(teff(iz,j)>teff(iz,j-1)) then
               itnum(iz)=itnum(iz)+1
               itemp(iz,itnum(iz))=j
            endif
         enddo
      enddo
      
      end subroutine phx_color_read

ccccccccccccccccccccccc subroutine phx_color_interp ccccccccccccccccccccc
c     interpolates the mags and colors based on T, g, [Fe/H], [a/Fe]    c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine phx_color_interp(nct,feh,tl0,gl0,color)
      integer, intent(in) :: nct
      double precision, intent(in) :: feh, tl0, gl0
      double precision, intent(out) :: color(maxnc)
      integer :: i,j,k,iii,jj,inewt,inewg,iz,izlo,izhi,tintp,gintp,zintp
      double precision :: fg(4),ft(4),tl,gl,t,qt(4),qg(4),colr(nz_max,maxnc),col(nz_max,maxnc,4),fz(4)
      double precision :: feh0

      iz=0
      tl=tl0
      gl=gl0
      t=dexp(cln*tl)
      colr=0d0
      if(gl>=5.5d0) gl=5.4999d0
      !if(gl<=0d0) gl=0d0
      inewt=0
      inewg=0
      tintp=4
      gintp=4
      zintp=4

      feh0 = min(feh,z(nz)-1d-3)
      if(phx_debug) write(0,*) 'phx: ', feh, feh0, z(nz)

c     [Fe/H] interpolation
      do i=1,nz-1
         if( feh0>=z(i) .and. feh0<z(i+1) ) iz=i
      enddo
      if(iz<zintp/2) iz=zintp/2
      if(iz>nz-zintp/2) iz=nz-zintp/2

      izlo=iz+1-zintp/2
      izhi=iz+zintp/2
      call interp(z(izlo:izhi),fz,feh0,zintp)

      do iii=izlo,izhi
c     locate T in the Teff array (default is linear interpolation in T)
         do i=1,itnum(iii)-1
            if(t>=teff(iii,itemp(iii,i)).and.t<teff(iii,itemp(iii,i+1))) inewt=i
         enddo
         if(inewt<tintp/2) then
            tintp=2
            inewt=tintp/2
         elseif(inewt>itnum(iii)-tintp/2) then
            tintp=2
            inewt=itnum(iii)-tintp/2
         endif
c     find interpolation coeff.'s in T
         do i=1,tintp
            qt(i)=teff(iii,itemp(iii,inewt+i-tintp/2))
         enddo
         call interp(qt,ft,t,tintp)

         if(phx_debug)then
            write(*,*) 'iz, izlo, izhi, z(iz)=', iz, izlo, izhi, z(iz)
            write(*,'(a20,3f12.4)') 'Teff, logg, [Fe/H]=', t, gl, feh
            write(*,'(a20,1p4e12.4)') 'fz(:)=', fz
            write(*,'(a20,1p4e12.4)') 'qt(:)=', qt
            write(*,'(a20,1p4e12.4)') 'ft(:)=', ft
         endif


c  locate Log G for each T (default is quadratic interpolation in log g)
         do j=1,tintp
            jj=inewt+j-tintp/2
            do k=itemp(iii,jj),itemp(iii,jj+1)-1
               if(gl>=ggl(iii,k).and.gl<ggl(iii,k+1)) inewg=k
            enddo

            if(inewg<itemp(iii,jj)+gintp/2)then
               gintp=2
               inewg=itemp(iii,jj)+gintp/2
            elseif(inewg>itemp(iii,jj+1)-1-gintp/2) then
               gintp=2
               inewg=itemp(iii,jj+1)-1-gintp/2
            endif

c     find interpolation coefficients in Log G
            do k=1,gintp
               qg(k)=ggl(iii,inewg+k-gintp/2)
            enddo
            call interp(qg,fg,gl,gintp)

c     g-interpolation
            do k=1,ncol(nct)
               col(iii,k,j)=0.0d0
               do i=1,gintp
                  col(iii,k,j)=col(iii,k,j) + 
     *                 fg(i)*coltbl(nct,k,inewg+i-gintp/2,iii)
               enddo            !i-loop
            enddo               !k-loop
         enddo                  !j-loop

c     T-interpolation
         do i=1,ncol(nct)
            colr(iii,i)=0.0d0
            do j=1,tintp
               colr(iii,i)=colr(iii,i)+ft(j)*col(iii,i,j)
            enddo               !j-loop
         enddo                  !i-loop       
      enddo                     !iii-loop
      
c     Z-interpolation
      do i=1,ncol(nct)
         color(i)=0.0d0
         do j=1,zintp
            color(i)=color(i)+fz(j)*colr(iz-zintp/2+j,i)
         enddo
      enddo

      if(phx_debug)then
         do i=1,ncol(nct)
            write(*,'(a20,i3,5f12.6)') 'i, color(i)=', i, color(i), colr(izlo:izhi,i)
         enddo
      endif

      end subroutine phx_color_interp

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                           END of PHX_COLOR.F                          c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end module phoenix
