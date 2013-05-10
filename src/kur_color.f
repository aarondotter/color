CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  CASTELLI & KURUCZ COLOR INTERPOLATION IN TEFF, LOG G, [Fe/H],& [a/Fe] C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

ccccccccccccccccccccccc subroutine kur_color_init cccccccccccccccccccccccc
c this subroutine determines which filter set will be used, reads in     c
c the appropriate tables, and stores them in coltbl                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module castelli_kurucz

      use color_def

      implicit none
      
      private
      integer, parameter :: nz=9
      integer, parameter :: nt=100, nfl=14, nmax=1000, kur_file=22
      double precision, parameter :: cln=dlog(1.0d1)
      double precision :: z(nz),coltbl(maxct,maxnc,nmax,nz),teff(nz,nmax),ggl(nz,nmax)
      integer :: itemp(nz,nt),itnum(nz),isize(nz),kcol(nfl)
      logical :: kur_debug = .false.
      character(len=128) :: kur_data_dir

      public :: kur_color_init, kur_color_interp

      contains

      subroutine kur_color_init(data_dir,afe,nct)
      character(len=128) :: data_dir
      double precision, intent(in) :: afe
      integer, intent(in) :: nct 
      integer :: iafe

      kur_data_dir = trim(data_dir) // '/kur/'

      z = [-4.0d0,-2.5d0,-2.0d0,-1.5d0,-1.0d0,-0.5d0,0.0d0,0.2d0,0.5d0]

c set [a/Fe] index parameter
      if(afe<=0.2d0) then
         iafe=1
      else
         iafe=2
      endif

c read in the required tables for a range of [Fe/H] at fixed [a/Fe]    
      call kur_color_read(iafe,nct)

c all set for interpolation in T_eff and log(g)
      end subroutine kur_color_init


ccccccccccccccccccccccccc subroutine kur_color_read ccccccccccccccccccccccc
c reads in color tables based on specified Z, [a/Fe], filter set          c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine kur_color_read(iafe,nct)

      integer :: iafe,nct,iz,j,ia,ierr
      character(len=128) :: filename
      character(len=5) :: zfile(nz)
      character(len=3) :: afile(2)
      double precision :: feh0

      afile = ['ap0','ap4']
      zfile = ['Zm4d0','Zm2d5','Zm2d0','Zm1d5','Zm1d0','Zm0d5','Zp0d0','Zp0d2','Zp0d5']

      do iz=1,nz
c set up filenames based on the input variables
         if(iz==1) then
            ia = 2
         else
            ia=iafe
         endif

         filename = trim(kur_data_dir) // zfile(iz) // afile(ia) // "." // trim(suffix(nct))
         if(kur_debug) write(0,*) filename
c open table for reading
         open(kur_file,file=trim(filename),status='old')
         read(kur_file,*) !skip header
         read(kur_file,'(17x,99a12)',iostat=ierr) filter_name(nct,1:ncol(nct))
         read(kur_file,*) !skip header
c read data in table
         do j=1,nmax
            read(kur_file,*,iostat=ierr) teff(iz,j),ggl(iz,j),feh0,coltbl(nct,1:ncol(nct),j,iz)
            if(ierr/=0) exit
            if(feh0/=z(iz)) write(0,*) 'WARNING table [Fe/H] incorrect'
         enddo
c close data file
         close(kur_file)
c isize counts number of data points, itnum counts number of distinct T's,
c itemp gives location of each T value in the teff-array
         isize(iz)=j-1
         itemp(iz,1)=1
         itnum(iz)=1
         do j=2,isize(iz)
            if(teff(iz,j).gt.teff(iz,j-1)) then
               itnum(iz)=itnum(iz)+1
               itemp(iz,itnum(iz))=j
            endif
         enddo
      enddo

      end subroutine kur_color_read

ccccccccccccccccccccccc subroutine kur_color_interp ccccccccccccccccccccccc
c     interpolates the mags and colors based on T, g, [Fe/H], [a/Fe]      c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine kur_color_interp(nct,feh,tl0,gl0,color)
      integer, intent(in) :: nct
      double precision, intent(in) :: feh, tl0, gl0
      double precision, intent(out) :: color(maxnc)
      integer :: i,j,k,iii,jj,inewt,inewg,iz,izlo,izhi,tintk,gintk,zintk
      double precision :: fg(4),ft(4),tl,gl,t,qt(4),qg(4),colr(nz,maxnc),col(nz,maxnc,4),fz(4)
      double precision :: feh0

      iz=0
      tl=tl0
      gl=gl0

      t=dexp(cln*tl)
      if(t.ge.4.9d4) t=4.8999d4
      if(gl.ge.5.5d0) gl=5.4999d0
      inewt=0
      inewg=0
      tintk=4
      gintk=4
      zintk=4

      feh0 = min(feh,z(nz)-1d-3)
      if(kur_debug) write(0,*) 'kur: ', feh, feh0, z(nz)

      do i=1,nz-1
         if( feh0 >=z(i) .and. feh0 <z(i+1) ) iz=i
      enddo
      if(iz<2) iz=zintk/2
      if(iz>nz-2) iz=nz-zintk/2

      izlo=iz+1-zintk/2
      izhi=iz+zintk/2
      call interp(z(izlo:izhi),fz,feh0,zintk)

      do iii=izlo,izhi
c     locate T in the Teff array
         inewt=0
         do i=1,itnum(iii)-1
            if(t.ge.teff(iii,itemp(iii,i)).and.
     *           t.lt.teff(iii,itemp(iii,i+1))) inewt=i
         enddo
         if(tintk.gt.itnum(iii)) tintk=itnum(iii)
         if(inewt.lt.2) inewt=tintk/2
         if(inewt.gt.itnum(iii)-tintk/2) inewt=itnum(iii)-tintk/2

c     find interpolation coeff.'s in T
         do i=1,tintk
            qt(i)=teff(iii,itemp(iii,inewt+i-tintk/2))
         enddo
         call interp(qt,ft,t,tintk)

         if(kur_debug)then
            write(*,*) 'Teff, logg, [Fe/H]=', t, gl, feh
            write(*,*) 'fz(:)=', fz
            write(*,*) 'ft(:)=', ft
         endif

c     locate Log G for each T
         do j=1,tintk
            jj=inewt+j-tintk/2
            inewg=0
            do k=itemp(iii,jj),itemp(iii,jj+1)-1
               if(gl.ge.ggl(iii,k)) inewg=k
            enddo

            if(inewg.lt.itemp(iii,jj)+gintk/2) inewg=itemp(iii,jj)+gintk/2
            if(inewg.gt.itemp(iii,jj+1)-1-gintk/2) 
     *           inewg=itemp(iii,jj+1)-1-gintk/2

c     find interpolation coefficients in Log G
            do k=1,gintk
               qg(k)=ggl(iii,inewg+k-gintk/2)
            enddo
            call interp(qg,fg,gl,gintk)

c     g-interpolation
            do k=1,ncol(nct)
               col(iii,k,j)=0.0d0
               do i=1,gintk
                  col(iii,k,j) = col(iii,k,j) + 
     *                 fg(i)*coltbl(nct,k,inewg+i-gintk/2,iii)
               enddo            !i-loop
            enddo               !k-loop
         enddo                  !j-loop
         
c     T-interpolation
         do i=1,ncol(nct)
            colr(iii,i)=0.0d0
            do j=1,tintk
               colr(iii,i)=colr(iii,i)+ft(j)*col(iii,i,j)
            enddo               !j-loop
         enddo                  !i-loop
      enddo                     !iii-loop
      
c     Z-interpolation
      do i=1,ncol(nct)
         color(i)=0.0d0
         do j=1,zintk
            color(i)=color(i)+fz(j)*colr(iz-zintk/2+j,i)
         enddo
      enddo

      end subroutine kur_color_interp

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                           END of KUR_COLOR.F                           c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end module castelli_kurucz
