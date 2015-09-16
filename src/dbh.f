c     Changes for 22-Jan-2003 (LMW)
c
c     This version of dbh.f will make a king model halo Mainly to see
c     how the code works before trying the nfw models
c     
c     CHANGES 25-FEB-1994: (KK) Core radius of halo is now calculated
c     from halo density only.  Before, the King radius was taken to
c     depend on the total density.  Makes no difference if coreparam is
c     set to zero.
c
c     R1 (the radius which appears for dimensional reasons in the halo
c     DF) has been reinstated (before it was set to 1.).  It is set as
c     follows: The radius ra is prompted for: it is such that the halo
c     rotation curve near r=0 is (r/ra)*v0 This radius is then used to
c     calculate the central density of the halo, which then gives
c     r1. (accounting only for the C-term with q=1 in the halo DF,
c     i.e. taking the isothermal sphere halo, and ignoring truncation
c     effects).
c
c     Convergence has been increased greatly by using a weighted mean of
c     old and new potential as start for the next step, rather than the
c     new potential only. .75*new+.25*old seems to be very effective in
c     damping waves which sometimes took a long time to settle down
c     before. The first time a particular harmonic is introduced its
c     full value is used, though.
c
c     Now also write out the second R-derivative of the potential, to
c     enable a good calculation of the epicycle frequency. This is read
c     by getfreqs3.
c
c     INCORPORATED WITH THE CHANGES OF JD 21-MAR-1994
c
c     CHANGES 21 MARCH 1994 (KK)
c  
c     disk cutoff radius implemented as a sharp cutoff in angular
c     momentum.  This corresponds to multiplying the disk density by a
c     complementary error function.  At every iteration the disk density
c     now changes slightly, though, since the density near the cutoff
c     depends on omega and kappa.
c
c     The output file now also includes a table of disk surface mass
c     density vs.  radius. The disk parameters as input are also to be
c     found there.  Getfreqs5 reads them OK.
c
c     CHANGES 11 MAY 1994 (KK):
c        
c     Changes the disk component to a dynamic one, which depends more on
c     the potential. Rather than enforcing a sech**2 disk, which is
c     difficult to realize selfconsistently in the presence of other
c     gravitating components (the disk is then no longer isothermal in
c     the vertical direction), the disk is now isothermal, with density
c     sur0/(2*zdisk)*exp(-r/rdisk)*exp(-(psi(r,z)-psi(r,0))/sigz2(r)).
c     This makes the vertical kinematics closer to isothermal.  The
c     vertical dispersion is chosen to give the correct half-density
c     scale height: sigz2=[psi(r,0)-psi(r,zdisk)]/ln[sech**2(1)] Instead
c     of writing the disk surface density, now write the disk midplane
c     density to the dbh.dat file.
c
c     Also write out the potential at every iteration on a series of 10
c     grid points, in the file 'refpoints.dat'. This can be used
c     afterwards to see how the different harmonics contribute.  The
c     grid lies in the disk equator, at half a scale height above it,
c     and at 2 scale height above it.
c
c     CHANGES FROM DBH6 TO DBH7: (KK 13 MAY 1994)
c
c     Try to reduce the load on the higher harmonics by subtracting off
c     a potential obtained by integrating a sech**2 law vertically.
c     Might give problems near the cutoff radius if the truncation width
c     is very narrow (on the order of the scale height).
c
c     Have also simplified the disk cutoff ---- simply use an error
c     function in radius. Makes life a lot easier for calculating the
c     density corresponding to the approximate disk potential...
c  
c     The parameter sigrtrunc has been turned into drtrunc, the gaussian
c     interval over which the truncation takes place.
c
c     FUNDAMENTAL SNAG: THE POTENTIAL HARMONICS NOW NO LONGER HAVE THE
c     CORRECT ASYMPTOTIC BEHAVIOUR AT INFINITY. HENCE THE SERIES
c     EXPANSION BECOMES KINDA BOGUS....
c
c     CHANGES 19 MAY 1994 (KK)
c
c     Can be remedied by writing the potential in spherical coordinates.
c     So psi=f(r) ln cosh r cos(theta)/zdisk which at large radii does
c     tend to zero.  f(r) again is such that the midplane density is
c     correct.
c
c     ALSO CHANGES THE FILE NAMES FOR THE OUTPUT, TO DBH.DAT, H.DAT AND
c     B.DAT
c
c     Output file dbh.dat now contains all input parameters, i.e. all
c     halo, disk and bulge parameters.

      program dbh

      common /potconstants/ apot(20,0:80000), fr(20,0:80000), 
     +     dr, nr, lmax, potcor
      common /legendre/ plcon(0:40)
      common /gparameters/  c, v0, a,
     +                      cbulge, v0bulge, abulge, 
     +                      psi0, haloconst, bulgeconst,
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /moreconstants/ v02, v03, rdisk2, diskconst, rhoh, rhob
      common /flags/ idiskflag, ibulgeflag, ihaloflag, ibhflag
      common /denspsiparameters/ npsi,minlog,
     +     tablepsi(1:1000),denspsihalo(1:1000),
     +     tablepsib(1:1000),denspsibulge(1:1000)
      common /cutoff/ psic_bulge, psic_halo, psi0_prime
      common /blackhole/ bhmass
      common /diskblackhole/ idiskbhflag

      parameter(pi=3.1415926535)
      real adens(20,0:80000),s1(0:80000),s2(0:80000)
      real rref(10),zref(10),pref(10),oldpref(10)
      real fr2(20,0:80000)
      real minlog,qc
      real psi0,psi0_prime
      real psic_halo,psic_halo_prime
      real psic_bulge,psic_bulge_prime 
      character ans*1
      
c  constants for legendre polynomials
      do l=0,40
         plcon(l) = sqrt((2*l+1)/(4.0*pi))
      enddo
c  set flags
      idiskflag = 0
      ibulgeflag = 0
      ihaloflag  = 0
      ibhflag = 0
      idiskbhflag = 0
c set default parameters
c disk: 
      rmdisk = 1.0
      rdisk  = 1.0
      zdisk  = .15
      outdisk = 5.0
      drtrunc = .3
c bulge:
      rho1 = 5
      psiout = -2.0
      sigbulge = .4
c halo:
      psi0 = -3
      q = 1.0
      v0 = 1.0
      rking = 1
c
c     Enter halo parameters
c
      write(*,*) 'Include halo?'
      read(*,'(a)') ans

      if(ans.eq.'y') then
         write(*,*) 'cutoff parameter, v0, a?'
         read(*,*) c, v0, a
         ihaloflag = 1
         v02 = v0*v0
         v03 = v0**3.
         haloconst = 3./(64.*sqrt(2.)*pi*pi*a*a*v0)
         rhoh = v02/(4.*pi*a*a)
      endif
      
c Enter disk parameters
      write(*,*) 'Include disk?'
      read(*,'(a)') ans
      if( ans .eq. 'y' ) then
         write(*,*)
     $        'Disk mass, scale length, radius, scale height, trunc width?'
         read(*,*) rmdisk, rdisk, outdisk, zdisk, drtrunc
         idiskflag = 1

         rdisk2 = rdisk*rdisk
         diskconst = rmdisk/(4.0*pi*rdisk2*zdisk)
         do iref=1,4
           rref(iref)=iref*rdisk*2
           zref(iref)=0
           enddo
         do iref=5,7
           rref(iref)=(iref-5)*rdisk*4
           zref(iref)=zdisk/2.
           enddo
         do iref=8,10
           rref(iref)=rref(iref-3)
           zref(iref)=zdisk*2.
           enddo
         open(18,file='refpoints.dat',status='unknown')
      else
         rmdisk = 0.
      endif

c Enter bulge parameters
      write(*,*) 'Include bulge?'
      read(*,'(a)') ans
      if( ans .eq. 'y' ) then
         write(*,*) 'Bulge cutoff, v0_bulge, a_bulge?'
         read(*,*) cbulge, v0bulge, abulge
         bulgeconst = 1./8./sqrt(2.)/(pi**3.)/abulge/abulge/v0bulge
         rhob = v0bulge*v0bulge/2./pi/abulge/abulge
         ibulgeflag = 1
      else
         cbulge = 0.
         v0bulge = 0.
         abulge = 1.
         ibulgeflag = 0
      endif

c     Enter blackhole parameters
      write(*,*) 'Include blackhole?'
      read(*,'(a)') ans
      if(ans.eq.'y') then
         write(*,*) 'blackhole mass?'
         read(*,*) bhmass
         ibhflag = 1
      else
         bhmass = 0.
         ibhflag = 0
      endif
cc
c     something new to try
cc
      
      v02disk = rmdisk/rdisk
c      v02disk = 0.
      psi0 = -v02 - v0bulge**2. - v02disk
      psi0_prime = -v02
      psic_halo_prime = -(1.-c)*v02
      call eprime_to_e_halo(-psic_halo_prime,psic_halo)
      psic_halo = -psic_halo

      if(ihaloflag.eq.0) psic_halo = 0.

      if(ibulgeflag.eq.1) then
         psic_bulge_prime = -(1.-cbulge)*v0bulge*v0bulge
         call eprime_to_e_bulge(-psic_bulge_prime,psic_bulge)
         psic_bulge = -psic_bulge
      endif

      if( ihaloflag .eq. 1 ) then 
         dr=1
         apot(1,0)=psi0*sqrt(4.*pi)
         lmax=0
         nr=2
      endif
c  cold start
      write(*,*) 'radial step size, no. of bins?'
      read(*,*) dr,nr
      write(*,*) 'max. azimuthal harmonic l?'
      read(*,*) lmaxx
c  initial potential is a spherical approximation to the unlowered Evans model.
      if( ihaloflag .eq. 0 ) then
         q = 1
         psi00 = -3.0
         rcmid = 1.0 
         v02 = 1.0
      endif
      z = 0.0
      do ir=0,nr
         r=ir*dr
         if(r.eq.0.) then
            apot(1,ir)=psi0*sqrt(4.*pi)
         else
            apot(1,ir)=psi0*sqrt(4.*pi)*a/r*log(1.+r/a)
         endif
         do l=2,lmaxx,2
            apot(l/2+1,ir)=0
         enddo
      enddo
      niter=(2+lmaxx/2)*10
      lmaxstep=2
      lmax=0
      ntheta=lmax*10+2
      ntheta=max(10,ntheta)

c     now iterate. number of iterations and change in lmax depends on
c     initial conditions.  iterate on first cycle until tidal radius is
c     stable, then once for every harmonic added, and until convergence
c     once l has reached the desired maximum lmax.
      
      drtidal = 2*dr
      rtidalold = 1e10
      tidalcheck = rtidalold
      lmaxold=-2
      iteroutside = 0
      open(file='in.gendenspsi',unit=40,status='old')
      read(40,*) minlog,npsi
      close(40)
      call gendenspsihalo
      if(ibulgeflag.eq.1) call gendenspsibulge
      open(file='check.dat',status='replace',unit=10)
      eps = dr/1.e4
      frac=0.75
      do iter=0,150
         if( lmax .eq. 0 .or. lmax .eq. lmaxx ) then
            if(drtidal .lt. dr .and. iter.gt.10) then
               lmax=lmax+lmaxstep
               ntheta=lmax*4+2
            endif
         else
            lmax=lmax+lmaxstep
            ntheta=lmax*4+2
         endif
         if( lmax .eq. lmaxx+lmaxstep ) goto 199

c     Now get the harmonics of the density in this potential --> adens
c     NB that dens only gives the density without an approximate sech**2
c     component --- that part of the potential is not represented by the
c     harmonics.  The function dens(r,z) returns the density -
c     high-frequency cpt The full density is returned by totdens

         adens(1,0)=dens(eps,0.)*sqrt(4.*pi)
         do l=2,lmax,2
            adens(l/2+1,0)=0
         enddo
         do ir=1,nr
            adens(1,ir)=0
         enddo
c     nrmax will mark the outermost radial bin with non-zero density.
         nrmax=nr
         do l=0,lmax,2
c     integrate density * spherical harmonic function over quadrant use
c     cos(theta) as independent variable.
            do ir=1,nrmax
               r=ir*dr
               s=0
               dctheta=1.0/ntheta
               s=s+polardens(r,1.0,l)+polardens(r,0.0,l)
               do is=1,ntheta-1,2
                  ctheta=is*dctheta
                  s=s+4*polardens(r,ctheta,l)
               enddo
               do is=2,ntheta-2,2
                  ctheta=is*dctheta
                  s=s+2*polardens(r,ctheta,l)
               enddo
               s=s*dctheta/3.
               s=s*4*pi
               adens(l/2+1,ir)=s               
c     mark the first even radial bin on which the density has fallen to
c     zero.
               if (l.eq.0 .and. s.eq.0.) then
                  nrmax=nr
                  nrmax=nrmax+mod(nrmax,2)
                  goto 77
               endif
            enddo
 77      enddo
c     now get the potential harmonics of this new density. (BT 2-208)
c     Simpson's rule integration.
         do l=0,lmax,2
            s1(0)=0
            r = 2*dr
            s1(2)=(r*dr/3.)*(4*adens(l/2+1,1)*(1.0-dr/r)**(l+2)
     +           +adens(l/2+1,2))
            rold = r 
            do ir=4,nr,2
               r=ir*dr
               s1a = (r*dr/3.)*(adens(l/2+1,ir-2)*(1.0-2*dr/r)**(l+2)+
     &              4*adens(l/2+1,ir-1)*(1.0-dr/r)**(l+2)+adens(l/2+1,ir))
               s1(ir) = s1a + s1(ir-2)*(rold/r)**(l+1)
               rold = r
            enddo
            s2(nr)=0
            rold = nr*dr
            do ir=nr-2,2,-2
               r=ir*dr
               s2a = (r*dr/3.)*(adens(l/2+1,ir+2)*(1.0+2*dr/r)**(1-l)+
     &              4*adens(l/2+1,ir+1)*(1.0+dr/r)**(1-l)+adens(l/2+1,ir))
               s2(ir) = s2a + s2(ir+2)*(r/rold)**l
               rold = r
            enddo

c     replace the potential harmonics with a mean of the previous
c     iteration (25%) and the current one (75%). This damps out
c     oscillations that otherwise occur.  if this is the first time this
c     harmonic is calculated, use the entire new value.
            do ir=2,nr,2
               if (l.le.lmaxold) then
                  apot(l/2+1,ir)=frac*apot(l/2+1,ir)-
     +                 (1.-frac)*4*pi/(2.*l+1.)*(s1(ir)+s2(ir))
               else
                  apot(l/2+1,ir)=-4*pi/(2.*l+1.)*(s1(ir)+s2(ir))
               endif
            enddo
c     Calculate the 1st and 2nd-order radial gradients
  
            do ir=2,nr,2
               r = ir*dr
               fr(l/2+1,ir)=-4*pi/(2.*l+1.)*(-(l+1)*s1(ir) + l*s2(ir))/r
               fr2(l/2+1,ir)=-4*pi/(2.*l+1.)*
     +              ((l+1)*(l+2)*s1(ir)/r**2+ 
     +              l*(l-1)*s2(ir)/r**2 -(2*l+1)*adens(l/2+1,ir))
            enddo
         enddo
c     now interpolate the gaps first quadratically interpolate the
c     monopole back to the origin.  the remaining multipoles are zero
c     there.
         apot(1,0)=3*(apot(1,2)-apot(1,4))+apot(1,6)
         fr(1,0)=0.0
         fr2(1,0)=2*fr2(1,2)-fr2(1,4)
         do l=2,lmax,2
            apot(l/2+1,0)=0
            fr(l/2+1,0)=0
            fr2(l/2+1,0)=0
         enddo
c  then linearly interpolate other bins.
         do ir=1,nr-1,2
            do l=0,lmax,2
               apot(l/2+1,ir)=(apot(l/2+1,ir-1)+apot(l/2+1,ir+1))/2.
               fr(l/2+1,ir)=(fr(l/2+1,ir-1)+fr(l/2+1,ir+1))/2.
               fr2(l/2+1,ir)=(fr2(l/2+1,ir-1)+fr2(l/2+1,ir+1))/2.
            enddo
         enddo
c  finally reset the potential at the origin to psi0
c  Note that the fake disk potential is zero at the origin.

         if (ihaloflag.eq.0.and.ibulgeflag.eq.0) goto 199

         a00=apot(1,0)
         open(file='test.dat',unit=10) ! ,status='new')
         do ir=0,nr
c            write(10,*) ir*dr,apot(1,ir)/sqrt(4.*pi),adens(1,ir)
            apot(1,ir)=apot(1,ir)+psi0*sqrt(4.*pi)-a00
         enddo
         if (a00/sqrt(4.*pi)-psi0.gt.-psic_halo) then
            write(*,'(''Iter'',i4,'': lmax='',i4,
     +           '', tidal radius is infinite'')') iter, lmax 
            iteroutside = iteroutside + 1
            if( iteroutside .gt. 80 ) then
                write(*,'(''nr='',i4,'' is too small'',
     +           '', try larger number of radial bins  - exiting program'')') nr
                goto 12345
            endif
            drtidal = 2.0*dr
         else
            potr=psi0
            do ir=1,nr
c     look for new tidal radius of the model defined as the radius at
c     the equator where the potential is zero
               potrm1=potr
               potr=pot(ir*dr,0.)
               aa = potrm1 - psic_halo
               bb = potr - psic_halo
               if (aa*bb.le.0.) then
                  dpot = potr - potrm1
                  if( dpot .eq. 0.0 ) then
                     rtidal = (ir-1)*dr
                  else
                     rtidal=(ir-1-aa/dpot)*dr
                  endif
                  drtidal = abs(rtidal - rtidalold)
                  tidalcheck = abs(rtidal - rtidalold)/rtidal
                  write(*,'(''Iter'',i4,'': lmax='',i4,
     +                 '', tidal radius is '',g15.6)') iter,lmax,rtidal
                  rtidalold = rtidal
                  goto 9
               endif
            enddo
            write(*,'(''Iter'',i4,'': lmax='',i4,
     +           '', tidal radius is outside grid'')') iter,lmax
            drtidal = 2.0*dr
 9       endif
c write out the changes in the potential at the reference points 
c at this iteration.
         if( idiskflag .eq. 1 ) then 
            do iref=1,10
               oldpref(iref)=pref(iref)
               pref(iref)=pot(rref(iref),zref(iref))
            enddo
            write(18,'(2i3,10g12.4)') iter,lmax,
     +           (pref(iref)-oldpref(iref),iref=1,10)
         endif
c     now repeat with this new potential!
         lmaxold=lmax
      enddo
c     end of the main loop
 199  continue
      call check

      
      open(file='rtidal.dat',unit=20,status='replace')
      if( idiskflag .eq. 1 ) close(18)
      totalmass = fr(1,nr)/sqrt(4*pi)*(dr*nr)**2
      write(*,*) 'Total mass=', totalmass
      write(20,*) rtidal
c  
c  Calculate force and potential for halo only
c  
      lmax = lmaxx
      halomass = 0.0
      bulgemass = 0.0
      haloedge = 0
      bulgeedge = 0
      write(*,*) 'ihaloflag,ibulgeflag=', ihaloflag,ibulgeflag
      if( ihaloflag .eq. 1 ) call halopotential(halomass, haloedge)
      if( ibulgeflag .eq. 1 ) THEN 
         call bulgepotential(bulgemass, bulgeedge)
      else
C           create an empty b.dat file so that make doesn't complain
         open(11,file='b.dat',status='unknown')
         close(11)
      endif
      if( ibhflag. eq. 1) call genblackhole(bhmass)

      diskmass = totalmass - halomass - bulgemass
      diskedge = outdisk + 2.0*drtrunc

      open(20,file='mr.dat',status='unknown')
      write(20,*) diskmass, diskedge
      write(20,*) bulgemass, bulgeedge
      write(20,*) halomass, haloedge
      close(20)

      open(11,file='dbh.dat',status='unknown')
      write(11,'('' # c,v0,a,cbulge,v0bulge,abulge,dr,nr,lmax='')')
      write(11,'('' #'',7g15.5,i6,i4)') c,v0,a,cbulge,v0bulge,abulge,
     +     dr,nr,lmaxx
      write(11,'('' # psi0, haloconst, bulgeconst:'')')
      write(11,'('' #'',3g15.5)') psi0,haloconst,bulgeconst
      write(11,'('' # Mdisk, rdisk, zdisk, outdisk, drtrunc'')')
      write(11,'('' #'',5g15.5)') rmdisk, rdisk, zdisk, outdisk, drtrunc
      write(11,'('' # psic_bulge,psic_halo,psi0_prime,bhmass'')')
      write(11,'('' #'',5g15.5)') psic_bulge,psic_halo,psi0_prime,bhmass
      write(11,'('' #'',4i5)') idiskflag, ibulgeflag, ihaloflag, ibhflag
      write(11,'('' #  OUTPUT FROM DBH8. TOTAL POTENTIAL.'')')

      do ir=0,nr
         write(11,'(8g16.8)') ir*dr,(adens(l/2+1,ir),l=0,lmaxx,2)
      enddo
      write(11,*)
      write(11,*)
      do ir=0,nr
         write(11,'(8g16.8)') ir*dr,(apot(l/2+1,ir),l=0,lmaxx,2)
      enddo
      write(11,*)
      write(11,*)
      do ir=0,nr
         write(11,'(8g16.8)') ir*dr,(fr(l/2+1,ir),l=0,lmaxx,2)
      enddo
      write(11,*)
      write(11,*)
      do ir=0,nr
         write(11,'(8g16.8)') ir*dr,(fr2(l/2+1,ir),l=0,lmaxx,2)
      enddo
      write(11,*)
      do ir=0,nr
         psi=pot(ir*dr,0.)
         write(11,'(2g16.8)') ir*dr,diskdens(ir*dr,0.,psi)
      enddo
      close(11)

      write(*,*) 'Final model written to file ''dbh.dat''.'
      call check
12345 continue

      open(file='tidalr.out',status='replace',unit=14)
      write(14,*) rtidal
      
      end

      subroutine check

      common /potconstants/ apot(20,0:80000), fr(20,0:80000), 
     +     dr, nr, lmax, potcor
      common /legendre/ plcon(0:40)
      common /gparameters/  c, v0, a,
     +                      cbulge, v0bulge, abulge, 
     +                      psi0, haloconst, bulgeconst,
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /moreconstants/ v02, v03, rdisk2, diskconst, rhoh, rhob
      common /flags/ idiskflag, ibulgeflag, ihaloflag, ibhflag

      z = 0.
      ebmass = 0.
      bmass = 0.
      open(file='check.dat',unit=10,status='replace')
      do ir=1,nr
         r = ir*dr
         psi = pot(r,z)
         ehdens = rhoh*a/r/((1.+r/a)**2.)
         hdens = halodens(r,z)
         bdens = bulgedenspsi(psi)
         tdens = totdens(r,z)
         ddens = diskdens(r,z,psi)
         call force(r,z,vc,fz,rpot)
         vc = sqrt(-vc*r)
         s = r/a
         w = v0*sqrt(log(1+s)/s - 1/(1+s))
         if(hdens.gt.0.) then
            write(10,88) log10(r),log10(ddens),log10(bdens),
     +           log10(hdens),log10(tdens),log10(ebdens),log10(ehdens)
c            write(10,87) log10(r),log10(hdens),log10(ehdens),r,vc,w
         endif
      enddo
 88   format(7f12.3)
 87   format(6f12.3)
      close(10)
      return
      end

      subroutine genblackhole(bhmass)

      open(file='blackhole',unit=20,status='replace')
      n = 1
      t = 0.
      write(20,*) n,t
      write(20,*) bhmass,t,t,t,t,t,t
      close(20)
      
      return
      end

