      subroutine gendenspsihalo
      parameter (pi=3.1415926535, sq8=2.828427125)

      real psi,psi_n,minlog
      common /denspsiparameters/ npsi,minlog,
     +     tablepsi(1:1000),denspsihalo(1:1000),
     +     tablepsib(1:1000),denspsibulge(1:1000)
      common /potconstants/ apot(20,0:20000), fr(20,0:20000), 
     +     dr, nr, lmax, potcor
      common /gparameters/  c, v0, a,
     +                      cbulge, v0bulge, abulge, 
     +                      psi0, haloconst, bulgeconst,
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /cutoff/ psic_bulge, psic_halo, psi0_prime
      common /dfcutoff/ fcut_bulge, fcut_halo

      open(file='denspsihalo.dat',unit=50,status='replace')
      write(50,*) npsi,minlog
      dlogpsi = (log10((-psi0+psic_halo)/2)-minlog)/float(npsi/2)
      fcut_halo = dfhalo(-psic_halo)
      do i=1,npsi/2
         tablepsi(i) = minlog+dlogpsi*float(i)
         psi = -psic_halo + 10.**tablepsi(i)
         call getdens(psi,denspsihalo(i))
         denspsihalo(i) = log10(denspsihalo(i))
         write(50,*) tablepsi(i),denspsihalo(i),psi
      enddo

      do i=npsi/2-1,0,-1
         tablepsi(npsi-i) = minlog + dlogpsi*float(i)
         psi = -psi0 - 10.**tablepsi(npsi-i)
         call getdens(psi,denspsihalo(npsi-i))
         denspsihalo(npsi-i) = log10(denspsihalo(npsi-i))
         write(50,*) tablepsi(npsi-i),denspsihalo(npsi-i),psi
      enddo
      close(50)
      return
      end

      subroutine getdens(psi,rho)
      common /gparameters/  c, v0, a,
     +                      cbulge, v0bulge, abulge, 
     +                      psi0, haloconst, bulgeconst,
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /cutoff/ psic_bulge, psic_halo, psi0_prime
      common /dfcutoff/ fcut_bulge, fcut_halo

      real pi,psi
      pi = 3.1415927
      m = 129
      vmin = -10.
      vmax = log(sqrt(2.*psi))
      dlogv = (vmax - vmin)/float(m-1)
      rho = 0.
      do j=1,4
         vlog = vmin + dlogv*float(j-1)
         v = exp(vlog)
         e = psi - v*v/2.
         dfl = dfhalo(e) - fcut_halo
         if(dfl.gt.0.) rho = rho + 4.*pi*dlogv*v*v*v*dfl*coef(j)
 77   enddo
      do j=5,m-4
         vlog = vmin + dlogv*float(j-1)
         v = exp(vlog)
         e = psi - v*v/2.
         if(e.lt.-psic_halo) goto 88
         dfl = dfhalo(e) - fcut_halo
         if(dfl.gt.0.) rho = rho + 4.*pi*dlogv*v*v*v*dfl
 88   enddo
      do jj=2,4
         j = m-jj+1
         vlog = vmin + dlogv*float(j-1)
         v = exp(vlog)
         e = psi - v*v/2.
         if(e.lt.-psic_halo) goto 99
         dfl = dfhalo(e) - fcut_halo
         if(dfl.gt.0.) rho = rho + 4.*pi*dlogv*v*v*v*dfl*coef(jj)
 99   enddo
      return
      end

      function coef(j)
      
      if(j.eq.1) coef=17./48.
      if(j.eq.2) coef=59./48.
      if(j.eq.3) coef=43./48.
      if(j.eq.4) coef=49./48.

      return
      end


