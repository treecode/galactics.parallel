      subroutine gendenspsibulge

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
      common /moreconstants/ v02, v03, rdisk2, diskconst
      common /cutoff/ psic_bulge, psic_halo, psi0_prime
      common /dfcutoff/ fcut_bulge, fcut_halo

      open(file='denspsibulge.dat',unit=50,status='replace')
      fcut_bulge = dfbulge(-psic_bulge)
      write(50,*) npsi,minlog

      dlogpsi = (log10((psic_bulge-psi0)/2)-minlog)/float(npsi/2)
      do i=1,npsi/2
         tablepsib(i) = minlog+dlogpsi*float(i)
         psi = -psi0 - 10.**tablepsib(i)
         call getbulgedens(psi,denspsibulge(i))
         denspsibulge(i) = log10(denspsibulge(i))
         write(50,*) tablepsib(i),denspsibulge(i),psi
      enddo
      
      do i=npsi/2-1,0,-1
         tablepsib(npsi-i) = minlog + dlogpsi*float(i)
         psi = -psic_bulge + 10.**tablepsib(npsi-i)
         call getbulgedens(psi,denspsibulge(npsi-i))
         denspsibulge(npsi-i) = log10(denspsibulge(npsi-i))
         write(50,*) tablepsib(npsi-i),denspsibulge(npsi-i),psi
      enddo
      close(50)
      return
      
      end

      subroutine getbulgedens(psi,rho)
      common /gparameters/  c, v0, a,
     +                      cbulge, v0bulge, abulge, 
     +                      psi0, haloconst, bulgeconst,
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /cutoff/ psic_bulge, psic_halo, psi0_prime
      common /dfcutoff/ fcut_bulge, fcut_halo
      
      real pi,psi
      pi = 3.1415927
      m = 257
      v02 = v0*v0
      vmin = -10.
      vmax = log(sqrt(2.*psi))
      dlogv = (vmax - vmin)/float(m-1)
      rho = 0.
      do j=1,4
         vlog = vmin + dlogv*float(j-1)
         v = exp(vlog)
         e = psi - v*v/2.
         df_lowered = dfbulge(e)-fcut_bulge
         if(df_lowered.gt.0.) 
     +         rho = rho + 4.*pi*dlogv*v*v*v*coef(j)*df_lowered
      enddo
      do j=5,m-4
         vlog = vmin + dlogv*float(j-1)
         v = exp(vlog)
         e = psi - v*v/2.
         df_lowered = dfbulge(e)-fcut_bulge
         if(df_lowered.gt.0.) 
     +         rho = rho + 4.*pi*dlogv*v*v*v*df_lowered
      enddo
      do jj=4,2,-1
         j = m-jj+1
         vlog = vmin + dlogv*float(j-1)
         v = exp(vlog)
         e = psi - v*v/2.
         df_lowered = dfbulge(e)-fcut_bulge
         if(df_lowered.gt.0.) 
     +         rho = rho + 4.*pi*dlogv*v*v*v*coef(jj)*df_lowered
      enddo
      return
      end

