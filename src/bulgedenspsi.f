      function bulgedenspsi(psi_n)
      real psi,psi_n,minlog
      common /denspsiparameters/ npsi,minlog,
     +     tablepsi(1:1000),denspsihalo(1:1000),
     +     tablepsib(1:1000),denspsibulge(1:1000)
      common /gparameters/  c, v0, a,
     +                      cbulge, v0bulge, abulge, 
     +                      psi0, haloconst, bulgeconst,
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /cutoff/ psic_bulge, psic_halo, psi0_prime
      
      psi = -psi_n

      dlogpsi = (log10((psic_bulge-psi0)/2)-minlog)/float(npsi/2)
      if(psi.lt.-psic_bulge) then
         bulgedenspsi = 0.
         return
      endif
      if(psi.ge.abs(-psic_bulge-psi0)/2) then
         if(-psi-psi0.le.0.) then
             psi = -psi0*(1.-1.e-7)
         endif
         j = int((log10(-psi-psi0)-minlog)/dlogpsi)
         if(j.lt.1) j=1
         frac = (log10(-psi-psi0)-tablepsib(j))
     +        /(tablepsib(j+1)-tablepsib(j))
         del = denspsibulge(j+1)-denspsibulge(j)
      else
         j = npsi-int((log10(psi+psic_bulge)-minlog)/dlogpsi)
         if(j.ge.npsi.or.j.le.0) j=npsi-1
         frac = (log10(psi+psic_bulge)-tablepsib(j))
     +        /(tablepsib(j+1)-tablepsib(j))
         del = denspsibulge(j+1)-denspsibulge(j)
      endif
      bulgedenspsi = 10.**(denspsibulge(j)+frac*del)

      return
      end

