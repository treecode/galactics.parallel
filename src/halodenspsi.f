      function halodenspsi(psi_n)
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

      dlogpsi = (log10((-psi0+psic_halo)/2)-minlog)/float(npsi/2)
      psi = -psi_n
      if(psi.gt.-psi0.or.psi.lt.-psic_halo) then
         halodenspsi = 0.
         return
      endif
      if(psi.le.(-psic_halo-psi0)/2) then
         j = int((log10(psic_halo+psi)-minlog)/dlogpsi)
         if(j.lt.1) j=1
         frac = (log10(psi+psic_halo)-tablepsi(j))
     +        /(tablepsi(j+1)-tablepsi(j))
         del = denspsihalo(j+1)-denspsihalo(j)
      else 
         j = npsi-1-int((log10(-psi0-psi)-minlog)/dlogpsi)
         if(j.ge.npsi .or. j .le. 0) j=npsi-1
         frac = (log10(-psi0-psi)-tablepsi(j))
     +        /(tablepsi(j+1)-tablepsi(j))
         del = denspsihalo(j+1)-denspsihalo(j)
      endif
      halodenspsi = 10.**(denspsihalo(j)+frac*del)
      return
      end

