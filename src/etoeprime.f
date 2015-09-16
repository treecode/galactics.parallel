      subroutine eprime_to_e_halo(eprime,e)

      common /gparameters/  c, v0, a,
     +                      cbulge, v0bulge, abulge, 
     +                      psi0, haloconst, bulgeconst,
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /cutoff/ psic_bulge, psic_halo, psi0_prime

      aa = -psi0
      bb = -psi0_prime
      alpha = (aa/bb-1)/bb
      e = eprime*(1.+alpha*eprime)

      return
      end

      subroutine e_to_eprime_halo(e,eprime)

      common /gparameters/  c, v0, a,
     +                      cbulge, v0bulge, abulge, 
     +                      psi0, haloconst, bulgeconst,
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /cutoff/ psic_bulge, psic_halo, psi0_prime

      aa = -psi0
      bb = -psi0_prime
      alpha = (aa/bb-1)/bb
      eprime = (sqrt(1.+4.*e*alpha)-1)/2./alpha      
      return
      end

      subroutine eprime_to_e_bulge(eprime,e)

      common /gparameters/  c, v0, a,
     +                      cbulge, v0bulge, abulge, 
     +                      psi0, haloconst, bulgeconst,
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1

      e = eprime - psi0 - v0bulge**2.

      return
      end

      subroutine e_to_eprime_bulge(e,eprime)

      common /gparameters/  c, v0, a,
     +                      cbulge, v0bulge, abulge, 
     +                      psi0, haloconst, bulgeconst,
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      
      eprime = e + psi0 + v0bulge**2.

      return
      end




