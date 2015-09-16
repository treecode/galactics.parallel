      function dfbulge(e)
      real qc,q
      common /gparameters/  c, v0, a,
     +                      cbulge, v0bulge, abulge, 
     +                      psi0, haloconst, bulgeconst,
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1

      dfbulge = 0.
      call e_to_eprime_bulge(e,eprime)
      if(eprime.lt.0.) return
      q = sqrt(eprime/v0bulge/v0bulge)
      if(q.ge.1) return
      temp = 8.*(q**4.)-8.*q*q-3.
      df = 3.*asin(q) + q*sqrt(1.-q*q)*
     +     (1.-2.*q*q)*temp
      dfbulge = bulgeconst/((1.-q*q)**2.5)*df        
      return
      end

