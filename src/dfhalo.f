      function dfhalo(e)

      common /gparameters/  c, v0, a,
     +                      cbulge, v0bulge, abulge, 
     +                      psi0, haloconst, bulgeconst,
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /cutoff/ psic_bulge, psic_halo, psi0_prime
      common /flags/ idiskflag, ibulgeflag, ihaloflag, ibhflag

      if(e.lt.-psic_halo.or.e.gt.-psi0) then
         dfhalo = 0.
      else
         if(ibulgeflag.eq.1.or.idiskflag.eq.1) then
            call e_to_eprime_halo(e,eprime)
         else
            eprime = e
         endif
         ev02 = eprime/v0/v0
         f1 = ev02**1.5
         f2 = (1-ev02)**(-2.5)
         f3 = ((1-ev02)/(-log(ev02)))**2.715
         f4 = 0.779+0.362*ev02-0.5639*ev02*ev02-0.0859*(ev02**3.)-
     *        0.4912*(ev02**4.)
         dfhalo = haloconst*f1*f2*f3*exp(f4)
      endif      
      return
      end



