      subroutine dfmaximumhalo(emax,dfmax)

      common /gparameters/  c, v0, a,
     +                      cbulge, v0bulge, abulge, 
     +                      psi0, haloconst, bulgeconst,
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /cutoff/ psic_bulge, psic_halo, psi0_prime
      common /blackhole/ bhmass
      
      real psi

      dfmax = 0.
      emax = -psi0
      if(bhmass.le.0.) return
      ne = 1000
      de = (-psi0 + psic_halo)/ne
      call dfcorrectionhalo(-psic_halo,fac)
      fcut = dfhalo(-psic_halo)
      do i=1,ne
         e = -psic_halo + de*float(i)
         call dfcorrectionhalo(e,fac)
         df = (dfhalo(e)-fcut)/fac
         if(df.gt.dfmax) then
            dfmax = df
            emax = e - de
         endif
      enddo
      dfmax = 1.05*dfmax
      
      return
      end


      subroutine dfmaximumbulge(emax,dfmax)

      common /gparameters/  c, v0, a,
     +                      cbulge, v0bulge, abulge, 
     +                      psi0, haloconst, bulgeconst,
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /cutoff/ psic_bulge, psic_halo, psi0_prime
      common /blackhole/ bhmass
      real psi

      
      dfmax = 0.
      emax = -psi0
      if(bhmass.le.0.) return
      ne = 1000
      de = (-psi0 + psic_bulge)/ne
      call dfcorrectionbulge(-psic_bulge,fac)
      fcut = dfbulge(-psic_bulge)
      do i=1,ne
         e = -psic_bulge + de*float(i)
         call dfcorrectionbulge(e,fac)
         df = (dfbulge(e)-fcut)/fac
         if(df.gt.dfmax) then
            dfmax = df
            emax = e - de
         endif
      enddo
      dfmax = 1.05*dfmax
      
      return
      end

