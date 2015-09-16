      subroutine omekapfactors(R,omega,kappa,omegafac,kappafac)

      common /blackhole/ bhmass
      real omega,kappa,omegafac,kappafac

      omegafac = 1.
      kappafac = 1.
      if(bhmass.gt.0.) then
         omegafac = sqrt(1 + bhmass/(R**3.)/(omega**2.))
         kappafac = sqrt(1 - 2.*bhmass/(R**3.)/(kappa**2.))
      endif
      return
      end
      
