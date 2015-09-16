      function densrpsi(r,psi_n)
      parameter (pi=3.1415926535, sq8=2.828427125)

      real psi,pi,psi00,psi_n,psicut
      common /potconstants/ apot(20,0:20000), fr(20,0:20000), 
     +     dr, nr, lmax, potcor
      common /gparameters/  a, b, c, v0, q, psi00, 
     +                      psiout, rho1, sigbulge2, 
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /moreconstants/ v02, v03, rdisk2, diskconst, bulgea
      common /cutoff/ psicut
      common /myconstants/ rho0,ra,eps,a00

      m = 33
      vmin = -10.
      psi = -psi_n
      if(psi.lt.0.) then
         rho=0.
         goto 99
      endif
      vmax = log(sqrt(2.*psi))
      dlogv = (vmax - vmin)/float(m-1)
      rho = 0.
      do j=1,4
         vlog = vmin + dlogv*float(j-1)
         v = exp(vlog)
         e = psi - v*v/2.
         rho = rho + 4.*pi*dlogv*v*v*v*df(e)*coef(j)
      enddo
      do j=5,m-4
         vlog = vmin + dlogv*float(j-1)
         v = exp(vlog)
         e = psi - v*v/2.
         rho = rho + 4.*pi*dlogv*v*v*v*df(e)
      enddo
      do jj=2,4
         j = m-jj+1
         vlog = vmin + dlogv*float(j-1)
         v = exp(vlog)
         e = psi - v*v/2.
         rho = rho + 4.*pi*dlogv*v*v*v*df(e)*coef(jj)
      enddo
 99   densrpsi = rho
      return
      end

      function df(e)
      real psiv02,psicut,pi
      common /gparameters/  a, b, c, v0, q, psi00, 
     +                      psiout, rho1, sigbulge2, 
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /moreconstants/ v02, v03, rdisk2, diskconst, bulgea
      common /cutoff/ psicut
      common /myconstants/ rho0,ra,eps,a00

      pi = 3.1415927
      ev02 = e/v02
      psiv02 = (v02+psi00)/v02
      if(e.le.0..or.ev02.ge.1-psiv02) then
         df = 0.
      else
         f1 = ev02**1.5
         f2 = (1-psiv02-ev02)**(-2.5)
         f3 = ((1-psiv02-ev02)/(-log(ev02+psiv02)))**2.715
         f4 = 0.362*ev02-0.5639*ev02*ev02-0.0859*(ev02**3.)-
     *        0.4912*(ev02**4.)
         df = 0.09197*a*f1*f3*f2*exp(f4)
clmw         f2 = (1-psiv02-ev02)**(-1.)
clmw         f4 = -1.02*ev02 + 23.8*ev02*ev02-98.3*(ev02**3.)+
clmw     *        194.5*(ev02**4.)-180*(ev02**5.)+63.3*(ev02**6.)
clmw         df = 0.0246*a*f1*f3*f2*exp(f4)
      endif      
      return
      end

      function coef(j)
      
      if(j.eq.1) coef=17./48.
      if(j.eq.2) coef=59./48.
      if(j.eq.3) coef=43./48.
      if(j.eq.4) coef=49./48.

      return
      end

