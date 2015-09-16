c--------------------------------
      function sigz2(r)
      parameter(nrmax=1000)
      common /gparameters/  c, v0, a,
     +                      cbulge, v0bulge, abulge, 
     +                      psi0, haloconst, bulgeconst,
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /moreconstants/ v02, v03, rdisk2, diskconst
      common /diskpars/ sigr0, disksr, nrdisk
      common /splines/ rr(0:nrmax),fdrat(0:nrmax),drat2(0:nrmax),
     +              fszrat(0:nrmax), szrat2(0:nrmax), nrspl
      psizh=pot(r,zdisk)
      psi0=pot(r,0.0)
      truesigz2=(psi0-psizh)/log(0.419974)

      call splintd(rr(0),fszrat(0),szrat2(0),nrspl+1,r,fcor)
      sigz2=truesigz2*fcor
      if(sigz2.lt.0.) sigz2 = 1.e-10
      return
      end
