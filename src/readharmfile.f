c     this file contains the following subroutines:
c
c     readharmfile: reads a file with the harmonic expansion of the
c     potential and density.  all the relevant stuff ends up in a common
c     block that only needs to be seen by the routines in this file (I
c     think!). This routine must be called first.
c
c     dens(r,z): density at point (r,z)
c
c     densrpsi(r,psi): density as a function of r and potential psi.
c
c     pot(r,z): potential at point (r,z)
c
c     all the rest is (slightly nobbled in places) numerical recipes
c     routines.
c

      subroutine readharmfile(filename,ibuf1,jbuf1,kbuf1,lbuf1)
      parameter(pi=3.1415926535)
      common /potconstants/ apot(20,0:80000), frad(20,0:80000), 
     +     dr, nr, lmax, potcor
      common /legendre/ plcon(0:40)
      common /gparameters/  c, v0, a,
     +                      cbulge, v0bulge, abulge, 
     +                      psi0, haloconst, bulgeconst,
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /moreconstants/ v02, v03, rdisk2, diskconst
      common /diskpars/ sigr0, disksr, nrdisk
      common /flags/ idiskflag, ibulgeflag, ihaloflag, ibhflag
      common /cutoff/ psic_bulge, psic_halo, psi0_prime
      common /blackhole/ bhmass

      real adens(20,0:80000)
      integer*4 ibuf(15), ibuf1(15)
      integer*4 jbuf(3), jbuf1(3)
      integer*4 kbuf(1), kbuf1(1)
      integer*4 lbuf(4), lbuf1(4)
c this allows the gparameters to be passes as a C structure
      character*30 filename
      equivalence (c, ibuf)
      equivalence (psic_bulge, jbuf)
      equivalence (bhmass, kbuf)
      equivalence (idiskflag, lbuf)

c  constants for legendre polynomials
      do l=0,40
         plcon(l) = sqrt((2*l+1)/(4.0*pi))
      enddo
c
      open(14,file=filename,status='old')
      read(14,*)
      read(14,'(2x,7g15.5,i6,i4)') c,v0,a,cbulge,v0bulge,abulge,
     +     dr,nr,lmax
      read(14,*)
      read(14,'(2x,3g15.5)') psi0, haloconst, bulgeconst
      read(14,*)
      read(14,'(2x,5g15.5)') rmdisk, rdisk, zdisk, outdisk, drtrunc
      read(14,*)
      read(14,'(2x,4g15.5)') psic_bulge,psic_halo,psi0_prime,bhmass
      read(14,'(2x,4i5)') idiskflag, ibulgeflag, ihaloflag, ibhflag
      read(14,*)
      do ir=0,nr
         read(14,'(8g16.8)') rr,(adens(l/2+1,ir),l=0,lmax,2)
         enddo
      read(14,*)
      read(14,*)
      do ir=0,nr
         read(14,'(8g16.8)') rr,(apot(l/2+1,ir),l=0,lmax,2)
         enddo
      read(14,*)
      read(14,*)
      do ir=0,nr
         read(14,'(8g16.8)') rr,(frad(l/2+1,ir),l=0,lmax,2)
         enddo
      close(14)
c
c these are needed for the bulgeconstants common block
c
      sigbulge2=sigbulge*sigbulge
c
c additional disk constants for efficiency - common diskconstants
c
      nrdisk= int((outdisk + 2.0*drtrunc)/dr) + 10
      rdisk2 = rdisk*rdisk
      diskconst = rmdisk/(4.0*pi*rdisk2*zdisk)
c
c and some more constants...
c
      v02=v0*v0
      v03=v0*v02
c
c calculate potential correction
c
      redge = nr*dr
      potcor = 0
      do l=0,lmax,2
          potcor = potcor + apot(l/2+1,nr) + frad(l/2+1,nr)*redge/(l+1)
      enddo
      potcor = potcor*plcon(0)
      potcor1 = potcor
c
c transfer gparams buffer for output
c

      do i=1,15
         ibuf1(i) = ibuf(i)
      enddo
      do j=1,3
         jbuf1(j) = jbuf(j)
      enddo
      do k=1,1
         kbuf1(k) = kbuf(k)
      enddo
      do k=1,4
         lbuf1(k) = lbuf(k)
      enddo
        
      return
      end
