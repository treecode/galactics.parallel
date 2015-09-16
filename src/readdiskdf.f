c to fix: interpolation in rcirc does not seem to work well. after 
c that, hope that things come together!

      subroutine readdiskdf(filename1,ibuf)
      parameter(nrmax=1000)
      parameter (pi=3.1415926535)
      integer*4 ibuf(15),jbuf(15),kbuf(15),lbuf(4)
      character*80 filename, filename1

      common /fileroot/ filename
      common /potconstants/ apot(20,0:80000), frad(20,0:80000), 
     +     dr, nr, lmax, potcor
      common /legendre/ plcon(0:40)
      common /gparameters/  c, v0, a,
     +                      cbulge, v0bulge, abulge, 
     +                      psi0, haloconst, bulgeconst,
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /moreconstants/ v02, v03, rdisk2, diskconst, bulgea
      common /diskpars/ sigr0, disksr, nrdisk
      common /splines/ rr(0:nrmax), fdrat(0:nrmax), drat2(0:nrmax), 
     +                 fszrat(0:nrmax), szrat2(0:nrmax), nrspl
      common /diskblackhole/ idiskbhflag

      idiskbhflag = 0
      filename1 = 'dbh.dat'
      call readharmfile(filename1,ibuf,jbuf,kbuf,lbuf)
      filename = filename1
c
c Read in the correction functions
c
      open(17,file='cordbh.dat',status='old')
      read(17,'(2x, 2g17.7,x,i4)') sigr0,disksr,nrspl
      do i=0,nrspl
         read(17,*) rr(i),fdrat(i),fszrat(i)
c         write(0,*) rr(i),fdrat(i),fszrat(i)
      enddo
      close(17)
      call splined(rr(0),fdrat(0),nrspl+1,1.e32,1.e32,drat2(0))
      call splined(rr(0),fszrat(0),nrspl+1,1.e32,1.e32,szrat2(0))

      rfid = 2.5*rdisk
      call omekap(rfid,fom,fka)
      sigr = sqrt(sigr2(rfid))
      sigden = diskdensf(rfid,0.0)*2.0*zdisk

c
c a close estimate of the disk surface density
c
      sigrcrit = 3.36*sigden/fka
      qtoomre = sigr/sigrcrit
c      write(0,*) 'Toomre Q = ',qtoomre, ' at R = 2.5 R_d'

c      open(12,file='toomre.dat',status='unknown')
c      do i=0, nrspl
c         call omekap(rr(i),fom,fka)
c         sigr = sqrt(sigr2(rr(i)))
c         sigden = diskdens(rr(i),0.0)*2.0*zdisk
c         sigrcrit = 3.36*sigden/fka
c         qtoomre = sigr/sigrcrit
c          write(12,*) rr(i), qtoomre
c      enddo
c      close(12)

      return

      end
