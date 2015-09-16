      parameter (mf=100)
      real rf(mf),vf(mf),vc(mf),sig(mf)
      character*60 filename
      integer*4 ibuf1(15)
      
      open(10,file='new_rotation.dat',status='old')
      do i=1,999
         read(10,*,end=111) rf(i),vf(i),err,sig(i)
         rf(i) = rf(i)/4.465
      enddo
 111  nf=i-1
      close(10)
      write(*,*) nf,' data points in rotation curve'

      filename='dbh.dat'
      call readharmfile(filename,ibuf1)
      do i=1,nf
         call force(rf(i),0.,vc(i),fz,p)
         vc(i)=sqrt(-vc(i)*rf(i))
      enddo
      top=0
      bot=0
      do i=1,nf
         top=top+vc(i)*vf(i)
         bot=bot+vc(i)*vc(i)
      enddo
      vscale=top/bot
      rms=0
      sum = 0.
      do i=1,nf
         rms=rms+((vc(i)*vscale -vf(i))**2)/sig(i)/sig(i)
      enddo
      rms=sqrt(rms/nf)

c write out masses in Msun  
c                 (assumes radial unit is 1 kpc, rotation curve in km/s)
      open(10,file='mr.dat',status='old')
      read(10,*) dm,dr, bm,br, hm,hr
      close(10)

      dm=dm*vscale**2/4.3e-6
      hm=hm*vscale**2/4.3e-6
      bm=bm*vscale**2/4.3e-6

      write(*,*) 'RMS VC fit, disk M, bulge M, halo M, halo R, v scale'
      write(*,'(f6.3,3e10.3,f8.1,f9.2)') rms,dm,bm,hm,hr,vscale

      end
