      subroutine readdenspsibulge

      real minlog
      common /denspsiparameters/ npsi,minlog,
     +     tablepsi(1:1000),denspsihalo(1:1000),
     +     tablepsib(1:1000),denspsibulge(1:1000)
      common /diskblackhole/ idiskbhflag

      idiskbhflag = 0

      open(file='denspsibulge.dat',unit=50,status='old')

      read(50,*) npsi,minlog
      do i=1,npsi
         read(50,*) tablepsib(i),denspsibulge(i)
      enddo

      close(50)
      return
      end
