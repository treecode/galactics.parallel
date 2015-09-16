      function totdens(r,z)
      common /flags/ idiskflag, ibulgeflag, ihaloflag, ibhflag
      common /blackhole/ bhmass
      totdens = 0
      psi=pot(r,z)
      if( idiskflag .eq. 1 ) then
         totdens = totdens + diskdens(r,z,psi)
      endif      
      if( ihaloflag .eq. 1 ) then
         totdens = totdens + halodenspsi(psi)
      endif
      if( ibulgeflag .eq. 1 ) then
         totdens = totdens + bulgedenspsi(psi)
      endif
      return
      end
