      function bulgedens(r,z)
      psi = pot(r,z)
      bulgedens = bulgedenspsi(psi)
      write(0,*) psi,bulgedens,r,z
      stop
      return
      end
