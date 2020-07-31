      program hello
      REAL*4 actual, tarray(2)
      actual = -999.
      call ETIME( tarray, actual )
      if ( actual.ge.0. ) then
        print *, 0, tarray, actual
      else
        print *, 1, tarray, actual
      endif
      end
