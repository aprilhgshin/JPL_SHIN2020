      program hello
      REAL*4 actual, tarray(2)
      EXTERNAL ETIME
      REAL*4 ETIME
      actual = ETIME( tarray )
      print *, tarray
      end
