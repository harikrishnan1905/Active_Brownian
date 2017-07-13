      DOUBLE PRECISION FUNCTION gasdev(iseed)
      IMPLICIT NONE
c
c     normally distributed with variance sigma and mean zero
c
      DOUBLE PRECISION r, v1, v2, fac, gset, ran2
      INTEGER iset, iseed
 
      SAVE gset, iset
      DATA iset/0/
 100  IF (iset.EQ.0) THEN
         v1 = 2.d0*ran2(iseed) - 1.d0
         v2 = 2.d0*ran2(iseed) - 1.d0
         r = v1**2 + v2**2
         IF (r.GE.1) GOTO 100
         fac = sqrt(-2.d0*log(r)/r)
         gset = v1*fac
         gasdev = v2*fac
         iset = 1
      ELSE
         gasdev = gset
         iset = 0
      END IF
      Return
      End