      SUBROUTINE HETVAL(CMNAME,TEMP,TIME,DTIME,STATEV,FLUX,PREDEF,
     1                  DPRED)
C
C    INCLUDE 'ABA_PARAM.INC'
      implicit real(8) (a-h,o-z)
C
      CHARACTER*80 CMNAME
C
      DIMENSION TEMP(2),STATEV(*),PREDEF(*),TIME(2),FLUX(2),
     1 DPRED(*)

      TEMP(1) = 10.0
      TEMP(2) = 10.0

      RETURN
      END
