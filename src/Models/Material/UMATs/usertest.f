      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
C    INCLUDE 'ABA_PARAM.INC'
      implicit real(8) (a-h,o-z)
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),
     2 DDSDDT(NTENS),DRPLDE(NTENS),
     3 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     4 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     5 JSTEP(4)
C
      I = 0.0
      DO K1=1,NTENS
            DO K2=1,NTENS
                  DDSDDE(K1,K2) = I*PROPS(1)
                  I=I+1
            END DO
      END DO
      DO K1=1,NTENS
            STRESS(K1) = STRESS(K1)*PROPS(1)
            STRAN(K1) = STRAN(K1)*PROPS(1)
            DSTRAN(K1) = DSTRAN(K1)*PROPS(1)
            DDSDDT(K1) = DDSDDT(K1)*PROPS(1)
            DRPLDE(K1) = DRPLDE(K1)*PROPS(1)
      END DO
      DO K1=1,NSTATV
            STATEV(K1) = STATEV(K1)*PROPS(1)
      END DO
      DO K1=1,3
            COORDS(K1) = COORDS(K1)*PROPS(1)
      END DO
      DO K1=1,3
            DO K2=1,3
                  DROT(K1,K2) = DROT(K1,K2)*PROPS(1)
                  DFGRD0(K1,K2) = DFGRD0(K1,K2)*PROPS(1)
                  DFGRD1(K1,K2) = DFGRD1(K1,K2)*PROPS(1)
            END DO
      END DO
      RETURN
      END
