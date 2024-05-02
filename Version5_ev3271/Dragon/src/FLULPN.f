*DECK FLULPN
      SUBROUTINE FLULPN(IPMACR,NUNKNO,OPTION,TYPE,NGRP,NREG,NMAT,NIFIS,
     1 VOL,MATCOD,NMERG,IMERG,KEYFLX,FLUX,IMPX,DIFHET,AKEFF,B2,OLDBIL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of heterogeneous leakage coefficients in non-fundamental
* mode condition using the Todorova approximation.
*
*Copyright:
* Copyright (C) 2022 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPMACR  pointer to the macrolib LCM object (L_MACROLIB signature).
* NUNKNO  number of flux/current unknowns.
* OPTION  type of leakage coefficients; can be 'P0' (P-0), 'P1' (P-1),
*         'P0TR' (P-0 with transport correction).
* TYPE    type of buckling iteration.
*         Can be 'DIFF' (do a P0 calculation of DIFHET and exit);
*                'K' (do a P-n calculation with keff search).
* NGRP    number of groups.
* NREG    number of volumes.
* NMAT    number of material mixtures.
* NIFIS   maximum number of fission spectrum assigned to a mixture.
* VOL     volumes.
* MATCOD  mixture number of each volume.
* NMERG   number of leakage zones.
* IMERG   leakage zone index in each material mixture zone.
* KEYFLX  position of each flux in the unknown vector.
* FLUX    direct unknown vector.
* IMPX    print flag.
* B2      buckling.
* OLDBIL  previous norm of the flux.
*
*Parameters: output
* DIFHET  diffusion coefficients.
* AKEFF   effective multiplication factor.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER*4 OPTION,TYPE
      TYPE(C_PTR) IPMACR
      INTEGER NUNKNO,NGRP,NREG,NMAT,NIFIS,MATCOD(NREG),NMERG,
     1 IMERG(NMAT),KEYFLX(NREG),IMPX
      REAL VOL(NREG),FLUX(NUNKNO,NGRP),DIFHET(NMERG,NGRP),B2
      DOUBLE PRECISION AKEFF,OLDBIL
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPMACR,KPMACR
      CHARACTER HSMG*131
      DOUBLE PRECISION DDELN1,DDELD1
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ,IPOS
      REAL, ALLOCATABLE, DIMENSION(:) :: WORK
      REAL, ALLOCATABLE, DIMENSION(:,:) :: ST,FLXIN
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SFNU,SCAT1
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: STOD
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(FLXIN(NMAT,NGRP),ST(NMAT,NGRP),SFNU(NMAT,NIFIS,NGRP),
     1 SCAT1(NMAT,NGRP,NGRP))
*----
*  INITIALIZATION
*----
      IF((IMPX.GT.0).AND.(TYPE.EQ.'DIFF')) THEN
         WRITE (6,100)
      ELSE IF(IMPX.GT.0) THEN
         WRITE (6,110) OPTION,TYPE,B2
      ENDIF
      ST(:NMAT,:NGRP)=0.0
      SFNU(:NMAT,:NIFIS,:NGRP)=0.0
      SCAT1(:NMAT,:NGRP,:NGRP)=0.0
*----
*  RECOVER MACROSCOPIC CROSS SECTIONS
*----
      ALLOCATE(IJJ(NMAT),NJJ(NMAT),IPOS(NMAT),WORK(NMAT*NGRP))
      JPMACR=LCMGID(IPMACR,'GROUP')
      DO IGR=1,NGRP
        KPMACR=LCMGIL(JPMACR,IGR)
        CALL LCMGET(KPMACR,'NTOT0',ST(1,IGR))
        IF(NIFIS.GT.0) CALL LCMGET(KPMACR,'NUSIGF',SFNU(1,1,IGR))
        CALL LCMLEN(KPMACR,'SCAT01',ILONG,ITYLCM)
        IF((ILONG.NE.0).AND.(OPTION.NE.'P0')) THEN
          CALL LCMGET(KPMACR,'IJJS01',IJJ)
          CALL LCMGET(KPMACR,'NJJS01',NJJ)
          CALL LCMGET(KPMACR,'IPOS01',IPOS)
          CALL LCMGET(KPMACR,'SCAT01',WORK)
          DO IBM=1,NMAT
            IPO=IPOS(IBM)
            DO JGR=IJJ(IBM),IJJ(IBM)-NJJ(IBM)+1,-1
              SCAT1(IBM,IGR,JGR)=WORK(IPO) ! IGR <-- JGR
              IPO=IPO+1
            ENDDO
          ENDDO
        ENDIF
      ENDDO
      DEALLOCATE(WORK,IPOS,NJJ,IJJ)
*----
*  RECOVER INTEGRATED FLUX
*----
      FLXIN(:NMAT,:NGRP)=0.0
      DO IGR=1,NGRP
        DO IBM=1,NMAT
          DO I=1,NREG
            IND=KEYFLX(I)
            IF((MATCOD(I).EQ.IBM).AND.(IND.GT.0)) THEN
              FLXIN(IBM,IGR)=FLXIN(IBM,IGR)+FLUX(IND,IGR)*VOL(I)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      IF((OPTION.EQ.'LKRD').OR.(OPTION.EQ.'RHS')) GO TO 10
*----
*  MAIN LOOP OVER LEAKAGE ZONES
*----
      DO INM=1,NMERG
        IF(OPTION.EQ.'P0') THEN
*         P0 approximation
          DO IGR=1,NGRP
            DDELN1=0.D0
            DDELD1=0.D0
            DO IBM=1,NMAT
              IF(IMERG(IBM).EQ.INM) THEN
                DDELN1=DDELN1+ST(IBM,IGR)*FLXIN(IBM,IGR)
                DDELD1=DDELD1+FLXIN(IBM,IGR)
              ENDIF
            ENDDO
            STR=REAL(DDELN1/DDELD1)
            DIFHET(INM,IGR)=1.0/(3.0*STR)
          ENDDO
        ELSE IF((OPTION.EQ.'P0TR').OR.(TYPE.EQ.'DIFF')) THEN
*         Outscatter approximation
          DO IGR=1,NGRP
            DDELN1=0.D0
            DDELD1=0.D0
            DO IBM=1,NMAT
              IF(IMERG(IBM).EQ.INM) THEN
                DDELN1=DDELN1+ST(IBM,IGR)*FLXIN(IBM,IGR)
                DO JGR=1,NGRP
                  DDELN1=DDELN1-SCAT1(IBM,JGR,IGR)*FLXIN(IBM,IGR)
                ENDDO
                DDELD1=DDELD1+FLXIN(IBM,IGR)
              ENDIF
            ENDDO
            STR=REAL(DDELN1/DDELD1)
            DIFHET(INM,IGR)=1.0/(3.0*STR)
          ENDDO
        ELSE IF(OPTION.EQ.'P1') THEN
*         Inscatter approximation
          ALLOCATE(STOD(NGRP,NGRP+1))
          STOD(:NGRP,:NGRP+1)=0.0D0
          DO IGR=1,NGRP
            DO IBM=1,NMAT
              IF(IMERG(IBM).EQ.INM) THEN
                STOD(IGR,IGR)=STOD(IGR,IGR)+ST(IBM,IGR)*FLXIN(IBM,IGR)
                DO JGR=1,NGRP
                  STOD(IGR,JGR)=STOD(IGR,JGR)-SCAT1(IBM,IGR,JGR)*
     1            FLXIN(IBM,JGR)
                ENDDO
                STOD(IGR,NGRP+1)=STOD(IGR,NGRP+1)+FLXIN(IBM,IGR)/3.0D0
              ENDIF
            ENDDO
          ENDDO
          CALL ALSBD(NGRP,1,STOD,IER,NGRP)
          IF(IER.NE.0) CALL XABORT('FLULPN: SINGULAR MATRIX.')
          DO IGR=1,NGRP
            DIFHET(INM,IGR)=REAL(STOD(IGR,NGRP+1))
          ENDDO
          DEALLOCATE(STOD)
        ELSE
          WRITE(HSMG,'(15HFLULPN: OPTION ,A,23H IS INVALID WITH TODORO,
     1    17HVA APPROXIMATION.)') OPTION
          CALL XABORT(HSMG)
        ENDIF
      ENDDO
      IF(TYPE.EQ.'DIFF') GO TO 20
*----
*  COMPUTE THE EFFECTIVE MULTIPLICATION FACTOR.
*----
   10 IF(TYPE.NE.'K') CALL XABORT('FLULPN: TYPE K EXPECTED.')
      PROD=0.0D0
      DO IGR=1,NGRP
        DO NF=1,NIFIS
          DO IBM=1,NMAT
            PROD=PROD+SFNU(IBM,NF,IGR)*FLXIN(IBM,IGR)
          ENDDO
        ENDDO
      ENDDO
      AKEFF=AKEFF*PROD/OLDBIL
      OLDBIL=PROD
      IF(IMPX.GT.0) WRITE (6,120) B2,AKEFF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
   20 DEALLOCATE(SCAT1,SFNU,ST,FLXIN)
      RETURN
*
  100 FORMAT(/54H FLULPN: OUTSCATTER DIFFUSION COEFFICIENT CALCULATION.)
  110 FORMAT(/21H FLULPN: SOLUTION OF ,A4,21H EQUATIONS WITH TYPE ,A4/
     1 9X,17HINITIAL BUCKLING=,1P,E13.5)
  120 FORMAT(/18H FLULPN: BUCKLING=,1P,E13.5,15H K-EFFECTIVE  =,E13.5)
      END
