*DECK HEADRV
      SUBROUTINE HEADRV(IPDEP,NPART,IPMAC,NBMIX,NGRP,ZNORM,IMPX,ESUM,
     1 CSUM,IBC,RHO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute energy and charge deposition from many particles.
*
*Copyright:
* Copyright (C) 2020 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPDEP   L_DEPOSITION pointer to the deposition information object.
* NPART   number of particles contributing to energy and charge
*         deposition.
* IPMAC   L_MACROLIB pointers to the extended macrolibs.
* NBMIX   number of material mixtures.
* NGRP    total number of energy groups.
* ZNORM   flux normalization factor.
* IMPX    print parameter.
*
*Parameters: output
* ESUM    total energy deposition (MeV/cc).
* CSUM    total charge deposition (electron/cc). 
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NPART,NBMIX,NGRP,IMPX
      TYPE(C_PTR) IPDEP,IPMAC(NPART)
      REAL RHO(NBMIX)
      DOUBLE PRECISION ZNORM,ESUM,CSUM
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPMAC,KPMAC
      PARAMETER(NSTATE=40,IOUT=6)
      INTEGER ISTATE(NSTATE)
      CHARACTER HSMG*131,TEXT1*1
      DOUBLE PRECISION VTOT
      LOGICAL LCHARG
      REAL FLUXC(NBMIX)
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: VOL,SGD,FLIN
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: EDEPOT,CDEPOT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: EDEPO,CDEPO
      CHARACTER(LEN=1), ALLOCATABLE, DIMENSION(:) :: TEXT1V
      CHARACTER(LEN=6), ALLOCATABLE, DIMENSION(:) :: SNAME
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: FUNA8
*----
*  MEMORY ALLOCATION
*----
      ALLOCATE(EDEPO(NBMIX,NPART),CDEPO(NBMIX,NPART),EDEPOT(NBMIX),
     1 CDEPOT(NBMIX),VOL(NBMIX),SGD(NBMIX),FLIN(NBMIX))
      ALLOCATE(FUNA8(3*NPART),SNAME(3*NPART))
*----
*  RECOVER ENERGY AND CHARGE DEPOSITION
*----
      CALL LCMLEN(IPDEP,'EDEPOS',LENGT,ITYLCM)
      IF(LENGT.NE.0) THEN
        IF(LENGT.NE.NBMIX*NPART) CALL XABORT('HEADRV: INVALID EDEPOS R'
     1  //'ECORD LENGTH.')
        CALL LCMGET(IPDEP,'EDEPOS',EDEPO)
        CALL LCMGET(IPDEP,'FLUX-NORM',ZNORM)
      ELSE
        CALL XDDSET(EDEPO,NBMIX*NPART,0.0D0)
      ENDIF
      CALL LCMLEN(IPDEP,'CDEPOS',LENGT,ITYLCM)
      IF(LENGT.NE.0) THEN
        IF(LENGT.NE.NBMIX*NPART) CALL XABORT('HEADRV: INVALID CDEPOS R'
     1  //'ECORD LENGTH.')
        CALL LCMGET(IPDEP,'CDEPOS',CDEPO)
      ELSE
        CALL XDDSET(CDEPO,NBMIX*NPART,0.0D0)
      ENDIF
      CALL XDDSET(EDEPOT,NBMIX,0.0D0)
      CALL XDDSET(CDEPOT,NBMIX,0.0D0)
      DO I=1,NPART
        CALL LCMLEN(IPMAC(I),'FLUXC',IBC2,ITYLCM)
        IF(IBC.EQ.1.AND.IBC2.NE.0) THEN
        CALL LCMGET(IPMAC(I),'FLUXC',FLUXC)
        CALL LCMGET(IPMAC(I),'ECUTOFF',ECUTOFF)
        ENDIF
        CALL LCMGET(IPMAC(I),'VOLUME',VOL)
        CALL LCMGTC(IPMAC(I),'PARTICLE',1,1,TEXT1)
        SNAME(I)=TEXT1
        FUNA8(I)='ENERGDEP'
        FUNA8(NPART+I+1)='CHARGDEP'
        SNAME(NPART+I+1)=SNAME(I)
        JPMAC=LCMGID(IPMAC(I),'GROUP')
        DO IGR=1,NGRP
          KPMAC=LCMGIL(JPMAC,IGR)
          CALL LCMLEN(KPMAC,'H-FACTOR',LENGT,ITYLCM)
          IF(LENGT.EQ.0) THEN
            WRITE(HSMG,'(42HHEADRV: NO H-FACTOR FOUND IN MACROLIB NUMB,
     1      2HER,I3,1H.)') I
            CALL XABORT(HSMG)
          ENDIF
          CALL LCMGET(KPMAC,'FLUX-INTG',FLIN)
          CALL LCMGET(KPMAC,'H-FACTOR',SGD)
          DO IBM=1,NBMIX
            EDEPO(IBM,I)=EDEPO(IBM,I)+FLIN(IBM)*SGD(IBM)*ZNORM/
     1      (VOL(IBM)*RHO(IBM))
            IF(IBC.EQ.1.AND.IBC2.NE.0) THEN
              EDEPO(IBM,I)=EDEPO(IBM,I)+ECUTOFF*FLUXC(IBM)*ZNORM/
     1        RHO(IBM)
            ENDIF
          ENDDO
          CALL LCMLEN(KPMAC,'C-FACTOR',LENGT,ITYLCM)
          IF(LENGT.GT.0) THEN
            LCHARG=.TRUE.
            CALL LCMGET(KPMAC,'C-FACTOR',SGD)
            DO IBM=1,NBMIX
              CDEPO(IBM,I)=CDEPO(IBM,I)+FLIN(IBM)*SGD(IBM)*ZNORM/
     1        (VOL(IBM)*RHO(IBM))
              IF(IBC.EQ.1.AND.IBC2.NE.0) THEN
                CDEPO(IBM,I)=CDEPO(IBM,I)-FLUXC(IBM)*ZNORM/RHO(IBM)
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      FUNA8(NPART+1)='TOTENDEP'
      FUNA8(2*NPART+2)='TOTCHDEP'
      SNAME(NPART+1)='EDEPO'
      SNAME(2*NPART+2)='CDEPO'
      VTOT=0.0D0
      ESUM=0.0D0
      CSUM=0.0D0
      DO IBM=1,NBMIX
        DO I=1,NPART
          EDEPOT(IBM)=EDEPOT(IBM)+EDEPO(IBM,I)
          CDEPOT(IBM)=CDEPOT(IBM)+CDEPO(IBM,I)
        ENDDO
        VTOT=VTOT+VOL(IBM)
        ESUM=ESUM+EDEPOT(IBM)*VOL(IBM)
        CSUM=CSUM+CDEPOT(IBM)*VOL(IBM)
      ENDDO
      ESUM=ESUM/VTOT
      CSUM=CSUM/VTOT
*----
*  PRINT ENERGY AND CHARGE DEPOSITION
*----
      IF(IMPX.GT.0) THEN
        WRITE(IOUT,1001) '   VOLUME  ',(FUNA8(J),'_',SNAME(J),J=1,
     1  2*NPART+2)
        DO IBM=1,NBMIX
          WRITE(IOUT,1002) VOL(IBM),(EDEPO(IBM,I),I=1,NPART),
     1    EDEPOT(IBM),(CDEPO(IBM,I),I=1,NPART),CDEPOT(IBM)
        ENDDO
        WRITE(IOUT,'(/14H TOTAL VOLUME:,14X,1P,E12.4)') VTOT
        WRITE(IOUT,'(25H TOTAL ENERGY DEPOSITION:,3X,1P,E12.4)') ESUM
        WRITE(IOUT,'(25H TOTAL CHARGE DEPOSITION:,3X,1P,E12.4)') CSUM
      ENDIF
*----
*  SAVE ENERGY AND CHARGE DEPOSITION
*----
      CALL LCMPUT(IPDEP,'VOLUME',NBMIX,2,VOL)
      CALL LCMPUT(IPDEP,'EDEPOS',NBMIX*NPART,4,EDEPO)
      CALL LCMPUT(IPDEP,'EDEPOS_TOT',NBMIX,4,EDEPOT)
      IF(LCHARG) THEN
        CALL LCMPUT(IPDEP,'CDEPOS',NBMIX*NPART,4,CDEPO)
        CALL LCMPUT(IPDEP,'CDEPOS_TOT',NBMIX,4,CDEPOT)
      ENDIF
      CALL LCMPUT(IPDEP,'FLUX-NORM',1,4,ZNORM)
*----
*  PROCESS STATE-VECTOR
*----
      CALL LCMLEN(IPDEP,'STATE-VECTOR',LENGT,ITYLCM)
      IF(LENGT.NE.0) THEN
        CALL LCMGET(IPDEP,'STATE-VECTOR',ISTATE)
        ALLOCATE(TEXT1V(NPART))
        CALL LCMGTC(IPDEP,'PARTICLE-NAM',1,NPART,TEXT1V)
        DO I=1,NPART
          IF(TEXT1V(I).NE.SNAME(I)(:1)) THEN
            WRITE(HSMG,'(22HHEADRV: PARTICLE NAMES,2A2,
     1      16H ARE INCOHERENT.)') TEXT1V(I),SNAME(I)(:1)
            CALL XABORT(HSMG)
          ENDIF
        ENDDO
        DEALLOCATE(TEXT1V)
      ELSE
        CALL XDISET(ISTATE,NSTATE,0)
        ISTATE(1)=NBMIX
        ISTATE(2)=NPART
        IF(LCHARG) ISTATE(3)=1
        ALLOCATE(TEXT1V(NPART))
        DO I=1,NPART
          TEXT1V(I)=SNAME(I)(:1)
        ENDDO
        CALL LCMPTC(IPDEP,'PARTICLE-NAM',1,NPART,TEXT1V)
        DEALLOCATE(TEXT1V)
      ENDIF
      ISTATE(4)=ISTATE(4)+1
      CALL LCMPUT(IPDEP,'STATE-VECTOR',NSTATE,1,ISTATE)
*----
*  MEMORY DEALLOCATION
*----
      DEALLOCATE(SNAME,FUNA8)
      DEALLOCATE(FLIN,SGD,VOL,CDEPOT,EDEPOT,CDEPO,EDEPO)
      RETURN
*
 1001 FORMAT(/1X,A11,21(1X,A8,A1,A6))
 1002 FORMAT(1X,1P,E11.4,21E16.4)
      END
