*DECK MCRSX2
      SUBROUTINE MCRSX2(IPMPO,HEDIT,RECNAM,NREA,NGRP,NMGF,NL,ISO,
     1 NOMREA,NOMISO,DEN,FACT,WEIGHT,SPH,FLUXS,IREAB,IREAF,LPURE,
     2 IGYELD,LXS,XS,SIGS,SS2D,TAUXFI,TAUXGF)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover the cross sections of an elementary calculation and single
* mixture in an MPO file and perform multiparameter interpolation.
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
* IPMPO   pointer to the MPO file.
* HEDIT   name of output group for a (multigroup mesh, output geometry)
*         couple (generally equal to 'output_0').
* RECNAM  character identification of calculation.
* NREA    number of reactions in the MPO file.
* NGRP    number of energy groups.
* NMGF    number of macrogroups for the fission yields.
* NL      maximum Legendre order (NL=1 is for isotropic scattering).
* ISO     isotope index.
* NOMREA  names of reactions in the MPO file.
* NOMISO  name of isotope ISO.
* DEN     number density of isotope.
* FACT    number density ratio for the isotope.
* WEIGHT  interpolation weight.
* SPH     SPH factors.
* FLUXS   averaged flux.
* IREAB   position of 'Absorption' reaction in NOMREA array.
* IREAF   position of 'NuFission' reaction in NOMREA array.
* LPURE   =.true. if the interpolation is a pure linear interpolation 
*         with TERP factors.
* IGYELD  yield macrogroup limits.
*
*Parameters: input/output
* LXS     existence flag of each reaction.
* XS      interpolated cross sections per reaction
* SIGS    interpolated scattering cross sections
* SS2D    interpolated scattering matrix
* TAUXFI  interpolated fission rate
* TAUXGF  interpolated fission rate in macrogroups
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE hdf5_wrap
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMPO
      CHARACTER(LEN=12) HEDIT
      CHARACTER(LEN=80) RECNAM
      INTEGER NREA,NGRP,NMGF,NL,ISO,IREAB,IREAF,IGYELD(NMGF)
      REAL DEN,FACT,WEIGHT,SPH(NGRP),FLUXS(NGRP),SS2D(NGRP,NGRP,NL),
     1 SIGS(NGRP,NL),XS(NGRP,NREA),TAUXFI,TAUXGF(NMGF)
      LOGICAL LXS(NREA),LPURE
      CHARACTER NOMREA(NREA)*24,NOMISO*24
*----
*  LOCAL VARIABLES
*----
      INTEGER IREA,IOF,IL,IGR,JGR,IGRC,IGRDEB,IGRFIN,ADDRZX,ADDRZI,
     1 IPROF,ISOM,JOFS,NISO,NL1,NL2,RANK,TYPE,NBYTE,DIMSR(5)
      REAL FLOTT,TAUXF,ZIL,B2
      CHARACTER RECNAM2*80
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDATAP,FAG,ADR,ADDRISO
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: ADDRXS
      REAL, ALLOCATABLE, DIMENSION(:) :: RDATAX,DIFF
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SIGSB,XSB
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SS2DB
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(SIGSB(NGRP,NL),SS2DB(NGRP,NGRP,NL),XSB(NGRP,NREA),
     1 FAG(NGRP),ADR(NGRP))
*----
*  FIND THE ISOTOPE INDEX IN ADDRXS
*----
      WRITE(RECNAM2,'(8H/output/,A,6H/info/)') TRIM(HEDIT)
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM2)//"NISO",NISO)
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM2)//"ADDRXS",ADDRXS)
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM2)//"ADDRISO",ADDRISO)
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM2)//"TRANSPROFILE",IDATAP)
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"ADDRZI",ADDRZI)
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"ADDRZX",ADDRZX)
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"CROSSECTION",RDATAX)
      ISOM=ISO-ADDRISO(ADDRZI+1)
      IF((ISOM.LE.0).OR.(ISOM.GT.NISO)) CALL XABORT('MCRSX2: ADDRXS OV'
     1 //'ERFLOW.')
      NL1=ADDRXS(NREA-1,ISOM,ADDRZX+1)
      NL2=ADDRXS(NREA,ISOM,ADDRZX+1)
      IF((NL1.GT.NL).OR.(NL2.GT.NL)) CALL XABORT('MCRSX2: NL OVERFLOW.')
*----
*  LOOP OVER REACTIONS
*----
      SIGSB(:NGRP,:NL)=0.0
      SS2DB(:NGRP,:NGRP,:NL)=0.0
      XSB(:NGRP,:NREA)=0.0
      DO IREA=1,NREA-2
        IOF=ADDRXS(IREA,ISOM,ADDRZX+1)
        IF(IOF.LT.0) CYCLE
        LXS(IREA)=.TRUE.
        IF(NOMREA(IREA).EQ.'Diffusion') THEN
          DO IL=1,NL1
            DO IGR=1,NGRP
              FLOTT=RDATAX(IOF+(IL-1)*NGRP+IGR)
              SIGSB(IGR,IL)=SIGSB(IGR,IL)+FLOTT
            ENDDO
          ENDDO
        ELSE IF(NOMREA(IREA).EQ.'Scattering') THEN
          IPROF=ADDRXS(NREA+1,ISOM,ADDRZX+1)
          DO IGR=1,NGRP
            FAG(IGR)=IDATAP(IPROF+IGR)+1
            ADR(IGR)=IDATAP(IPROF+NGRP+IGR)
          ENDDO
          ADR(NGRP+1)=IDATAP(IPROF+1+2*NGRP)
          JOFS=0
          DO IL=1,NL2
            ZIL=REAL(2*IL-1)
            DO IGR=1,NGRP
              DO JGR=FAG(IGR),FAG(IGR)+(ADR(IGR+1)-ADR(IGR))-1
                IF(JGR.GT.NGRP) CALL XABORT('MCRSX2: SS2D OVERFLOW.')
                FLOTT=RDATAX(IOF+JOFS+1)/ZIL
                SS2DB(JGR,IGR,IL)=SS2DB(JGR,IGR,IL)+FLOTT ! JGR <-- IGR
                JOFS=JOFS+1
              ENDDO
            ENDDO
          ENDDO
        ELSE
          XSB(:NGRP,IREA)=RDATAX(IOF+1:IOF+NGRP)
        ENDIF
      ENDDO ! end of loop over reactions
      DEALLOCATE(IDATAP,RDATAX,ADDRISO,ADDRXS)
      LXS(NREA-1)=.TRUE.
*----
*  RECOVER DIFFUSION COEFFICIENT INFORMATION
*----
      IF(NOMISO.EQ.'TotalResidual_mix') THEN
        IF(hdf5_group_exists(IPMPO,TRIM(RECNAM)//"leakage")) THEN
          CALL hdf5_info(IPMPO,TRIM(RECNAM)//"leakage/DIFFCOEF",RANK,
     1    TYPE,NBYTE,DIMSR)
          IF(TYPE.NE.99) THEN
            LXS(NREA)=.TRUE.
            CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"leakage/DIFFCOEF",
     1      DIFF)
            XSB(:NGRP,NREA)=DIFF(:NGRP)*DEN
            DEALLOCATE(DIFF)
            GO TO 10
          ENDIF
          CALL hdf5_info(IPMPO,TRIM(RECNAM)//"leakage/DB2",RANK,TYPE,
     1    NBYTE,DIMSR)
          IF(TYPE.NE.99) THEN
            LXS(NREA)=.TRUE.
            CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"leakage/BUCKLING",
     1      B2)
            CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"leakage/DB2",DIFF)
            DO IGR=1,NGRP
              XSB(IGR,NREA)=DIFF(IGR)*DEN/B2
            ENDDO
            DEALLOCATE(DIFF)
          ENDIF
        ENDIF
      ENDIF
*----
*  COMPUTE FISSION RATE FOR AN ELEMENTARY CALCULATION
*----
   10 TAUXF=0.0
      TAUXGF(:NMGF)=0.0
      IF(IREAF.GT.0) THEN
        DO IGR=1,NGRP
          TAUXF=TAUXF+XSB(IGR,IREAF)*FLUXS(IGR)
        ENDDO
        TAUXFI=TAUXFI+WEIGHT*FACT*TAUXF
        IGRFIN=0
        DO IGRC=1,NMGF
          IGRDEB=IGRFIN+1
          IGRFIN=IGYELD(IGRC)
          DO IGR=IGRDEB,IGRFIN
            TAUXGF(IGRC)=TAUXGF(IGRC)+XSB(IGR,IREAF)*FLUXS(IGR)
          ENDDO
          TAUXGF(:NMGF)=WEIGHT*FACT*TAUXGF(:NMGF)
        ENDDO
      ENDIF
*----
*  WEIGHT MICROSCOPIC CROSS SECTION DATA IN AN INTERPOLATED MICROLIB
*----
      DO IGR=1,NGRP
        DO IREA=1,NREA
          IF(.NOT.LXS(IREA)) CYCLE
          IF(NOMREA(IREA).EQ.'Total') THEN
            XS(IGR,IREA)=XS(IGR,IREA)+FACT*SPH(IGR)*WEIGHT*
     1      (XSB(IGR,IREAB)+SIGSB(IGR,1))
          ELSE IF(LPURE.AND.NOMREA(IREA).EQ.'FissionSpectrum') THEN
            XS(IGR,IREA)=XS(IGR,IREA)+WEIGHT*XSB(IGR,IREA)
          ELSE IF(NOMREA(IREA).EQ.'FissionSpectrum') THEN
            IF(IREAF.EQ.0) CALL XABORT('MCRSX2: IREAF=0.')
            XS(IGR,IREA)=XS(IGR,IREA)+WEIGHT*FACT*TAUXF*XSB(IGR,IREA)
          ELSE
            XS(IGR,IREA)=XS(IGR,IREA)+FACT*SPH(IGR)*WEIGHT*XSB(IGR,IREA)
          ENDIF
        ENDDO
        DO IL=1,NL
          IF(MOD(IL,2).EQ.1) THEN
            SIGS(IGR,IL)=SIGS(IGR,IL)+FACT*SPH(IGR)*WEIGHT*SIGSB(IGR,IL)
          ELSE
            DO JGR=1,NGRP
              SIGS(IGR,IL)=SIGS(IGR,IL)+FACT*WEIGHT*SS2DB(JGR,IGR,IL)
     1        /SPH(JGR)
            ENDDO
          ENDIF
        ENDDO
        DO JGR=1,NGRP
          DO IL=1,NL
            IF(MOD(IL,2).EQ.1) THEN
              SS2D(IGR,JGR,IL)=SS2D(IGR,JGR,IL)+FACT*SPH(JGR)*WEIGHT*
     1        SS2DB(IGR,JGR,IL)
            ELSE
              SS2D(IGR,JGR,IL)=SS2D(IGR,JGR,IL)+FACT*WEIGHT*
     1        SS2DB(IGR,JGR,IL)/SPH(IGR)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(ADR,FAG,XSB,SS2DB,SIGSB)
      RETURN
      END
