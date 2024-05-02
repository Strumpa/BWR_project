*DECK EDIWCU
      SUBROUTINE EDIWCU(IPFLUX,IPRINT,NGROUP,NUN,NREGIO,NDIM,NLIN,
     > NFUNL,NGCOND,NMERGE,KEYANI,VOLUME,IGCOND,IMERGE,COUWP1)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Homogenize the currents based on spherical harmonic moments of the
* flux.
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
* IPFLUX  pointer to the flux LCM object.
* IPRINT  print parameter.
* NGROUP  number of energy groups.
* NUN     number of unknowns in flux array.
* NREGIO  number of regions.
* NDIM    number of dimensions.
* NLIN    number of polynomial components in flux.
* NFUNL   number of spherical harmonic components in flux.
* NGCOND  number of merged energy groups.
* NMERGE  number of merged regions.
* KEYANI  position of spherical harmonic components in unknown vector.
* VOLUME  volumes.
* IGCOND  limit condensed groups.
* IMERGE  region merging matrix.
*
*Parameters: input/output
* COUWP1  homogenized and condensed currents.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPFLUX
      INTEGER    IPRINT,NREGIO,NGROUP,NUN,NDIM,NLIN,NFUNL,NGCOND,NMERGE,
     >           KEYANI(NREGIO,NLIN,NFUNL),IGCOND(NGCOND),
     >           IMERGE(NREGIO)
      REAL       VOLUME(NREGIO),COUWP1(NMERGE,NGCOND,NDIM)
*----
*  LOCAL VARIABLES
*----
      PARAMETER  (IUNOUT=6)
      TYPE(C_PTR) JPFLUX
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) ::  WORKF
      REAL, ALLOCATABLE, DIMENSION(:,:,:) ::  FLUXES
*----
*  INITIALIZATION
*----
      IF(NFUNL.EQ.1) CALL XABORT('EDIWCU: ANIS.GE.2 EXPECTED IN TRACKI'
     > //'NG.')
      JPFLUX=LCMGID(IPFLUX,'FLUX')
      ALLOCATE(WORKF(NUN),FLUXES(NREGIO,NGROUP,NDIM))
      COUWP1(:NMERGE,:NGCOND,:NDIM)=0.0
*----
*  PROCESS TRIVIAL 1D CASE
*----
      IF(NDIM.EQ.1) THEN
        DO IGR=1,NGROUP
          CALL LCMGDL(JPFLUX,IGR,WORKF)
          DO IREG=1,NREGIO
            FLUXES(IREG,IGR,1)=WORKF(KEYANI(IREG,1,2))
          ENDDO
        ENDDO
      ELSE IF(NDIM.GT.1) THEN
*----
*  PROCESS 2D AND 3D CASES
*----
        IL=1
        IOF0=1
        DO IGR=1,NGROUP
          CALL LCMGDL(JPFLUX,IGR,WORKF)
          DO IREG=1,NREGIO
            IOF=IOF0
            DO IM=-IL,IL
              IF((NDIM.EQ.2).AND.(MOD(IL+IM,2).EQ.1)) CYCLE
              IOF=IOF+1
              IF(IOF.GT.NFUNL) CALL XABORT('EDIWCU: KEYANI OVERFLOW.')
              IF(IM.EQ.-1) THEN
                FLUXES(IREG,IGR,2)=WORKF(KEYANI(IREG,1,IOF))
              ELSE IF(IM.EQ.0) THEN
                FLUXES(IREG,IGR,3)=WORKF(KEYANI(IREG,1,IOF))
              ELSE IF(IM.EQ.1) THEN
                FLUXES(IREG,IGR,1)=WORKF(KEYANI(IREG,1,IOF))
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF
*----
*  CONDENSATION AND HOMOGENIZATION OF SPHERICAL HARMONIC MOMENTS
*----
      IGRFIN=0
      DO IGRC=1,NGCOND
        IGRDEB=IGRFIN+1
        IGRFIN=IGCOND(IGRC)
        DO IGR=IGRDEB,IGRFIN
          DO IREG=1,NREGIO
            IRA=IMERGE(IREG)
            IF(IRA.EQ.0) CYCLE
            DVOL=VOLUME(IREG)
            DO ID=1,NDIM
              COUWP1(IRA,IGRC,ID)=COUWP1(IRA,IGRC,ID)+
     >                            FLUXES(IREG,IGR,ID)*DVOL
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(FLUXES,WORKF)
*----
*  PRINTOUTS
*----
      IF(IPRINT.GT.0) THEN
        WRITE(6,'(/42H EDIWCU: INCLUDE CURRENTS IN THE MACROLIB.)')
        DO IDIM=1,NDIM
          DO IGR=1,NGCOND
            IF(IDIM.EQ.1) WRITE(IUNOUT,6010) IGR,'X'
            IF(IDIM.EQ.2) WRITE(IUNOUT,6010) IGR,'Y'
            IF(IDIM.EQ.3) WRITE(IUNOUT,6010) IGR,'Z'
            WRITE(IUNOUT,6012) (COUWP1(IKK,IGR,IDIM),IKK=1,NMERGE)
          ENDDO
        ENDDO
      ENDIF
      RETURN
*
 6010 FORMAT(/' G R O U P   :',I4/' REGION INTEGRATED ',A1,'-CURRENT')
 6012 FORMAT(1P,7(3X,E15.7))
      END
