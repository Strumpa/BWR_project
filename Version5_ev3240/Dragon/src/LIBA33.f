*DECK LIBA33
      SUBROUTINE LIBA33(NG,NANI,TT,NT0,NPSN0,FGTD,TEMP,IAFAG,IFAGR,
     1 PSN0,SCAT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly and temperature interpolation of a transfer matrix stored
* in the APOLIB-3 format.
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
* NG      number of energy groups.
* NANI    anisotropy level. NANI=1 for isotropic scattering.
* TT      temperature of isotope.
* NT0     number of tabulated temperatures.
* NPSN0   size of vector PSN0.
* FGTD    first temperature-dependent group.
* TEMP    tabulated temperatures.
* IAFAG   address for the first arrival group XS
* IFAGR   first arrival group index.
* PSN0    input cross section data in APOLIB-3 compressed format.
*
*Parameters: output
* SCAT    interpolated transfer matrix (JG<-IG,ITEMP,IL).
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NG,NANI,NT0,NPSN0,FGTD,IAFAG(NG+1,NANI),IFAGR(NG,NANI)
      REAL TT,TEMP(NT0),PSN0(NPSN0),SCAT(NG,NG,NANI)
*----
*  LOCAL VARIABLES
*----
      CHARACTER HSMG*131
      PARAMETER (NINT=2,DTMIN=1.0)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DTEMP,WEIJHT,S
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: DSCATT
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(DTEMP(NT0),WEIJHT(NT0),DSCATT(NG,NG,NT0,NANI),S(NG))
*
      IF(NT0.EQ.1) THEN
        IPROX=1
        IGTFIX=1
      ELSE
        DO I=1,NT0
          DTEMP(I)=TEMP(I)
        ENDDO
        CALL LIBA28(TT,DTEMP,NT0,NINT,WEIJHT,IORD,IPROX,I0)
        IF(ABS(TT-TEMP(IPROX)).LE.DTMIN) THEN
          IGTFIX=1
        ELSE IF((TT.LT.TEMP(1)).OR.(TT.GT.TEMP(NT0))) THEN
          WRITE(HSMG,'(A,F8.2,A,F8.2,A,F8.2)')
     1    'LIBA33: A TEMPERATURE', TT,'K IS NOT INCLUDED BETWEEN ',
     2    TEMP(1),' AND ',TEMP(NT0)
          WRITE(6,'(/1X,A)') HSMG
          IGTFIX=2
        ELSE
          IGTFIX=0
        ENDIF
      ENDIF
*----
*  SCATTERING MATRIX RECONSTRUCTION
*----
      DSCATT(:NG,:NG,:NT0,:NANI)=0.D0
      NV=0
      DO IL=1,NANI
        DO IG=1,NG ! departure group
          JG1=IFAGR(IG,IL)+1
          ISIZE=IAFAG(IG+1,IL)-IAFAG(IG,IL)
          JG2=JG1+ISIZE-1
          IF(JG2.GT.NG) CALL XABORT('LIBA33: NG OVERFLOW(1)')
          IF(NV+ISIZE.GT.NPSN0) CALL XABORT('LIBA33: NPSN0 OVERFLOW(1)')
          DSCATT(JG1:JG2,IG,1,IL)=PSN0(NV+1:NV+ISIZE)/REAL(2*IL-1)
          NV=NV+ISIZE
        ENDDO
        IF(FGTD.GE.1) THEN
          DO IT=2,NT0
            DO IG=1,FGTD-1 ! departure group
              DSCATT(:NG,IG,IT,IL)=DSCATT(:NG,IG,1,IL)
            ENDDO
            DO IG=FGTD,NG ! departure group
              JG1=IFAGR(IG,IL)+1
              ISIZE=IAFAG(IG+1,IL)-IAFAG(IG,IL)
              JG2=JG1+ISIZE-1
              IF(JG2.GT.NG) CALL XABORT('LIBA33: NG OVERFLOW(2)')
              IF(NV+ISIZE.GT.NPSN0) CALL XABORT('LIBA33: NPSN0 OVERFLO'
     1        //'W(2)')
              DSCATT(JG1:JG2,IG,IT,IL)=PSN0(NV+1:NV+ISIZE)/REAL(2*IL-1)
              NV=NV+ISIZE
            ENDDO
          ENDDO
        ENDIF
      ENDDO
*----
*  TEMPERATURE INTERPOLATION
*----
      SCAT(:NG,:NG,:NANI)=0.0
      IF(FGTD.GE.1) THEN
        DO IL=1,NANI
          SCAT(:NG,:FGTD-1,IL)=REAL(DSCATT(:NG,:FGTD-1,1,IL))
        ENDDO
      ELSE
        DO IL=1,NANI
          SCAT(:NG,:NG,IL)=REAL(DSCATT(:NG,:NG,1,IL))
        ENDDO
        RETURN
      ENDIF
      IDIS=NG+1-FGTD
      DO IL=1,NANI
        IF(IGTFIX.EQ.1) THEN
          DO I=1,IDIS
            SCAT(:NG,FGTD+I-1,IL)=REAL(DSCATT(:NG,FGTD+I-1,IPROX,IL))
          ENDDO
        ELSE
          DO IG=FGTD,NG ! departure group
            S(:NG)=0.D0
            DO J=1,IORD ! temperature weighting
              S(:NG)=S(:NG)+WEIJHT(J)*DSCATT(:NG,IG,I0+J,IL)
            ENDDO
            IF(IGTFIX.EQ.2) THEN
              DO JG=1,NG ! arrival group
                IF(DSCATT(JG,IG,IPROX,IL).GE.0.) THEN
                  S(JG)=MAX(0.D0,S(JG))
                ELSE
                  S(JG)=MIN(S(JG),0.D0)
                ENDIF
              ENDDO
            ENDIF
            SCAT(:NG,IG,IL)=REAL(S(:NG))
          ENDDO
        ENDIF
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(S,DSCATT,WEIJHT,DTEMP)
      RETURN
      END
