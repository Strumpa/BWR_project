*DECK LIBADD
      SUBROUTINE LIBADD (IPLIB,NBISO,MASKI,IMPX,NGRO,NL,ITRANC,NDEPL,
     1 ISONAM,ISONRF,IPISO,NIR,GIR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Add transport correction, Goldstein-Cohen and H-FACTOR data to a
* /microlib/ directory.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*  09/2026 proposed modification by R. Guasch :
*    - flag to exclude energy from lumped isotopes and garantee 
*      edepomode 0 Serpent2 formulation from R. Tuominen et al., ANE 129 (2019) 224–232.
*
*Parameters: input
* IPLIB   pointer to the lattice microscopic cross section library
*         (L_LIBRARY signature).
* NBISO   number of isotopes present in the calculation domain.
* MASKI   isotopic mask. Isotope with index I is processed if
*         MASKI(I)=.true.
* IMPX    print flag.
* NGRO    number of energy groups.
* NL      number of Legendre orders required in the calculation
*         NL=1 (for isotropic scattering) or higher.
* ITRANC  transport correction option (=0: no correction; =1: Apollo-
*         type; =2: recover TRANC record; =3: Wims-type; =4: leakage
*         correction alone).
* NDEPL   number of depleting isotopes.
* ISONAM  alias name of each isotope.
* ISONRF  library reference name of each isotope.
* IPISO   pointer array towards microlib isotopes.
* NIR     group index with an imposed IR slowing-down model (=0 for no
*         IR model).
* GIR     value of the imposed Goldstein-Cohen parameter for groups
*         with an IR model.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB,IPISO(NBISO)
      INTEGER NBISO,IMPX,NGRO,NL,ITRANC,NDEPL,ISONAM(3,NBISO),
     1 ISONRF(3,NBISO),NIR(NBISO),ELUMPFLAG
      LOGICAL MASKI(NBISO)
      REAL GIR(NBISO)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (IOUT=6,NSTATE=40)
      INTEGER ISTATE(NSTATE)
      TYPE(C_PTR) JPLIB,KPLIB
      CHARACTER HSONAM*12,HSONRF*12,HSMG*131
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: WORK,WR2,DELTA
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SCAT,RER
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: HREAC
      CHARACTER(LEN=12), ALLOCATABLE, DIMENSION(:) :: HGAR
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(WORK(NGRO),WR2(NGRO),SCAT(NGRO,NGRO),DELTA(NGRO+1))
*----
*  RECOVER THE ENERGY GRID.
*----
      CALL LCMGET(IPLIB,'STATE-VECTOR',ISTATE)
      ELUMPFLAG=ISTATE(28)
      CALL LCMLEN(IPLIB,'ENERGY',LENGT,ITYLCM)
      IF(LENGT.EQ.0) CALL XABORT('LIBADD: NO GROUP STRUCTURE AVAILABLE')
      CALL LCMGET(IPLIB,'ENERGY',DELTA)
      NGX=0
      DO 10 IGR=1,NGRO
      IF((NGX.EQ.0).AND.(DELTA(IGR+1).LT.4.0)) NGX=IGR-1
  10  CONTINUE
      DO 15 IGR=1,NGRO
      DELTA(IGR)=LOG(DELTA(IGR)/DELTA(IGR+1))
  15  CONTINUE
*----
*  RECOVER DEPLETION DATA.
*----
      NREAC=0
      IF(NDEPL.NE.0) THEN
        CALL LCMLEN(IPLIB,'DEPL-CHAIN',LENGTH,ITYLCM)
        IF(LENGTH.EQ.0) THEN
          CALL LCMLIB(IPLIB)
          CALL XABORT('LIBADD: MISSING DEPL-CHAIN DATA.')
        ENDIF
        CALL LCMSIX(IPLIB,'DEPL-CHAIN',1)
        CALL LCMGET(IPLIB,'STATE-VECTOR',ISTATE)
        IF(ISTATE(1).NE.NDEPL) CALL XABORT('LIBADD: INVALID NUMBER OF '
     1   //'DEPLETING ISOTOPES.')
        NREAC=ISTATE(8)
        ALLOCATE(HGAR(NDEPL),RER(NREAC,NDEPL),HREAC(NREAC))
        CALL LCMGTC(IPLIB,'ISOTOPESDEPL',12,NDEPL,HGAR)
        CALL LCMGET(IPLIB,'DEPLETE-ENER',RER)
        CALL LCMGTC(IPLIB,'DEPLETE-IDEN',8,NREAC,HREAC)
        CALL LCMSIX(IPLIB,' ',2)
      ENDIF
*
      DO 110 ISO=1,NBISO
      IF(MASKI(ISO)) THEN
        WRITE(HSONAM,'(3A4)') (ISONAM(I,ISO),I=1,3)
        WRITE(HSONRF,'(3A4)') (ISONRF(I,ISO),I=1,3)
        KPLIB=IPISO(ISO) ! set ISO-th isotope
        IF(.NOT.C_ASSOCIATED(KPLIB)) GO TO 110
        CALL LCMLEN(KPLIB,'NTOT0',ILENG,ITYLCM)
        IF(ILENG.EQ.0) THEN
          JPLIB=LCMGID(IPLIB,'ISOTOPESLIST')
          CALL LCMLIB(JPLIB)
          WRITE(HSMG,'(17H LIBADD: ISOTOPE ,A12,6H (ISO=,I6,
     1    17H) IS NOT DEFINED.)')  HSONAM,ISO
          CALL XABORT(HSMG)
        ENDIF
*
*       REDIFINE THE GOLDSTEIN-COHEN PARAMETERS.
        IF(NIR(ISO).GT.0) THEN
           DO 20 IGR=1,MIN(NGRO,NIR(ISO)-1)
           WORK(IGR)=1.0
   20      CONTINUE
           DO 30 IGR=NIR(ISO),NGRO
           WORK(IGR)=GIR(ISO)
   30      CONTINUE
           CALL LCMPUT(KPLIB,'NGOLD',NGRO,2,WORK)
           IF(IMPX.GT.1) THEN
              IF(GIR(ISO).EQ.-998.0) THEN
                 WRITE(IOUT,210) HSONAM,'PT',NIR(ISO)
              ELSE IF(GIR(ISO).EQ.-999.0) THEN
                 WRITE(IOUT,210) HSONAM,'PTSL',NIR(ISO)
              ELSE IF(GIR(ISO).EQ.-1000.0) THEN
                 WRITE(IOUT,210) HSONAM,'PTMC',NIR(ISO)
              ELSE
                 WRITE(IOUT,200) HSONAM,GIR(ISO),NIR(ISO)
              ENDIF
           ENDIF
        ENDIF
*
*       COMPUTE OR RECOVER THE TRANSPORT CORRECTION.
        IF(ITRANC.EQ.2) THEN
*          RECOVER THE TRANSPORT CORRECTION FROM THE LIBRARY.
           CALL LCMLEN(KPLIB,'TRANC',ILENG,ITYLCM)
           IF(ILENG.EQ.0) THEN
              WORK(:NGRO)=0.0
              CALL LCMPUT(KPLIB,'TRANC',NGRO,2,WORK)
           ENDIF
        ELSE IF(ITRANC.NE.0) THEN
           WORK(:NGRO)=0.0
           CALL LCMLEN(KPLIB,'NTOT1',ILENG,ITYLCM)
           IF(ILENG.NE.0) THEN
*             LEAKAGE CORRECTION.
              CALL LCMGET(KPLIB,'NTOT1',WORK)
              CALL LCMGET(KPLIB,'NTOT0',WR2)
              DO 40 IG1=1,NGRO
              WORK(IG1)=WR2(IG1)-WORK(IG1)
   40         CONTINUE
           ENDIF
           IF((NL.GE.2).AND.(ITRANC.NE.4)) THEN
              CALL LCMLEN(KPLIB,'SCAT-SAVED',ILENG,ITYLCM)
              IF(ILENG.EQ.0) THEN
                 WRITE(HSMG,'(37H LIBADD: NO SCAT-SAVED RECORD FOR ISO,
     1           5HTOPE ,A12,1H.)') HSONAM
                 CALL XABORT(HSMG)
              ENDIF
              CALL XDRLGS(KPLIB,-1,0,1,1,1,NGRO,WR2,SCAT,ITY)
              IF(ITRANC.EQ.1) THEN
*                APOLLO-TYPE TRANSPORT CORRECTION. USE THE MICRO-
*                REVERSIBILITY PRINCIPLE AT ALL ENERGIES.
                 DO 50 IG1=1,NGRO
                 WORK(IG1)=WORK(IG1)+WR2(IG1)
   50            CONTINUE
              ELSE IF(ITRANC.EQ.3) THEN
*                WIMS-TYPE TRANSPORT CORRECTION. USE THE MICRO-
*                REVERSIBILITY PRINCIPLE BELOW 4 EV AND A 1/E SPECTRUM
*                ABOVE.
                 DO 65 IG1=1,MIN(NGRO,NGX)
                 DO 60 IG2=1,NGRO
                 WORK(IG1)=WORK(IG1)+SCAT(IG1,IG2)*DELTA(IG2)/DELTA(IG1)
   60            CONTINUE
   65            CONTINUE
                 DO 70 IG1=NGX+1,NGRO
                 WORK(IG1)=WORK(IG1)+WR2(IG1)
   70            CONTINUE
              ELSE
                 CALL XABORT('LIBADD: UNKNOWN TYPE OF CORRECTION.')
              ENDIF
           ENDIF
*          ***CAUTION*** 'TRANC' CONTAINS BOTH TRANSPORT AND LEAKAGE
*          CORRECTIONS.
           CALL LCMPUT(KPLIB,'TRANC',NGRO,2,WORK)
        ENDIF
*
*       ADD OR CORRECT H-FACTOR INFORMATION IN THE MICROLIB.
        IF(NDEPL.NE.0) THEN
          JDEPL=0
          DO IDEPL=1,NDEPL
           JDEPL=IDEPL
            IF(HSONRF.EQ.HGAR(IDEPL)) GO TO 80
          ENDDO
          CYCLE
   80     WORK(:NGRO)=0.0
          CALL LCMLEN(KPLIB,'H-FACTOR',LENGTH,ITYLCM)
          IF(LENGTH.NE.0) CALL LCMGET(KPLIB,'H-FACTOR',WORK)
          DO IREA=2,NREAC
            CALL LCMLEN(KPLIB,HREAC(IREA),LENGTH,ITYLCM)
            IF(LENGTH.EQ.0) CYCLE
            IF(LENGTH.GT.NGRO) CALL XABORT('LIBADD: WR2 OVERFLOW.')
            WR2(:NGRO)=0.0
            CALL LCMGET(KPLIB,HREAC(IREA),WR2)
            IF(ELUMPFLAG.EQ.0) THEN
             DO IG=1,LENGTH
              WORK(IG)=WORK(IG)+RER(IREA,JDEPL)*WR2(IG)*1.0E6
             ENDDO
            ENDIF
          ENDDO ! IREA
          CALL LCMPUT(KPLIB,'H-FACTOR',NGRO,2,WORK)
          IF(IMPX.GT.1) THEN
            WRITE(IOUT,'(42H LIBADD: ADD H-FACTOR INFORMATION TO ISOTO,
     1      3HPE ,A,1H.)') TRIM(HSONRF)
          ENDIF
        ENDIF
      ENDIF
  110 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      IF(NDEPL.NE.0) DEALLOCATE(HREAC,RER,HGAR)
      DEALLOCATE(DELTA,SCAT,WR2,WORK)
      RETURN
*
  200 FORMAT(/51H LIBADD: THE GOLDSTEIN-COHEN PARAMETER OF ISOTOPE ',
     1 A12,12H' WAS SET TO,F5.2,33H FOR GROUPS WITH INDEX GREATER OR,
     2 9H EQUAL TO,I4,1H.)
  210 FORMAT(/18H LIBADD: ISOTOPE ',A12,20H' IS PROCESSED WITH ,A,
     1 48H METHOD IN GROUPS WITH INDEX GREATER OR EQUAL TO,I4,1H.)
      END
