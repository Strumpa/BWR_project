*DECK AUTDRV
      SUBROUTINE AUTDRV(IPLI0,IPTRK,IPLIB,IFTRAK,INDREC,CDOOR,IMPX,
     1 IGRMIN,IGRMAX,NGRP,NBMIX,NREG,NUN,NBISO,NL,NED,NDEL,LEAKSW,
     2 ITRANC,IPHASE,TITR,KSPH,NRES,NPASS,ICALC,IALTER,MAXTRA,ISEED,
     3 DIL,DELI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Driver for a resonance self-shielding calculation with the Autosecol
* method.
*
*Copyright:
* Copyright (C) 2023 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPLI0   pointer to the internal microscopic cross section library
*         builded by the self-shielding module (L_LIBRARY signature).
* IPTRK   pointer to the tracking (L_TRACK signature).
* IPLIB   pointer to the internal microscopic cross section library
*         with subgroups (L_LIBRARY signature).
* IFTRAK  unit number of the sequential binary tracking file.
* INDREC  access flag for the internal microscopic cross section library
*         builded by the self-shielding module (=1 IPLI0 access in
*         creation mode; =2 in modification mode).
* CDOOR   name of the geometry/solution operator.
* IMPX    print flag (equal to zero for no print).
* IGRMIN  first group where the self-shielding is applied.
* IGRMAX  most thermal group where the self-shielding is applied.
* NGRP    number of energy groups.
* NBMIX   number of mixtures in the internal library.
* NREG    number of regions.
* NUN     number of unknowns per energy group.
* NBISO   number of isotopes specifications in the internal library.
* NL      number of Legendre orders required in the calculation
*         (NL=1 or higher).
* NED     number of extra vector edits.
* NDEL    number of delayed neutron precursor groups.
* LEAKSW  leakage flag (LEAKSW=.TRUE. if neutron leakage through
*         external boundary is present).
* ITRANC  type of transport correction.
* IPHASE  type of flux solution (=1 use a native flux solution door;
*         =2 use collision probabilities).
* TITR    title.
* KSPH    SPH equivalence flag (=0 no SPH correction; =1 SPH correction
*         in the fuel).
* NRES    number of self-shielding zones, as given by LIB:.
* NPASS   number of outer iterations.
* ICALC   simplified self-shielding flag (=1 IPLI0 is containing ICALC
*         data. =0 no ICALC data).
* IALTER  type of elastic slowing-down kernel (=0: use exact kernel;
*         =1: use an approximate kernel for the resonant isotopes).
* MAXTRA  maximum number of down-scattering terms.
* ISEED   the seed for the generation of random numbers in the
*         unresolved energy domain.
* DIL     microscopic dilution cross section of each isotope.
*
*Parameters: output
* DELI    elementary lethargy width used by the elastic kernel.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLI0,IPTRK,IPLIB
      INTEGER IFTRAK,INDREC,IMPX,IGRMIN,IGRMAX,NGRP,NBMIX,NREG,NUN,
     1 NBISO,NL,NED,NDEL,ITRANC,IPHASE,KSPH,NRES,NPASS,ICALC,IALTER,
     2 MAXTRA,ISEED
      REAL DIL(NBISO),DELI
      CHARACTER CDOOR*12,TITR*72
      LOGICAL LEAKSW
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40,MAXRSS=300,MAXESP=4)
      TYPE(C_PTR) JPLI0,KPLI0,JPLIB,KPLIB
      CHARACTER HSMG*131,TEXT4*4,NAM1*4,FNAM1*4,NAM2*12,FNAM2*12,
     1 TEXT8*8
      INTEGER IPAR(NSTATE),IRSS(MAXRSS),IESP(MAXESP+1)
      REAL TMPDAY(3),EESP(MAXESP+1)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MAT,KEYFLX,MIX,IEVOL,ITYPE,
     1 LSHI,IHSUF,ILLIB,JCEDM,NBIN,NBIN_AU
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISONAM,IHLIB
      REAL, ALLOCATABLE, DIMENSION(:) :: VOL,TN,DEN,ENER,DELTAU,EBIN,GS,
     1 VOLMIX
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: MASK,MASKL
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: IPISO1
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(MAT(NREG),KEYFLX(NREG),ISONAM(3,NBISO),MIX(NBISO),
     1 IEVOL(NBISO),ITYPE(NBISO),LSHI(NBISO),IHSUF(NBISO),
     2 IHLIB(2,NBISO),ILLIB(NBISO),IPISO1(NBISO))
      ALLOCATE(VOL(NREG),TN(NBISO),DEN(NBISO),ENER(NGRP+1))
*----
*  RECOVER USEFUL INFORMATION FROM TRACKING OBJECT.
*----
      CALL LCMGET(IPTRK,'MATCOD',MAT)
      CALL LCMGET(IPTRK,'VOLUME',VOL)
      CALL LCMGET(IPTRK,'KEYFLX',KEYFLX)
*----
*  RECOVER USEFUL INFORMATION FROM LIBRARY OBJECTS.
*----
      CALL LCMGET(IPLIB,'ISOTOPESUSED',ISONAM)
      CALL LCMGET(IPLIB,'ISOTOPESMIX',MIX)
      CALL LCMGET(IPLIB,'ISOTOPESTODO',IEVOL)
      CALL LCMGET(IPLIB,'ISOTOPESTYPE',ITYPE)
      CALL LCMGET(IPLIB,'ISOTOPESTEMP',TN)
*
      CALL LCMPUT(IPLI0,'ISOTOPESMIX',NBISO,1,MIX)
      CALL LCMPUT(IPLI0,'ISOTOPESTODO',NBISO,1,IEVOL)
      CALL LCMPUT(IPLI0,'ISOTOPESTYPE',NBISO,1,ITYPE)
      CALL LCMPUT(IPLI0,'ISOTOPESTEMP',NBISO,2,TN)
      IF(INDREC.EQ.1) THEN
         CALL LCMGET(IPLIB,'ISOTOPESDENS',DEN)
         CALL LCMPUT(IPLI0,'ISOTOPESDENS',NBISO,2,DEN)
      ELSE IF(INDREC.EQ.2) THEN
         CALL LCMGET(IPLI0,'ISOTOPESDENS',DEN)
      ENDIF
      CALL LCMGET(IPLIB,'ISOTOPESSHI',LSHI)
      CALL LCMLEN(IPLIB,'ISOTOPESDSN',NELSN,ITYLCM)
      IF(NELSN.GT.0) THEN
        NGIS=NGRP*NBISO
        ALLOCATE(GS(NGIS))
        CALL LCMGET(IPLIB,'ISOTOPESDSN',GS)
        CALL LCMPUT(IPLI0,'ISOTOPESDSN',NGIS,2,GS)
        CALL LCMGET(IPLIB,'ISOTOPESDSB',GS)
        CALL LCMPUT(IPLI0,'ISOTOPESDSB',NGIS,2,GS)
        DEALLOCATE(GS)
      ENDIF
      CALL LCMGET(IPLIB,'DELTAU',ENER)
      CALL LCMPUT(IPLI0,'DELTAU',NGRP,2,ENER)
      CALL LCMGET(IPLIB,'ENERGY',ENER)
      CALL LCMPUT(IPLI0,'ENERGY',NGRP+1,2,ENER)
      CALL LCMLEN(IPLIB,'CHI-LIMITS',NBESP,ITYLCM)
      IF(NBESP.GT.0) THEN
         NBESP=NBESP-1
         IF(NBESP.GT.MAXESP) CALL XABORT('AUTDRV: MAXESP OVERFLOW.')
         CALL LCMGET(IPLIB,'CHI-LIMITS',IESP)
         CALL LCMPUT(IPLI0,'CHI-LIMITS',NBESP+1,1,IESP)
         CALL LCMGET(IPLIB,'CHI-ENERGY',EESP)
         CALL LCMPUT(IPLI0,'CHI-ENERGY',NBESP+1,2,EESP)
      ENDIF
*----
*  COMPUTE MIXTURESVOL.
*----
      ALLOCATE(VOLMIX(NBMIX))
      VOLMIX(:NBMIX)=0.0
      DO I=1,NREG
        IBM=MAT(I)
        IF(IBM.GT.0) VOLMIX(IBM)=VOLMIX(IBM)+VOL(I)
        CALL LCMPUT(IPLI0,'MIXTURESVOL',NBMIX,2,VOLMIX)
      ENDDO
      DEALLOCATE(VOLMIX)
*----
*  RECOVER BIN TYPE INFORMATION (IF AVAILABLE).
*  ASSUME THAT THE ELEMENTARY LETHARGY WIDTH DELI IS A RATIONAL FRACTION
*  OF THE LETHARGY UNIT. CHECK THIS ASSUMPTION.
*----
      CALL LIBIPS(IPLIB,NBISO,IPISO1)
      ALLOCATE(NBIN(NGRP),NBIN_AU(NGRP))
      NBIN(:NGRP)=-1
      DO ISO=1,NBISO
        KPLIB=IPISO1(ISO) ! set ISO-th isotope
        IF(.NOT.C_ASSOCIATED(KPLIB)) THEN
          WRITE(HSMG,'(17HAUTDRV: ISOTOPE '',3A4,17H'' IS NOT ASSOCIAT,
     1    3HED.)') ISONAM(:3,ISO)
          CALL XABORT(HSMG)
        ENDIF
        CALL LCMLEN(KPLIB,'BIN-NFS',LENGT,ITYLCM)
        IF(LENGT.EQ.0) CYCLE
        CALL LCMGET(KPLIB,'BIN-NFS',NBIN_AU)
        DO 10 IGRP=NGRP,1,-1
        IF(NBIN(IGRP).EQ.-1) THEN
          NBIN(IGRP)=NBIN_AU(IGRP)
        ELSE IF(NBIN(IGRP).NE.NBIN_AU(IGRP)) THEN
          WRITE(HSMG,'(38HAUTDRV: INCONSISTENT BIN DATA IN GROUP,I4,
     1    13H OF ISOTOPE '',3A4,2H''.)') IGRP,ISONAM(:3,ISO)
          CALL XABORT(HSMG)
        ENDIF
   10   CONTINUE
      ENDDO
      DEALLOCATE(NBIN_AU)
      IGRRES=IGRMIN
      DO 20 IGRP=IGRMIN,IGRMAX
      IGRRES=IGRP
      IF(NBIN(IGRP).GT.0) GO TO 30
   20 CONTINUE
   30 DO 40 IGRP=1,NGRP
      IF(NBIN(IGRP).EQ.-1) CALL XABORT('AUTDRV: NBIN SETTING FAILURE.')
   40 CONTINUE
      LBIN=SUM(NBIN(:NGRP))
      IF(LBIN.EQ.0) CALL XABORT('AUTDRV: NO AUTOLIB DATA.')
      ALLOCATE(DELTAU(LBIN),EBIN(LBIN+1))
      DO ISO=1,NBISO
        KPLIB=IPISO1(ISO) ! set ISO-th isotope
        CALL LCMLEN(KPLIB,'BIN-NFS',LENGT,ITYLCM)
        IF(LENGT.EQ.0) CYCLE
        CALL LCMGET(KPLIB,'BIN-ENERGY',EBIN)
        EXIT
      ENDDO
      DELMIN=1.0E10
      IBIN=0
      DO 60 IGRP=1,NGRP
      DO 50 IGF=1,NBIN(IGRP)
      DELM=LOG(EBIN(IBIN+IGF)/EBIN(IBIN+IGF+1))
      DELMIN=MIN(DELMIN,DELM)
      DELTAU(IBIN+IGF)=DELM
   50 CONTINUE
      IBIN=IBIN+NBIN(IGRP)
   60 CONTINUE
      CALL LCMLEN(KPLIB,'BIN-DELI',LENGT,ITYLCM)
      IF((LENGT.EQ.1).AND.(ITYLCM.EQ.2)) THEN
        CALL LCMGET(KPLIB,'BIN-DELI',DELI)
      ELSE
        DELI=1.0/REAL(INT(1.00001/DELMIN))
      ENDIF
      IBIN=0
      ERR=0.0
      DO 80 IGRP=1,NGRP
      DO 70 IGF=1,NBIN(IGRP)
      LARGH=INT(DELTAU(IBIN+IGF)/DELI+0.1)
      ERR=MAX(ERR,ABS(DELTAU(IBIN+IGF)/DELI-REAL(LARGH)))
   70 CONTINUE
      IBIN=IBIN+NBIN(IGRP)
   80 CONTINUE
      IF(ERR.GT.0.05) THEN
         WRITE(HSMG,'(45HAUTDRV: UNABLE TO SET THE ELEMENTARY LETHARGY,
     1   7H WIDTH.)')
         WRITE(6,'(A)') HSMG
      ENDIF
      DEALLOCATE(EBIN,DELTAU)
*----
*  RECOMPUTE THE AUTOLIB ENERGY MESH BETWEEN IGRMIN AND IGRMAX.
*----
      DO 90 IGRP=IGRMIN,IGRRES-1
      IF(NBIN(IGRP).EQ.0) THEN
        DELM=LOG(ENER(IGRP)/ENER(IGRP+1))
        NBIN(IGRP)=INT(DELM/DELI+0.1)
      ENDIF
   90 CONTINUE
      LBIN=SUM(NBIN(IGRMIN:NGRP))
      IF(IMPX.GT.0) THEN
        WRITE(6,'(/32H AUTDRV: NUMBER OF AUTOLIB BINS=,I9)') LBIN
        WRITE(6,'(35H AUTDRV: FIRST SELF-SHIELDED GROUP=,I6)') IGRMIN
        WRITE(6,'(30H AUTDRV: FIRST RESOLVED GROUP=,I11)') IGRRES
        WRITE(6,'(34H AUTDRV: LAST SELF-SHIELDED GROUP=,I7)') IGRMAX
        WRITE(6,'(35H AUTDRV: ELEMENTARY LETHARGY WIDTH=,1P,E9.2)') DELI
        WRITE(6,'(33H AUTDRV: ERROR ON LETHARGY WIDTH=,1P,E11.2)') ERR
      ENDIF
      ALLOCATE(EBIN(LBIN+1))
      EBIN(:LBIN+1)=0.0
      LLL=0
      DO IGRP=1,IGRRES-1
        DUUU=0.0D0
        EBIN(LLL+1)=ENER(IGRP)
        DO IBIN=LLL+1,LLL+NBIN(IGRP)
          DUUU=DUUU+DELI
          EBIN(IBIN+1)=REAL(ENER(IGRP)*EXP(-DUUU))
        ENDDO
        LLL=LLL+NBIN(IGRP)
      ENDDO
      DO ISO=1,NBISO
        KPLIB=IPISO1(ISO) ! set ISO-th isotope
        CALL LCMLEN(KPLIB,'BIN-ENERGY',LENGT,ITYLCM)
        IF(LENGT.EQ.0) CYCLE
        IF(LLL+LENGT.GT.LBIN+1) CALL XABORT('AUTDRV: EBIN OVERFLOW.')
        CALL LCMGET(KPLIB,'BIN-ENERGY',EBIN(LLL+1))
        IF(EBIN(LBIN+1).EQ.0.0) EBIN(LBIN+1)=1.0E-5
        EXIT
      ENDDO
*
      DO 100 ISO=1,NBISO
      TEXT8='MICROLIB'
      READ(TEXT8,'(2A4)') IHLIB(1,ISO),IHLIB(2,ISO)
      ILLIB(ISO)=1
  100 CONTINUE
*
      JPLIB=LCMGID(IPLIB,'ISOTOPESLIST')
      JPLI0=LCMLID(IPLI0,'ISOTOPESLIST',NBISO)
      IF(INDREC.EQ.1) THEN
*        COPY THE NON RESONANT ISOTOPES.
         CALL KDRCPU(TK1)
         DO 130 ISO=1,NBISO
         IF((LSHI(ISO).EQ.0).OR.(DEN(ISO).EQ.0.0)) THEN
            CALL LCMLEL(JPLIB,ISO,ILEN,ITYLCM)
            IF(ILEN.EQ.0) THEN
              DO JSO=1,ISO-1
                CALL LCMLEL(JPLIB,JSO,ILEN,ITYLCM)
                IF(ILEN.EQ.0) CYCLE
                IF((ISONAM(1,ISO).EQ.ISONAM(1,JSO)).AND.(ISONAM(2,ISO)
     1          .EQ.ISONAM(2,JSO)).AND.(ISONAM(3,ISO).EQ.ISONAM(3,JSO)))
     2          THEN
                  IF(LSHI(JSO).GT.0) THEN
                    KPLIB=LCMGIL(JPLIB,JSO) ! set JSO-th isotope
                    GO TO 120
                  ELSE
                    GO TO 130
                  ENDIF
                ENDIF
              ENDDO
            ELSE
              KPLIB=LCMGIL(JPLIB,ISO) ! set ISO-th isotope
              GO TO 120
            ENDIF
            GO TO 130
  120       CALL LCMLEL(JPLI0,ISO,ILEN,ITYLCM)
            IF(ILEN.NE.0) GO TO 130
            KPLI0=LCMDIL(JPLI0,ISO) ! set ISO-th isotope
            CALL LCMEQU(KPLIB,KPLI0)
         ENDIF
  130    CONTINUE
         CALL KDRCPU(TK2)
         IF(IMPX.GT.1) WRITE(6,'(/33H AUTDRV: CPU TIME SPENT TO COPY T,
     1   26HHE NON-RESONANT ISOTOPES =,F8.1,8H SECOND.)') TK2-TK1
*
*        WRITE THE OUTPUT INTERNAL LIBRARY PARAMETERS.
         CALL LCMGET(IPLIB,'STATE-VECTOR',IPAR)
         IPAR(8)=0
         IPAR(17)=0
         CALL LCMPUT(IPLI0,'STATE-VECTOR',NSTATE,1,IPAR)
         IF(NED.GT.0) THEN
            ALLOCATE(JCEDM(2*NED))
            CALL LCMGET(IPLIB,'ADDXSNAME-P0',JCEDM)
            CALL LCMPUT(IPLI0,'ADDXSNAME-P0',2*NED,3,JCEDM)
            DEALLOCATE(JCEDM)
         ENDIF
         CALL LCMLEN(IPLIB,'DEPL-CHAIN',ILENG,ITYLCM)
         IF(ILENG.NE.0) THEN
            CALL LCMSIX(IPLIB,'DEPL-CHAIN',1)
            CALL LCMSIX(IPLI0,'DEPL-CHAIN',1)
            CALL LCMEQU(IPLIB,IPLI0)
            CALL LCMSIX(IPLI0,' ',2)
            CALL LCMSIX(IPLIB,' ',2)
         ENDIF
      ENDIF
      IF(NRES.EQ.0) THEN
        CALL LCMEQU(IPLIB,IPLI0)
        GO TO 310
      ENDIF
*----
* FIND THE ISOTOPE-NAME SUFFIX VALUES.
*----
      TEXT4=' '
      READ(TEXT4,'(A4)') IHBLK
      DO 140 ISO=1,NBISO
      IF((LSHI(ISO).NE.0).AND.(DEN(ISO).NE.0.0)) THEN
         WRITE(TEXT4,'(I4.4)') MIX(ISO)
         READ(TEXT4,'(A4)') IHSUF(ISO)
      ELSE
         IHSUF(ISO)=IHBLK
      ENDIF
  140 CONTINUE
      IF(ICALC.EQ.1) THEN
         CALL LCMSIX(IPLI0,'SHIBA_SG',1)
         CALL LCMSIX(IPLI0,'-DATA-CALC-',1)
         NAM1=' '
         CALL LCMNXT(IPLI0,NAM1)
         FNAM1=NAM1
  150    CALL LCMSIX(IPLI0,NAM1,1)
            NAM2=' '
            CALL LCMNXT(IPLI0,NAM2)
            FNAM2=NAM2
  160       CALL LCMLEN(IPLI0,NAM2,NRSS,ITYLCM)
            CALL LCMGET(IPLI0,NAM2,IRSS)
            READ(NAM2,'(2A4)') IN1,IN2
            DO 180 ISO=1,NBISO
            IF((ISONAM(1,ISO).EQ.IN1).AND.(ISONAM(2,ISO).EQ.IN2).AND.
     1      (LSHI(ISO).NE.0)) THEN
               IF((NRSS.EQ.1).AND.(IRSS(1).EQ.-999)) THEN
                  READ(NAM1,'(A4)') IHSUF(ISO)
               ELSE
                  DO 170 I=1,NRSS
                  IF(IRSS(I).EQ.MIX(ISO)) READ(NAM1,'(A4)') IHSUF(ISO)
  170             CONTINUE
               ENDIF
            ENDIF
  180       CONTINUE
            CALL LCMNXT(IPLI0,NAM2)
            IF(NAM2.EQ.FNAM2) GO TO 190
            GO TO 160
  190    CALL LCMSIX(IPLI0,' ',2)
         CALL LCMNXT(IPLI0,NAM1)
         IF(NAM1.EQ.FNAM1) THEN
            CALL LCMSIX(IPLI0,' ',2)
            CALL LCMSIX(IPLI0,' ',2)
            GO TO 200
         ENDIF
         GO TO 150
      ENDIF
*
  200 NPASS2=NPASS
      DO 300 IPASS=1,NPASS
      IF((IMPX.GT.0).AND.(NPASS2.GT.1)) WRITE (6,'(/15H AUTDRV: SELF S,
     1 25HHIELDING ITERATION NUMBER,I4,7H  NRES=,I4,1H.)') IPASS,NRES
      DO 290 INRS=1,NRES
*----
*  PERFORM A SELF-SHIELDING CALCULATION IN RESONANT REGION INRS.
*----
      CALL AUTONE(IPLI0,IPTRK,IPLIB,IFTRAK,CDOOR,IMPX,INRS,IGRMIN,
     1 IGRRES,IGRMAX,NGRP,NBMIX,NREG,NUN,NBISO,NL,NED,NDEL,ISONAM,
     2 IHSUF,DEN,LSHI,DIL,MIX,MAT,VOL,KEYFLX,LEAKSW,ITRANC,IPHASE,
     3 TITR,KSPH,IALTER,DELI,LBIN,NBIN,EBIN,MAXTRA,ISEED)
  290 CONTINUE
  300 CONTINUE
  310 IF(IMPX.GE.3) CALL LCMLIB(IPLI0)
*----
*  BUILD THE MACROLIB IN THE OUTPUT INTERNAL LIBRARY.
*----
      ALLOCATE(MASK(NBMIX))
      DO 330 IBM=1,NBMIX
      MASK(IBM)=.TRUE.
      DO 320 I=1,NREG
      IF(MAT(I).EQ.IBM) GO TO 330
  320 CONTINUE
      MASK(IBM)=.FALSE.
  330 CONTINUE
      ALLOCATE(MASKL(NGRP))
      DO 340 I=1,NGRP
      MASKL(I)=.TRUE.
  340 CONTINUE
*
      ITSTMP=0
      TMPDAY(1)=0.0
      TMPDAY(2)=0.0
      TMPDAY(3)=0.0
      CALL KDRCPU(TK1)
      CALL LCMLEN(IPLI0,'ISOTOPESUSED',ILENG,ITYLCM)
      IF(ILENG.EQ.0) CALL XABORT('AUTDRV: MISSING ISOTOPESUSED RECORD.')
      CALL LCMGET(IPLI0,'ISOTOPESUSED',ISONAM)
      CALL LIBMIX(IPLI0,NBMIX,NGRP,NBISO,ISONAM,MIX,DEN,MASK,MASKL,
     1 ITSTMP,TMPDAY)
      CALL KDRCPU(TK2)
      IF(IMPX.GT.1) WRITE(6,'(/37H AUTDRV: CPU TIME SPENT TO BUILD THE ,
     1 19HEMBEDDED MACROLIB =,F8.1,8H SECOND.)') TK2-TK1
      DEALLOCATE(MASKL,MASK)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(NBIN,EBIN,ENER,DEN,TN,VOL)
      DEALLOCATE(IPISO1,ILLIB,IHLIB,IHSUF,LSHI,ITYPE,IEVOL,MIX,ISONAM,
     1 KEYFLX,MAT)
      RETURN
      END
