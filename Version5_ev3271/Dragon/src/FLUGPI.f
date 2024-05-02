*DECK FLUGPI
      SUBROUTINE FLUGPI(IPFLUX,IPMACR,ITYPEC,MAXOUT,MAXINR,EPSOUT,
     >                  EPSUNK,EPSINR,IREBAL,IFRITR,IACITR,COPTIO,
     >                  ILEAK,B2,NGROUP,NREGIO,NMAT,NIFISS,LEAKSW,
     >                  REFKEF,ITPIJ,IPRINT,REC,INITFL,NMERG,IMERG)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read data for flux solution operator.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input
* IPFLUX  pointer to the flux LCM object.
* IPMACR  pointer to the macrolib LCM object.
* NGROUP  number of energy groups.
* NREGIO  number of regions.
* NMAT    number of material mixtures.
* NIFISS  number of fissile isotopes.
* LEAKSW  leakage flag (=.true. if leakage is present on the outer
*         surface).
* ITPIJ   type of collision probability information available:
*         =1 scattering modified pij (wij);
*         =2 standard pij;
*         =3 scattering modified pij+pijk (wij,wijk);
*         =4 standard pij+pijk.
* REC     flux recovery flag:
*         =.true. recover the existing solution as initial estimate;
*         =.false. use a new initial estimate.
*
*Parameters: output
* ITYPEC  type of flux evaluation:
*         =-1 skip the flux calculation;
*         =0  fixed sources;
*         =1  fixed source eigenvalue problem (GPT type);
*         =2  fission sources/K-eff convergence;
*         =3  fission sources/K-eff convergence/db2 buckling evaluation;
*         =4  fission sources/db2 buckling convergence;
*         =5  b2 sources/db2 buckling convergence;
* MAXOUT  maximum number of outer iterations.
* MAXINR  maximum number of thermal iterations.
* EPSOUT  outer iterations eigenvalue epsilon.
* EPSUNK  outer iterations eigenvector epsilon.
* EPSINR  thermal iterations epsilon.
* IREBAL  flux rebalancing flag (=1: perform rebalancing).
* IFRITR  number of free iterations in an acceleration cycle.
* IACITR  number of accelerated iterations in an acceleration cycle.
* COPTIO  type of leakage coefficients: 
*         'LKRD' (recover leakage coefficients in Macrolib);
*         'RHS' (recover leakage coefficients in RHS flux object);
*         'B0' (B-0), 'P0' (P-0), 'B1' (B-1),
*         'P1' (P-1), 'B0TR' (B-0 with transport correction) or 'P0TR'
*         (P-0 with transport correction).
* ILEAK   method used to include db2 effect:
*         =1 the scattering modified cp matrix is multiplied by PNLR;
*         =2 the reduced cp matrix is multiplied by PNL;
*         =3 sigs0-db2 approximation;
*         =4 albedo approximation;
*         =5 Ecco-type isotropic streaming model;
*         =16,26,36,46,56,66 for heterogeneous method with pijk and
*         fixed b (16) or search b_x(26), b_y (36), b_z (46),
*         b_r (56) or b_x=b_y=b_z (66).
* B2      initial or imposed directional bucklings.
* REFKEF  target effective multiplication factor (K-eff).
* IPRINT  print selection for FLU: module (= 0/1/2/3 no print/short
*         print/long print).
* INITFL  flux initialisation flags:
*         = 0 flux initialisation (=0.0 or 1.0);
*         = 1 flux read on LCM;
*         = 2 initialization from DSA flux read on LCM.
* NMERG   number of leakage zones.
* IMERG   leakage zone index in each material mixture zone.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPFLUX,IPMACR
      INTEGER     ITYPEC,MAXOUT,MAXINR,IREBAL,IFRITR,IACITR,ILEAK,
     >            NGROUP,NREGIO,NMAT,NIFISS,ITPIJ,IPRINT,INITFL,NMERG,
     >            IMERG(NMAT)
      REAL        EPSOUT,EPSUNK,EPSINR,B2(4)
      CHARACTER   COPTIO*4
      LOGICAL     LEAKSW,REC
      DOUBLE PRECISION REFKEF
*----
*  LOCAL VARIABLES
*----
      PARAMETER  (NBUCKN=7,NLEAK=6,NSDIR=5,ROUMIN=1.0E-5,RINMIN=5.0E-5,
     >            NSTATE=40)
      CHARACTER   CARLIR*4,CTYPEC*1,CBUCKN(0:NBUCKN)*4,CLEAK(NLEAK)*4,
     >            CSDIR(NSDIR)*1
      SAVE        CBUCKN,CLEAK,CSDIR
      INTEGER     ITYPLU,INTLIR,IRSDIR(NSDIR),ISTATE(NSTATE)
      REAL        REALIR,BSDIR(NSDIR),EPSCON(5)
      DOUBLE PRECISION DBLINP
      DATA       (CBUCKN(JJ),JJ=0,NBUCKN)
     >           /'LKRD','RHS','B0','P0','B1','P1','B0TR','P0TR'/
      DATA       (CLEAK(JJ),JJ=1,NLEAK)
     >           /'PNLR','PNL ','SIGS','ALBS','ECCO','HETE'/
      DATA       (CSDIR(III),III=1,NSDIR)
     >           /'X','Y','Z','R','G'/
*----
*  INITIALIZE TO DEFAULT VALUE
*----
      ISDIR=0
      COPTIO='B0'
      IF((ITPIJ.EQ.1).OR.(ITPIJ.GE.3)) ISDIR=5
      DO 10 III=1,NSDIR
        BSDIR(III)=0.0
        IRSDIR(III)=0
 10   CONTINUE
      REFKEF=1.0D0
      IPRINT=1
      IF(REC) THEN
         CALL LCMGET(IPFLUX,'STATE-VECTOR',ISTATE)
         ITYPEC=ISTATE(6)
         ILEAK=ISTATE(7)
         IFRITR=ISTATE(8)
         IACITR=ISTATE(9)
         IREBAL=ISTATE(10)
         MAXOUT=ISTATE(12)
         NMERG=ISTATE(18)
         CALL LCMGET(IPFLUX,'EPS-CONVERGE',EPSCON)
         EPSINR=EPSCON(1)
         EPSUNK=EPSCON(2)
         EPSOUT=EPSCON(3)
         INITFL=1
         CALL LCMLEN(IPFLUX,'B2  HETE',ILCML1,ITYLCM)
         CALL LCMLEN(IPFLUX,'B2  B1HOM',ILCML2,ITYLCM)
         IF(ILCML1.EQ.3) THEN
           CALL LCMGET(IPFLUX,'B2  HETE',BSDIR)
           IRSDIR(1)=1
           IRSDIR(2)=1
           IRSDIR(3)=1
         ELSE IF(ILCML2.EQ.1) THEN
           CALL LCMGET(IPFLUX,'B2  B1HOM',BSDIR(5))
           IRSDIR(5)=1
         ENDIF
         IF(NMERG.GT.0) CALL LCMGET(IPFLUX,'IMERGE-LEAK',IMERG)
      ELSE
         ITYPEC=-99
         ILEAK=0
         IFRITR=3
         IACITR=3
         IREBAL=1
         MAXOUT=0
         EPSINR=RINMIN
         EPSUNK=0.0
         EPSOUT=ROUMIN
         INITFL=0
         NMERG=1
         IMERG(:NMAT)=1
         CALL LCMPUT(IPFLUX,'IMERGE-LEAK',NMAT,1,IMERG)
      ENDIF
      IF(NGROUP.EQ.1) MAXINR=1
      IF(MOD(ITPIJ,2).EQ.0) THEN
        MAXINR=4*NGROUP
      ELSE
        MAXINR=2*NGROUP
      ENDIF
*----
*  READ OPTION NAME
*----
  20  CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLINP)
      IF(ITYPLU.EQ.10) GO TO 140
  30  IF(ITYPLU.NE.3)
     >  CALL XABORT('FLUGPI: READ ERROR - CHARACTER VARIABLE EXPECTED')
      IF(CARLIR.EQ.';') THEN
        GO TO 140
      ELSE IF(CARLIR.EQ.'EDIT') THEN
        CALL REDGET(ITYPLU,IPRINT,REALIR,CARLIR,DBLINP)
        IF(ITYPLU.NE.1) CALL XABORT('FLUGPI: READ ERROR - INTEGER VA'
     >    //'RIABLE EXPECTED')
      ELSE IF(CARLIR.EQ.'TYPE') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CTYPEC,DBLINP)
        IF(ITYPLU.NE.3) CALL XABORT('FLUGPI: READ ERROR - CHARACTER '
     >    //'VARIABLE EXPECTED')
        IF(CTYPEC.EQ.'N') THEN
          ITYPEC=-1
        ELSE IF(CTYPEC.EQ.'S') THEN
          ITYPEC=0
        ELSE IF(CTYPEC.EQ.'P') THEN
          ITYPEC=1
        ELSE IF(CTYPEC.EQ.'K') THEN
          ITYPEC=2
        ELSE IF(CTYPEC.EQ.'B') THEN
          ITYPEC=4
          ILEAK=3
        ELSE IF(CTYPEC.EQ.'L') THEN
          ITYPEC=5
          ILEAK=3
        ELSE IF(CTYPEC.EQ.'F') THEN
          ITYPEC=-2
          ILEAK=0
          MAXOUT=1
          MAXINR=1
        ELSE
          CALL XABORT('FLUGPI: READ ERROR - INVALID TYPE KEYWORD= '
     >      //CTYPEC//' -- ONLY VALUES ALLOWED ARE: N,S,K,B,L OR F')
        ENDIF
        IF(ITYPEC.GE.2) THEN
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLINP)
          IF(ITYPLU.NE.3) CALL XABORT('FLUGPI: READ ERROR - CHARACTE'
     >      //'R VARIABLE EXPECTED')
          DO 40 JBUC=0,NBUCKN
            IF(CARLIR.EQ.CBUCKN(JBUC)) THEN
              COPTIO=CARLIR
              GO TO 50
            ENDIF
  40      CONTINUE
          GO TO 30
  50      IF(ITYPEC.EQ.2) ITYPEC=3
          ILEAK=3
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLINP)
          IF(ITYPLU.NE.3) CALL XABORT('FLUGPI: READ ERROR - CHARACTE'
     >      //'R VARIABLE EXPECTED')
          DO 70 JLEA=1,NLEAK
            IF(CARLIR.EQ.CLEAK(JLEA)) THEN
              ILEAK=JLEA
              IF(ILEAK.EQ.6) THEN
                CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLINP)
                IF(ITYPLU.NE.3) CALL XABORT('FLUGPI: READ ERROR - '
     >                         //'CHARACTER VARIABLE EXPECTED')
                DO 60 III=1,NSDIR
                  IF(CARLIR.EQ.CSDIR(III)) THEN
                    ISDIR=III
                    GO TO 20
                  ENDIF
  60            CONTINUE
                GO TO 30
              ENDIF
              CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLINP)
              GO TO 80
            ENDIF
  70      CONTINUE
  80      IF((ILEAK.EQ.3).AND.(COPTIO.EQ.'LKRD').AND.LEAKSW) THEN
*           recover the number of leakage zones from the macrolib.
            CALL LCMGET(IPMACR,'STATE-VECTOR',ISTATE)
            NMERG=ISTATE(2)
            DO IBM=1,NMAT
              IMERG(IBM)=IBM
            ENDDO
          ENDIF
          GO TO 30
        ENDIF
      ELSE IF(CARLIR.EQ.'REBA') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLINP)
        IF(ITYPLU.NE.3) CALL XABORT('FLUGPI: READ ERROR - CHARACTE'
     >    //'R VARIABLE EXPECTED')
        IF(CARLIR.EQ.'OFF ') THEN
          IREBAL=0
        ELSE
          IREBAL=1
          GO TO 30
        ENDIF
      ELSE IF(CARLIR.EQ.'INIT') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLINP)
        IF(ITYPLU.NE.3) CALL XABORT('FLUGPI: READ ERROR - CHARACTE'
     >  //'R VARIABLE EXPECTED')
        IF(CARLIR.EQ.'OFF') THEN
*         initial flat distribution
          INITFL=0
        ELSE IF(CARLIR.EQ.'ON') THEN
*         use LCM flux
          INITFL=1
        ELSE IF(CARLIR.EQ.'DSA') THEN
*         use DSA flux
          INITFL=2
        ELSE
          CALL XABORT('FLUGPI: OFF/ON/DSA KEYWORD EXPECTED')
        ENDIF
      ELSE IF(CARLIR.EQ.'EXTE') THEN
        DO 90 II=1,3
          CALL REDGET(ITYPLU,NITMA,EPSOUT,CARLIR,DBLINP)
          IF(ITYPLU.EQ.1) MAXOUT=NITMA
          IF(ITYPLU.EQ.3) GO TO 30
  90    CONTINUE
      ELSE IF(CARLIR.EQ.'UNKT') THEN
        CALL REDGET(ITYPLU,INTLIR,EPSUNK,CARLIR,DBLINP)
        IF(ITYPLU.NE.2)
     >     CALL XABORT('FLUGPI: REAL VALUE OF EPSUNK MUST FOLLOW'
     >               //' UNKT')
        GO TO 20
      ELSE IF(CARLIR.EQ.'THER') THEN
        DO 100 II=1,3
          CALL REDGET(ITYPLU,NITMA,EPSINR,CARLIR,DBLINP)
          IF(ITYPLU.EQ.1) MAXINR=NITMA
          IF(ITYPLU.EQ.3) GO TO 30
 100    CONTINUE
      ELSE IF(CARLIR.EQ.'ACCE') THEN
        CALL REDGET(ITYPLU,IFRITR,REALIR,CARLIR,DBLINP)
        IF(ITYPLU.NE.1) CALL XABORT('FLUGPI: READ ERROR - INTEGER VA'
     >    //'RIABLE EXPECTED')
        CALL REDGET(ITYPLU,IACITR,REALIR,CARLIR,DBLINP)
        IF(ITYPLU.NE.1) CALL XABORT('FLUGPI: READ ERROR - INTEGER VA'
     >    //'RIABLE EXPECTED')
      ELSE IF(CARLIR.EQ.'KEFF') THEN
         CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLINP)
         IF(ITYPLU.NE.2) CALL XABORT('FLUGPI: READ ERROR - REAL VA'
     >   //'RIABLE EXPECTED FOLLOWING KEFF KEYWORD')
         REFKEF=REALIR
      ELSE IF(CARLIR.EQ.'BUCK') THEN
         CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLINP)
         CALL XDISET(IRSDIR,NSDIR,0)
         IF(ITYPLU.EQ.2) THEN
           BSDIR(5)=REALIR
           IRSDIR(5)=1
           GO TO 20
         ELSE IF(ITYPLU.EQ.1) THEN
           CALL XABORT('FLUGPI: READ ERROR - INTEGER '
     >       //'VARIABLE FOUND FOLLOWING BUCK KEYWORD')
         ENDIF
 110     CONTINUE
         DO 120 III=1,NSDIR
           IF(CARLIR.EQ.CSDIR(III)) THEN
             CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLINP)
             IF(ITYPLU.NE.2)
     >         CALL XABORT('FLUGPI: READ ERROR - REAL VARIABLE '//
     >           'EXPECTED FOLLOWING BUCKLING DIRECTION KEYWORD')
             BSDIR(III)=REALIR
             IRSDIR(III)=1
             GO TO 130
           ENDIF
 120     CONTINUE
         GO TO 30
 130     CONTINUE
         CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLINP)
         IF(ITYPLU.NE.3)
     >     CALL XABORT('FLUGPI: READ ERROR - CHARACTER '
     >     //'VARIABLE EXPECTED')
         GO TO 110
      ELSE IF(CARLIR.EQ.'IDEM') THEN
         CALL XDISET(IRSDIR,NSDIR,0)
         CALL LCMLEN(IPMACR,'B2  HETE',ILCMLN,ITYLCM)
         IF(ILCMLN.EQ.3) THEN
            CALL LCMGET(IPMACR,'B2  HETE',BSDIR)
            IRSDIR(1)=1
            IRSDIR(2)=1
            IRSDIR(3)=1
         ENDIF
         CALL LCMLEN(IPMACR,'B2  B1HOM',ILCMLN,ITYLCM)
         IF(ILCMLN.EQ.1) THEN
           CALL LCMGET(IPMACR,'B2  B1HOM',BSDIR(5))
           IRSDIR(5)=1
         ENDIF
      ELSE IF(CARLIR.EQ.'MERG') THEN
         DO IBM=1,NMAT
           CALL REDGET(ITYPLU,IMERG(IBM),REALIR,CARLIR,DBLINP)
           IF(ITYPLU.NE.1) CALL XABORT('FLUGPI: READ ERROR - INTEGER V'
     >     //'ARIABLE EXPECTED')
           NMERG=MAX(NMERG,IMERG(IBM))
         ENDDO
         CALL LCMPUT(IPFLUX,'IMERGE-LEAK',NMAT,1,IMERG)
      ELSE
         CALL XABORT('FLUGPI: READ ERROR - ILLEGAL KEYWORD '//CARLIR)
      ENDIF
      GO TO 20
*----
 140  CONTINUE
      IF(ITYPEC.EQ.3) THEN
        ISDIR=0
        IF(ILEAK.EQ.6) THEN
          DO 150 III=1,NSDIR
            IF(IRSDIR(III).EQ.1) GO TO 160
 150      CONTINUE
        ELSE
          GO TO 160
        ENDIF
        CALL XABORT('FLUGPI: NO BUCKLING READ FOR TYPE K '//
     >    'CALCULATION WITH IMPOSED BUCKLING')
 160    CONTINUE
      ENDIF
      IF(ILEAK.EQ.6) THEN
        IF(IRSDIR(5).EQ.1) THEN
          DO 210 III=1,NSDIR-1
            IF(IRSDIR(III).EQ.1) THEN
              CALL XABORT('FLUGPI: GLOBAL INITIAL BUCKLING '//
     >        'INCONSISTENT WITH X, Y, Z, R BUCKLING')
            ENDIF
 210      CONTINUE
          B2(1)=BSDIR(5)/3.0
          B2(2)=B2(1)
          B2(3)=B2(1)
          B2(4)=BSDIR(5)
        ELSE IF(IRSDIR(4).EQ.1) THEN
          DO 220 III=1,2
            IF(IRSDIR(III).EQ.1) THEN
              CALL XABORT('FLUGPI: RADIAL INITIAL BUCKLING '//
     >        'INCONSISTENT WITH X, Y BUCKLING')
            ENDIF
 220      CONTINUE
          B2(1)=BSDIR(4)/2.0
          B2(2)=B2(1)
          B2(3)=BSDIR(3)
          B2(4)=BSDIR(3)+BSDIR(4)
        ELSE
          B2(1)=BSDIR(1)
          B2(2)=BSDIR(2)
          B2(3)=BSDIR(3)
          B2(4)=BSDIR(1)+BSDIR(2)+BSDIR(3)
        ENDIF
        ILEAK=(ISDIR+1)*10+ILEAK
      ELSE
        IF(IRSDIR(4).NE.0.OR.IRSDIR(3).NE.0.OR.
     >     IRSDIR(2).NE.0.OR.IRSDIR(1).NE.0)
     >    CALL XABORT('FLUGPI: FOR HOMOGENEOUS LEAKAGE METHOD'//
     >    'DIRECTIONS X ,Y, Z AND R BUCKLING ARE ILLEGAL')
        B2(4)=BSDIR(5)
      ENDIF
      IF(EPSOUT.LT.1.0E-10) THEN
        CALL XABORT('FLUGPI: ERROR -- EPSOUT MUST BE GREATER '//
     >                       'THAN 1.0E-10')
      ENDIF
      IF(EPSINR.LT.1.0E-10) THEN
        CALL XABORT('FLUGPI: ERROR -- EPSINR MUST BE GREATER '//
     >                       'THAN 1.0E-10')
      ENDIF
      IF(EPSUNK.LE.1.0E-10) EPSUNK=EPSINR
      IF(ITYPEC.EQ.-99) CALL XABORT('FLUGPI: TYPE NOT DEFINED.')
      IF(MAXOUT.EQ.0)THEN
        IF(ITYPEC.LE.2) THEN
          MAXOUT=MAX(2*NREGIO-1,2*NIFISS-1)
        ELSE
          MAXOUT=MAX(10*NREGIO+1,10*NIFISS+1)
        ENDIF
      ENDIF
      RETURN
      END
