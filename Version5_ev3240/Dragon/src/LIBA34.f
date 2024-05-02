*DECK LIBA34
      SUBROUTINE LIBA34(HNAMIS,NGRO,FGHOMO,NGHOMO,NSEQHO,NTEMPS,LFIS,
     1 L104,SEQHOM,TEMPS,TN,SN,ABSOHE,DIFFHE,FISSHE,FLUXHE,IMPX,TAUX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Temperature and dilution interpolation of self-shielded effective
* rates in the APOLIB-3 format.
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
* HNAMIS  name of the isotope.
* NGRO    number of energy groups.
* FGHOMO  first self-shielded energy group.
* NGHOMO  number of self-shielded energy groups.
* NSEQHO  number of tabulated dilutions.
* NTEMPS  number of tabulated temperatures.
* LFIS    fission reaction flag (=.true. if the fission reaction is
*         self-shielded).
* L104    resonance flux flag (=.true. if the apolib contains dilution
*         /temperature-dependent flux information). If this information
*         is not provided, it will be reconstructed from a balance
*         relation.
* SEQHOM  tabulated dilutions.
* TEMPS   tabulated temperatures.
* TN      temperature of isotope.
* SN      dilution of isotope.
* ABSOHE  tabulated absorption effective reaction rates.
* DIFFHE  tabulated diffusion effective reaction rates.
* FISSHE  tabulated nu*fission effective reaction rates 
*         (if LFIS=.true.).
* FLUXHE  tabulated self-shielded fluxes (if L104=.true.).
* IMPX    print flag.
*
*Parameters: output
* TAUX    interpolated effective rates:
*         TAUX(I,1) absorption effective rates;
*         TAUX(I,2) diffusion effective rates;
*         TAUX(I,3) nu*fission effective rates;
*         TAUX(I,4) pseudo-absorption effective rates used to
*                    reconstruct the self-shielded flux;
*         TAUX(I,5) infinite-dilution absorption x-s;
*         TAUX(I,6) infinite-dilution diffusion x-s;
*         TAUX(I,7) infinite-dilution fission x-s.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER HNAMIS*12
      INTEGER NGRO,FGHOMO,NSEQHO,NGHOMO,NTEMPS,IMPX
      LOGICAL LFIS,L104
      REAL SEQHOM(NSEQHO),TEMPS(NTEMPS),TN,SN(NGRO),
     1 ABSOHE(NSEQHO,NGHOMO,NTEMPS),DIFFHE(NSEQHO,NGHOMO,NTEMPS),
     2 FISSHE(NSEQHO,NGHOMO,NTEMPS),FLUXHE(NSEQHO,NGHOMO,NTEMPS),
     3 TAUX(NGHOMO,7)
*----
*  LOCAL VARIABLES
*----
      CHARACTER HSMG*131,TEXTE*80
      PARAMETER (NINT=2,NINTSS=3,DTMIN=1.0)
      LOGICAL LGONE
      DOUBLE PRECISION S1,S2,S3,S4,SUMA,SUMS,SUMF,SUM104,REL,RNTERP
      REAL, ALLOCATABLE, DIMENSION(:,:) :: ABSOH,DIFFH,FISSH,FLUXH
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: SQRTEM,SEQ2,WEIJHT,
     1 WEIGH
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(WEIJHT(NTEMPS),SQRTEM(NTEMPS),ABSOH(NSEQHO,NGHOMO),
     1 DIFFH(NSEQHO,NGHOMO),FISSH(NSEQHO,NGHOMO),
     2 FLUXH(NSEQHO,NGHOMO))
      CALL XDDSET(WEIJHT,NTEMPS,0.0D0)
*----
*  SQUARE ROOT OF TEMPERATURE INTERPOLATION.
*
*  IGTFIX=1 IF ONLY ONE TABULATED TEMPERATURE OR IF STT IS ONE OF THE
*  TABULATED TEMPERATURES. IGTFIX=2 IF STT IS OUTSIDE THE TABULATED
*  RANGE. IGTFIX=0 OTHERWISE.
*
*----
      DO 10 I=1,NTEMPS
      SQRTEM(I)=SQRT(TEMPS(I))
   10 CONTINUE
      IF(NTEMPS.EQ.1) THEN
        IGTFIX=1
        IPROX=1
      ELSE
        STT=SQRT(TN)
        CALL LIBA28(STT,SQRTEM,NTEMPS,NINT,WEIJHT,IORD,IPROX,I0)
        IF(ABS(TN-TEMPS(IPROX)).LE.DTMIN) THEN
          IGTFIX=1
        ELSEIF((STT.LT.SQRTEM(1)).OR.(STT.GT.SQRTEM(NTEMPS))) THEN
          WRITE(HSMG,'(A,F8.2,A,F8.2,A,F8.2,2A)')
     1    'LIBA34: A TEMPERATURE', TN,'K IS NOT INCLUDED BETWEEN ',
     2    TEMPS(1),' AND ',TEMPS(NTEMPS),' ISOTOPE:',HNAMIS
          WRITE(6,'(/1X,A)') HSMG
          IGTFIX=2
        ELSE
          IGTFIX=0
        ENDIF
      ENDIF
*
      IF(IGTFIX .EQ. 1) THEN
         DO 25 I=1,NGHOMO
         DO 20 J=1,NSEQHO
         ABSOH(J,I)=ABSOHE(J,I,IPROX)
         DIFFH(J,I)=DIFFHE(J,I,IPROX)
         IF(LFIS) FISSH(J,I)=FISSHE(J,I,IPROX)
         IF(L104) FLUXH(J,I)=FLUXHE(J,I,IPROX)
   20    CONTINUE
   25    CONTINUE
      ELSE
       DO 45 J=1,NSEQHO
       DO 40 I=1,NGHOMO
        S1=0.D0
        S2=0.D0
        S3=0.D0
        S4=0.D0
        DO 30 K=1,IORD
         S1=S1+WEIJHT(K)*ABSOHE(J,I,I0+K)
         S2=S2+WEIJHT(K)*DIFFHE(J,I,I0+K)
         IF(LFIS)S3=S3+WEIJHT(K)*FISSHE(J,I,I0+K)
         IF(L104)S4=S4+WEIJHT(K)*FLUXHE(J,I,I0+K)
   30   CONTINUE
        IF(IGTFIX.EQ.2) THEN
          IF(ABSOHE(J,I,IPROX).GE.0.) THEN
            S1=MAX(0.D0,S1)
          ELSE
            S1=MIN(S1,0.D0)
          ENDIF
          IF(DIFFHE(J,I,IPROX).GE.0.) THEN
            S2=MAX(0.D0,S2)
          ELSE
            S2=MIN(S2,0.D0)
          ENDIF
        ENDIF
        ABSOH(J,I)=REAL(S1)
        DIFFH(J,I)=REAL(S2)
        IF(LFIS) THEN
          IF(IGTFIX .EQ. 2) THEN
            IF(FISSHE(J,I,IPROX).GE.0.) THEN
              S3=MAX(0.D0,S3)
            ELSE
              S3=MIN(S3,0.D0)
            ENDIF
          ENDIF
          FISSH(J,I)=REAL(S3)
        ENDIF
        IF(L104) THEN
          IF(IGTFIX .EQ. 2) THEN
            IF(FLUXHE(J,I,IPROX).GE.0.) THEN
              S4=MAX(0.D0,S4)
            ELSE
              S4=MIN(S4,0.D0)
            ENDIF
          ENDIF
          FLUXH(J,I)=REAL(S4)
        ENDIF
   40  CONTINUE
   45  CONTINUE
      ENDIF
*----
*  SET INFINITE DILUTION VALUES.
*----
      DO 50 I=1,NGHOMO
      TAUX(I,5)=ABSOH(NSEQHO,I)
      TAUX(I,6)=DIFFH(NSEQHO,I)
      IF(LFIS) TAUX(I,7)=FISSH(NSEQHO,I)
   50 CONTINUE
*----
*  DILUTION INTERPOLATION.
*----
      LGONE=NSEQHO.EQ.1
      NSEQH1=0
      SEQHO1=0.0
      SEQHO0=0.0
      IF(.NOT.LGONE)THEN
        NSEQH1=NSEQHO-1
        SEQHO1=SEQHOM(NSEQH1)
        SEQHO0=SEQHOM(NSEQHO)
      ENDIF
      DO 110 IGG=FGHOMO,FGHOMO+NGHOMO-1
         IGSSC=IGG+1-FGHOMO
         BACK=SN(IGG)
         IF(LGONE) THEN
*----
*  UNIQUE TABULATED TEMPERATURE.
*----
             TAUX(IGSSC,1)=ABSOH(NSEQHO,IGSSC)
             TAUX(IGSSC,2)=DIFFH(NSEQHO,IGSSC)
             IF(LFIS) TAUX(IGSSC,3)=FISSH(NSEQHO,IGSSC)
             IF(L104) TAUX(IGSSC,4)=FLUXH(NSEQHO,IGSSC)
           GOTO 110
         ENDIF
*----
*  MANY TABULATED TEMPERATURES.
*----
         IF(BACK.GE.SEQHO1)THEN
*
*           ASYMPTOTIC BEHAVIOR: REACTION RATES VARY LINEARLY WITH
*           1/SEQHOM FOR THE LAST 2 POINTS OF THE TABULATION
*
            IF(BACK.GT.SEQHO0) BACK=SEQHO0
            AUX=1.0/(BACK*(SEQHO0-SEQHO1))
            AUX1=SEQHO1*(SEQHO0-BACK)*AUX
            AUX2=SEQHO0*(SEQHO1-BACK)*AUX
            TAUX(IGSSC,1)=ABSOH(NSEQH1,IGSSC)*AUX1
     1                     -ABSOH(NSEQHO,IGSSC)*AUX2
            TAUX(IGSSC,2)=DIFFH(NSEQH1,IGSSC)*AUX1
     1                     -DIFFH(NSEQHO,IGSSC)*AUX2
            IF(LFIS) TAUX(IGSSC,3)=FISSH(NSEQH1,IGSSC)*AUX1
     1               -FISSH(NSEQHO,IGSSC)*AUX2
            IF(L104) TAUX(IGSSC,4)=FLUXH(NSEQH1,IGSSC)*AUX1
     1               -FLUXH(NSEQHO,IGSSC)*AUX2
         ELSE
*
*           REACTION RATES VARY WITH THE SQRT OF THE BACKGROUND XSECT
*
            BACKH2=SQRT(BACK)
            ALLOCATE(SEQ2(NSEQHO),WEIGH(NINTSS))
            DO 60 I=1,NSEQHO
            SEQ2(I)=SQRT(SEQHOM(I))
   60       CONTINUE
            CALL LIBA28(BACKH2,SEQ2,NSEQHO,NINTSS,WEIGH,IORD,IPR,I0)
            DO 70 ISEQHO=1,NSEQHO
              IF(ABS(BACK-SEQHOM(ISEQHO)).LE.1.E-2) THEN
                TAUX(IGSSC,1)=ABSOH(ISEQHO,IGSSC)
                TAUX(IGSSC,2)=DIFFH(ISEQHO,IGSSC)
                IF(LFIS) TAUX(IGSSC,3)=FISSH(ISEQHO,IGSSC)
                IF(L104) TAUX(IGSSC,4)=FLUXH(ISEQHO,IGSSC)
                DEALLOCATE(WEIGH,SEQ2)
                GOTO 110
              ENDIF
   70       CONTINUE
            SUMA=0.D0
            SUMS=0.D0
            SUMF=0.D0
            SUM104=0.D0
            DO 80 I=1,IORD
            I1=I+I0
            SUMA=SUMA+WEIGH(I)*ABSOH(I1,IGSSC)
            SUMS=SUMS+WEIGH(I)*DIFFH(I1,IGSSC)
            IF(LFIS) SUMF=SUMF+WEIGH(I)*FISSH(I1,IGSSC)
            IF(L104) SUM104=SUM104+WEIGH(I)*FLUXH(I1,IGSSC)
   80       CONTINUE
            DO 90 I=1,IORD
              I1=I+I0
              IF(SEQHOM(I1).GT.BACK) THEN
                IF(I1-1.GT.0) THEN
*
*               ABSORPTION RATE CRITERION.
*
                YMIN=MIN(ABSOH(I1-1,IGSSC),ABSOH(I1,IGSSC))
                YMAX=MAX(ABSOH(I1-1,IGSSC),ABSOH(I1,IGSSC))
                IF((SUMA.GT.YMAX) .OR. (SUMA.LT.YMIN)) THEN
                  RNTERP=SUMA
                  SUMA=ABSOH(I1-1,IGSSC)+
     1            (ABSOH(I1,IGSSC)-ABSOH(I1-1,IGSSC))*
     1            (BACKH2-SEQ2(I1-1))/(SEQ2(I1)-SEQ2(I1-1))
                  REL = (RNTERP-SUMA)/SUMA
                  IF(REL.GE.0.1 .OR. IMPX .GT. 3) THEN
                    WRITE(TEXTE,10000)
     1              'ABS. G=',IGG,' DIL=',BACK,
     1              ' INT. LIN. --> ERR. RELA.=',REL
                      WRITE(6,'(/1X,A)') TEXTE
                    ENDIF
                  ENDIF
*
*               SCATTERING RATE CRITERION.
*
                YMIN = MIN(DIFFH(I1-1,IGSSC),DIFFH(I1,IGSSC))
                YMAX = MAX(DIFFH(I1-1,IGSSC),DIFFH(I1,IGSSC))
                IF((SUMS.GT.YMAX) .OR. (SUMS.LT.YMIN)) THEN
                  RNTERP=SUMS
                  SUMS=DIFFH(I1-1,IGSSC)+
     1                 (DIFFH(I1,IGSSC)-DIFFH(I1-1,IGSSC))*
     1                 (BACKH2-SEQ2(I1-1))/(SEQ2(I1)-SEQ2(I1-1))
                  REL = (RNTERP-SUMS)/SUMS
                  IF(REL.GE. 0.1 .OR. IMPX .GT. 3) THEN
                    WRITE(TEXTE,10000)
     1              'DIF. G=',IGG,' DIL=',BACK,
     1              ' INT. LIN. --> ERR. RELA.=',REL
                      WRITE(6,'(/1X,A)') TEXTE
                   ENDIF
                  ENDIF
*
*               PRODUCTION RATE CRITERION.
*
                IF(LFIS) THEN
                  YMIN = MIN(FISSH(I1-1,IGSSC),FISSH(I1,IGSSC))
                  YMAX = MAX(FISSH(I1-1,IGSSC),FISSH(I1,IGSSC))
                  IF((SUMF.GT.YMAX) .OR. (SUMF.LT.YMIN)) THEN
                    RNTERP=SUMF
                    SUMF=FISSH(I1-1,IGSSC)+
     1                   (FISSH(I1,IGSSC)-FISSH(I1-1,IGSSC))*
     1                   (BACKH2-SEQ2(I1-1))/(SEQ2(I1)-SEQ2(I1-1))
                    REL = (RNTERP-SUMF)/SUMF
                    IF(REL.GE.0.1 .OR. IMPX .GT. 3) THEN
                      WRITE(TEXTE,10000)
     1                'FIS. G=',IGG,' DIL=',BACK,
     1                ' INT. LIN. --> ERR. RELA.=',REL
                      WRITE(6,'(/1X,A)') TEXTE
                      ENDIF
                    ENDIF
                  ENDIF
*
*               TEST FLUX 104
*
                IF(L104) THEN
                  YMIN = MIN(FLUXH(I1-1,IGSSC),FLUXH(I1,IGSSC))
                  YMAX = MAX(FLUXH(I1-1,IGSSC),FLUXH(I1,IGSSC))
                  IF((SUM104.GT.YMAX) .OR. (SUM104.LT.YMIN)) THEN
                    RNTERP=SUM104
                    SUM104=FLUXH(I1-1,IGSSC)+
     1                   (FLUXH(I1,IGSSC)-FLUXH(I1-1,IGSSC))*
     1                   (BACKH2-SEQ2(I1-1))/(SEQ2(I1)-SEQ2(I1-1))
                    REL = (RNTERP-SUM104)/SUM104
                    IF(REL.GE.0.1 .OR. IMPX .GT. 3) THEN
                      WRITE(TEXTE,10000)
     1                'FIS. G=',IGG,' DIL=',BACK,
     1                ' INT. LIN. --> ERR. RELA.=',REL
                      WRITE(6,'(/1X,A)') TEXTE
                      ENDIF
                    ENDIF
                  ENDIF
*
                ELSE
                  SUMA=ABSOH(1,IGSSC)+
     1                 (ABSOH(2,IGSSC)-ABSOH(1,IGSSC))*
     1                 (BACKH2-SEQ2(1))/(SEQ2(2)-SEQ2(1))
                  IF(SUMA.LE.0.) THEN
                    SUMA=ABSOH(1,IGSSC)
                    WRITE(TEXTE,3000)
     1              ' DIL. : ',BACK, ' TROP PETITE ',
     2              'TAUX ABS. NON EXTRAPOLES GR. ',IGG
                    WRITE(6,'(/1X,A)') TEXTE
                    ENDIF
*
                  SUMS=DIFFH(1,IGSSC)+
     1                 (DIFFH(2,IGSSC)-DIFFH(1,IGSSC))*
     1                 (BACKH2-SEQ2(1))/(SEQ2(2)-SEQ2(1))
                  IF(SUMS.LE.0.) THEN
                    SUMS=DIFFH(1,IGSSC)
                    WRITE(TEXTE,3000)
     1              ' DIL. : ',BACK, ' TROP PETITE ',
     2              'TAUX DIFF. NON EXTRAPOLES GR. ',IGG
                    WRITE(6,'(/1X,A)') TEXTE
                    ENDIF
*
                  IF(LFIS) SUMF=FISSH(1,IGSSC)+
     1                 (FISSH(2,IGSSC)-FISSH(1,IGSSC))*
     1                 (BACKH2-SEQ2(1))/(SEQ2(2)-SEQ2(1))
                  IF(LFIS.AND.SUMF.LE.0.) THEN
                    SUMF=FISSH(1,IGSSC)
                    WRITE(TEXTE,3000)
     1              ' DIL. : ',BACK, ' TROP PETITE ',
     2              'TAUX PROD. NON EXTRAPOLES GR. ',IGG
                    WRITE(6,'(/1X,A)') TEXTE
                    ENDIF
*
                  IF(L104) SUM104=FLUXH(1,IGSSC)+
     1                 (FLUXH(2,IGSSC)-FLUXH(1,IGSSC))*
     1                 (BACKH2-SEQ2(1))/(SEQ2(2)-SEQ2(1))
                  IF(L104.AND.SUM104.LE.0.) THEN
                    SUM104=FLUXH(1,IGSSC)
                    WRITE(TEXTE,3000)
     1              ' DIL. : ',BACK, ' TROP PETITE ',
     2              'FLUX 104 NON EXTRAPOLES GR. ',IGG
                    WRITE(6,'(/1X,A)') TEXTE
                    ENDIF
                  ENDIF
                GOTO 100
              ENDIF
   90       CONTINUE
*
  100       TAUX(IGSSC,1)=REAL(SUMA)
            TAUX(IGSSC,2)=REAL(SUMS)
            IF(LFIS) TAUX(IGSSC,3)=REAL(SUMF)
            IF(L104) TAUX(IGSSC,4)=REAL(SUM104)
            IF(SUMA.LE.0.) THEN
              WRITE(TEXTE,1000)
     1        HNAMIS,'GROUPE ',IGG,' DIL. ',BACK,' ABS. <= 0.'
              CALL XABORT('LIBA34:'//TEXTE)
            ENDIF
            IF(SUMS.LE.0.) THEN
              WRITE(TEXTE,1000)
     1        HNAMIS,'GROUPE ',IGG,' DIL. ',BACK,' DIF. <= 0.'
              CALL XABORT('LIBA34:'//TEXTE)
            ENDIF
            IF(LFIS .AND. SUMF.LE.0.) THEN
              WRITE(TEXTE,1000)
     1        HNAMIS,'GROUPE ',IGG,' DIL. ',BACK,' FIS. <= 0.'
              CALL XABORT('LIBA34:'//TEXTE)
            ENDIF
            IF(L104 .AND. (1.-SUM104/BACK).LE.0.) THEN
              WRITE(TEXTE,1000)
     1        HNAMIS,'GROUPE ',IGG,' DIL. ',BACK,' FLU. <= 0.'
              CALL XABORT('LIBA34:'//TEXTE)
            ENDIF
            DEALLOCATE(WEIGH,SEQ2)
         ENDIF
  110 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(FLUXH,FISSH,DIFFH,ABSOH,SQRTEM,WEIJHT)
      RETURN
*
1000  FORMAT(9H ISOTOPE:,A12,2X,A,I3,A,E13.5,A)
3000  FORMAT(A,E13.5,A,A,I3)
10000 FORMAT(A,I3,A,1P,E13.5,A,E13.5)
      END
