*DECK CLMGET
      SUBROUTINE CLMGET(IPRINT,NBMIX,NBISO,ISONRF,ISOMIX,
     >                  NCLM  ,IDCLM,IACT ,DENRD )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read CLM module options.
*
*Copyright:
* Copyright (C) 2022 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input
* IPRINT  print index,
* NBMIX   maximum number of mixtures.
* NBISO   maximum number of isotopes.
* ISONRF  reference names of isotopes.
* ISOMIX  mixture associated with each isotope.
*
*Parameters: output
* NCLM    number of liquid mixtures to combine.
* IDCLM   liquid mixtures indices to combine.
* IACT    isotope identifier (IACT(1,ISO)) for mixture considered,
*         reference isotope (IACT(2,ISO)) 
*         and action on each isotope for which concentration 
*         is modified with: IACT(3,ISO)=0 no change; 
*         IACT(3,ISO)=-1 for ADDI ABS; IACT(3,ISO)=1 for ADDI REL;
*         IACT(3,ISO)=-2 for SETI ABS; IACT(3,ISO)=2 for SETI REL.
* DENRD   isotope concentration or relative concentration.
*
*Comments:
* Input data is of the form:
*    [ EDIT iprint ]
*    MIXCLM (IDCLM(ii),ii=1,NCLM) 
*    [ { ADDI | SETI } { ABS | REL } (isot(ii) dens(ii),ii=1,niso)] 
*-----------------------------------------------------------------------
*
      IMPLICIT         NONE
      INTEGER          IPRINT,NBMIX,NBISO,ISONRF(3,NBISO),ISOMIX(NBISO),
     >                 NCLM,IDCLM(NBMIX),IACT(3,NBISO)
      REAL             DENRD(NBISO)
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='CLMGET')
*----
*  REDGET variables
*----
      INTEGER          ITYPLU,INTLIR
      CHARACTER        CARLIR*8
      REAL             REALIR
      DOUBLE PRECISION DBLLIR
*----
*  LOCAL variables
*----
      INTEGER          IMIX,ISO,JSO,JACT,KSO,INAM(2)
*----
*  INITIALIZE MIXMER
*----
      CALL XDISET(IDCLM,NBMIX,0)
      CALL XDISET(IACT,3*NBISO,0)
      CALL XDRSET(DENRD,NBISO,0.0)
*----
*  READ OPTION NAME
*----
 10   CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
 20   IF(ITYPLU.EQ.10) GO TO 100
      IF(ITYPLU.NE.3) CALL XABORT(NAMSBR//': READ ERROR - '//
     >'Character variable expacted')
      IF(CARLIR.EQ.';') GO TO 100
      IF(CARLIR.EQ.'EDIT') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT(NAMSBR//': READ ERROR -'//
     >  'Integer variable expacted')
        IPRINT=INTLIR
      ELSE IF(CARLIR.EQ.'MIXCLM') THEN
        NCLM=0
        DO IMIX=1,NBMIX
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.NE.1) THEN
            DO ISO=1,NBISO
              WRITE(IOUT,'(2A4,2(5X,I10))') ISONRF(1,ISO),ISONRF(2,ISO),
     >        IACT(1,ISO),IACT(2,ISO)
            ENDDO
            GO TO 20
          ENDIF
          NCLM=NCLM+1
          IF(INTLIR .LE. 0 .OR. INTLIR .GT. NBMIX) CALL XABORT(NAMSBR//
     >    ': READ ERROR - Mixture < 0 or > NBMIX')
          IDCLM(NCLM)=INTLIR
*----
*  Associate isotopes mixture number to first mixture to process
*---- 
          IF(NCLM.EQ.1) THEN
            DO ISO=1,NBISO
              IF(ISOMIX(ISO).EQ.INTLIR) THEN
                IACT(1,ISO)=NCLM
              ENDIF
            ENDDO
          ELSE
*----
*  Test additional mixture number for coherent isotopic contents
*---- 
            DO ISO=1,NBISO
              IF(ISOMIX(ISO).EQ.INTLIR) THEN
                DO JSO=1,NBISO
                  IF(IACT(1,JSO).EQ.1) THEN
                    IF(ISONRF(1,ISO).EQ.ISONRF(1,JSO) .AND.
     >                 ISONRF(2,ISO).EQ.ISONRF(2,JSO)) THEN
                      IACT(1,ISO)=NCLM
                      IACT(2,ISO)=JSO
                      GO TO 110
                    ENDIF
                  ENDIF
                ENDDO
                CALL XABORT(NAMSBR//
     >          ': Mixtures do not have the same isotopic contents')
 110            CONTINUE
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ELSE IF(CARLIR.EQ.'ADDI' .OR. CARLIR.EQ.'SETI') THEN
        JACT=1
        IF(CARLIR.EQ.'SETI') JACT=2
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.EQ.3) THEN
          IF(CARLIR.EQ.'ABS') THEN
            JACT=-JACT
          ELSE IF (CARLIR.NE.'REL') THEN
            CALL XABORT(NAMSBR//
     >    ': READ ERROR - Invalid ADDI or SETI option.'//
     >    ' Only REL or ABS valid')
          ENDIF
        ELSE
          CALL XABORT(NAMSBR//
     >    ': READ ERROR - No ADDI or SETI option provided')
        ENDIF
*----
*  Read all isotopes for SETI and ADDI.
*----
        DO ISO=1,NBISO
          KSO=0
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.EQ.3) THEN
*----
* Test if valid isotopes associated to MIXCLM
*---- 
            READ(CARLIR,'(2A4)') INAM(1),INAM(2)
            DO JSO=1,NBISO
*----
*  Only need to check first mixture
*----
              IF(IACT(1,JSO).EQ.1) THEN
                IF(INAM(1).EQ.ISONRF(1,JSO) .AND.
     >             INAM(2).EQ.ISONRF(2,JSO)) THEN
                  IACT(3,JSO)=JACT
                  KSO=JSO
                  GO TO 120
                ENDIF
              ENDIF
            ENDDO
            GO TO 20
          ELSE
            CALL XABORT(NAMSBR//
     >    ': READ ERROR - Invalid isotope name')
          ENDIF 
 120      CONTINUE
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.EQ.2) THEN
            DENRD(KSO)=REALIR            
          ELSE
            CALL XABORT(NAMSBR//
     >    ': READ ERROR - Invalid isotopic density (REL or ABS)')
          ENDIF 
        ENDDO
      ELSE
        CALL XABORT(NAMSBR//': READ ERROR - '//
     >  'Illegal keyword')
      ENDIF
      GO TO 10
 100  CONTINUE
*----
*  RETURN
*----
*----
*  Print if required
*----
      DO ISO=1,NBISO
        WRITE(IOUT,'(2A4,3(5X,I5),1P,E20.9)') 
     >        ISONRF(1,ISO),ISONRF(2,ISO),
     >        IACT(1,ISO),IACT(2,ISO),IACT(3,ISO),DENRD(ISO)
      ENDDO
      RETURN
      END
