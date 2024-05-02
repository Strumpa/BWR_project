*DECK CLM
      SUBROUTINE CLM(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Combine liquid fuel mixtures from different fuel channels 
* and redistribute in channels.
*
*Copyright:
* Copyright (C) 2022 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): G. Marleau
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1): create or modification type(L_LIBRARY)
*         HENTRY(2): optional read-only type(L_MACROLIB) used to
*                    initialize a new lattice code library.
* IENTRY  type of each LCM object or file:
*         =1 LCM memory object; =2 XSM file; =3 sequential binary file;
*         =4 sequential ascii file.
* JENTRY  access of each LCM object or file:
*         =0 the LCM object or file is created;
*         =1 the LCM object or file is open for modifications;
*         =2 the LCM object or file is open in read-only mode.
* KENTRY  LCM object address or file unit number.
*
* Comment:
*  All the mixture must contain the same isotopes with possibly
*  different concentrations.
*-----------------------------------------------------------------------
*
      USE          GANLIB
      IMPLICIT     NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER      NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR)  KENTRY(NENTRY)
      CHARACTER    HENTRY(NENTRY)*12
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR)  IPLIB,IPTRK
      INTEGER      IOUT,NSTATE,ILCMUP,ILCMDN
      CHARACTER    NAMSBR*6
      PARAMETER   (IOUT=6,NSTATE=40,ILCMUP=1,ILCMDN=2,NAMSBR='CLM   ')
*----
*  LOCAL PARAMETERS
*----
      INTEGER      NBSL,NBST,IEN,ISTATE(NSTATE),NBMIX,NBISO,NBREG,IREG,
     >             IPRINT,NCLM,ICLM,ISO,JSO,NGRO,ITSTMP
      INTEGER      MIXI,MIXJ
      REAL         VOLTOT,TMPDAY(3)
      CHARACTER    HSIGN*12
      INTEGER, ALLOCATABLE, DIMENSION(:)   :: IDCLM,ISOMIX,MATCOD
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISONRF,IACT
      REAL, ALLOCATABLE, DIMENSION(:)      :: VOLUME,DENISO,VOLMIX,
     >                                        DENRD
      LOGICAL, ALLOCATABLE, DIMENSION(:)   :: MASK,MASKL
*----
*  PARAMETER VALIDATION.
*----
      IF(NENTRY .LT. 2) CALL XABORT(NAMSBR//
     >': At least 2 parameters expected.')
      IF(IENTRY(1) .NE. 1 .AND.
     >   IENTRY(1) .NE. 2) CALL XABORT(NAMSBR//
     >': LCM OBJECT OR XSM FILE EXPECTED AT LHS.')
      IF(JENTRY(1) .NE. 1) CALL XABORT(NAMSBR//': ENTRY'
     1 //' IN MODIFICATION MODE EXPECTED.')
      IPLIB=KENTRY(1)
*----
*  Find IPLIB and IPTRK structures
*----
      IPRINT=1
      NBSL=0
      NBST=0
      DO IEN=1,NENTRY
        HSIGN=' '
        IF(NBSL .EQ. 0) THEN
*----
*  Find Library to modify
*----
          IF(IENTRY(IEN) .LE. 2 .AND. JENTRY(IEN) .EQ. 1) THEN
            CALL LCMGTC(KENTRY(IEN),'SIGNATURE',12,1,HSIGN)
            IF(HSIGN.EQ.'L_LIBRARY') THEN
              IPLIB=KENTRY(IEN)
              WRITE(IOUT,6000) HENTRY(IEN)
              NBSL=1
              GO TO 100
            ENDIF
          ENDIF
        ENDIF
        IF(NBST .EQ. 0) THEN
*----
*  Find Tracking for volume of mixtures to combine
*----
          IF(IENTRY(IEN) .LE. 2 .AND. JENTRY(IEN) .EQ. 2) THEN
            CALL LCMGTC(KENTRY(IEN),'SIGNATURE',12,1,HSIGN)
            IF(HSIGN.EQ.'L_TRACK') THEN
              IPTRK=KENTRY(IEN)
              WRITE(IOUT,6001) HENTRY(IEN)
              NBST=1
              GO TO 100
            ENDIF
          ENDIF
        ENDIF
 100    CONTINUE
        IF(NBSL+NBST.EQ.2) GO TO 105
      ENDDO
  105 CONTINUE
*----
*  Get information about mixtures on IPLIB
*----
      CALL LCMGET(IPLIB,'STATE-VECTOR',ISTATE)
      NBMIX=ISTATE(1)
      NBISO=ISTATE(2)
      NGRO=ISTATE(3)
      ALLOCATE(ISONRF(3,NBISO),ISOMIX(NBISO),DENISO(NBISO))
      CALL LCMGET(IPLIB,'ISOTOPERNAME',ISONRF)
      CALL LCMGET(IPLIB,'ISOTOPESMIX ',ISOMIX)
      CALL LCMGET(IPLIB,'ISOTOPESDENS',DENISO)
      ALLOCATE(IDCLM(NBMIX))
      ALLOCATE(IACT(3,NBISO),DENRD(NBISO))
*----
*  Read proceessing option
*----
      CALL CLMGET(IPRINT,NBMIX,NBISO,ISONRF,ISOMIX,
     >            NCLM,IDCLM,IACT,DENRD)
*----
*  Get information about volumes for mixtures on IPTRK
*----
      ALLOCATE(VOLMIX(NCLM))
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      NBREG=ISTATE(1)
      ALLOCATE(MATCOD(NBREG),VOLUME(NBREG))
      CALL LCMGET(IPTRK,'MATCOD      ',MATCOD)
      CALL LCMGET(IPTRK,'VOLUME      ',VOLUME)
*----
*  Find volume of each mixtures to combine and total volume
*----
      VOLTOT=0.0
      DO ICLM=1,NCLM
        VOLMIX(ICLM)=0.0
        DO IREG=1,NBREG
          IF(IDCLM(ICLM) .EQ. MATCOD(IREG)) 
     >       VOLMIX(ICLM)=VOLMIX(ICLM)+VOLUME(IREG)
        ENDDO
        VOLTOT=VOLTOT+VOLMIX(ICLM)
      ENDDO
      DEALLOCATE(MATCOD,VOLUME)
*----
*  Find isotopes associated with first mixture to combine
*  with wame isotope from other mixtures
*----
      DO ISO=1,NBISO
        MIXI=IACT(1,ISO)
        IF(MIXI.EQ.1) THEN
          DENISO(ISO)=DENISO(ISO)*VOLMIX(MIXI)
          DO JSO=1,NBISO
            MIXJ=IACT(1,JSO)
            IF(MIXJ.GT.1 .AND. IACT(2,JSO).EQ.ISO) THEN
              DENISO(ISO)=DENISO(ISO)+DENISO(JSO)*VOLMIX(MIXJ)
            ENDIF
          ENDDO
          DENISO(ISO)=DENISO(ISO)/VOLTOT
        ENDIF
      ENDDO 
      DEALLOCATE(VOLMIX)
*----
*  correct mixture according to SETI or ADDI
*----
      DO ISO=1,NBISO
        MIXI=IACT(1,ISO)
        IF(MIXI.EQ.1) THEN
          IF(IACT(3,ISO).EQ. -2) THEN
            DENISO(ISO)=DENRD(ISO)
          ELSE IF(IACT(3,ISO).EQ. -1) THEN
            DENISO(ISO)=DENISO(ISO)+DENRD(ISO)
          ELSE IF(IACT(3,ISO).EQ.  1) THEN
            DENISO(ISO)=DENISO(ISO)*(1.0+DENRD(ISO))
          ELSE IF(IACT(3,ISO).EQ.  2) THEN
            DENISO(ISO)=DENISO(ISO)*DENRD(ISO)
          ENDIF
          DO JSO=1,NBISO
            MIXJ=IACT(1,JSO)
            IF(MIXJ.GT.1 .AND. IACT(2,JSO).EQ.ISO) THEN
              DENISO(JSO)=DENISO(ISO)
            ENDIF
          ENDDO
        ENDIF
      ENDDO 
      DEALLOCATE(DENRD,IACT)
*----
*  Replace new densities in adequate location in DESISO vector
*----
      ALLOCATE(MASK(NBMIX),MASKL(NGRO))
      CALL XDLSET(MASKL,NBMIX,.FALSE.)
      CALL XDLSET(MASKL,NGRO,.TRUE.)
      DO ICLM=1,NCLM
        MASK(IDCLM(ICLM))=.TRUE.
      ENDDO
      DEALLOCATE(IDCLM)
      CALL LCMPUT(IPLIB,'ISOTOPESDENS',NBISO,2,DENISO)
*----
*  Reset macrolib
*----
      ITSTMP=0
      TMPDAY(1)=0.0
      TMPDAY(2)=0.0
      TMPDAY(3)=0.0
      CALL LCMGET(IPLIB,'ISOTOPESUSED',ISONRF)
      CALL LIBMIX(IPLIB,NBMIX,NGRO,NBISO,ISONRF,ISOMIX,DENISO,MASK,
     >            MASKL,ITSTMP,TMPDAY)
      DEALLOCATE(MASKL,MASK)
      DEALLOCATE(ISONRF,ISOMIX,DENISO)
      RETURN
*----
*  FORMATS
*----
 6000 FORMAT('LIBRARY is identified as  : ',A12)
 6001 FORMAT('TRACKING is identified as : ',A12)
      END
