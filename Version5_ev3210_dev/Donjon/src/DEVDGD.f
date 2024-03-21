*DECK DEVDGD
      SUBROUTINE DEVDGD(IPDEV,NROD,DGRP,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Create rod-device group directories on the device data structure.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* D. Sekki
*
*Parameters: input
* IPDEV  pointer to device information.
* NROD   total number of rod-devices.
* DGRP   total number of rod-device groups.
* IMPX   printing index (=0 for no print).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPDEV
      INTEGER NROD,DGRP,IMPX
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOUT=6)
      CHARACTER TEXT*12
      INTEGER RODID(NROD)
      DOUBLE PRECISION DFLOT
      TYPE(C_PTR) JPDEV,KPDEV
*----
*  CREATE GROUPS
*----
      JPDEV=LCMLID(IPDEV,'ROD_GROUP',DGRP)
      IGRP=0
      IF(IMPX.GT.0)WRITE(IOUT,1001)
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.3)CALL XABORT('@DEVDGD: KEYWORD GROUP-ID EXPECTED.')
      IF(TEXT.NE.'GROUP-ID')CALL XABORT('@DEVDGD: KEYWORD GROUP-'
     1 //'ID EXPECTED.')
   10 IGRP=IGRP+1
      CALL REDGET(ITYP,JGRP,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.1)CALL XABORT('@DEVDGD: INTEGER GROUP-ID NUMBER'
     1 //' EXPECTED.')
      IF(JGRP.NE.IGRP)THEN
        WRITE(IOUT,*)'@DEVDGD: READ GROUP-ID NUMBER #',JGRP
        WRITE(IOUT,*)'@DEVDGD: EXPECTED GROUP-ID NUMBER #',IGRP
        CALL XABORT('@DEVDGD: WRONG GROUP-ID NUMBER.')
      ENDIF
      IF(JGRP.GT.DGRP)THEN
        WRITE(IOUT,*)'@DEVDGD: GIVEN TOTAL NUMBER OF GROUPS ',DGRP
        WRITE(IOUT,*)'@DEVDGD: READ GROUP-ID NUMBER #',JGRP
        CALL XABORT('@DEVDGD: WRONG GROUP-ID NUMBER.')
      ENDIF
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.3)CALL XABORT('@DEVDGD: KEYWORD EXPECTED.')
*----
*  OPTION ALL
*----
      IF(TEXT.EQ.'ALL')THEN
        KPDEV=LCMDIL(JPDEV,IGRP)
        DO 30 ID=1,NROD
        RODID(ID)=ID
   30   CONTINUE
        CALL LCMPUT(KPDEV,'GROUP-ID',1,1,IGRP)
        CALL LCMPUT(KPDEV,'NUM-ROD',1,1,NROD)
        CALL LCMPUT(KPDEV,'ROD-ID',NROD,1,RODID)
*
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.3)CALL XABORT('@DEVDGD: WRONG INPUT DATA.')
        IF(TEXT.EQ.';')THEN
          IF(IGRP.EQ.DGRP)THEN
            NDG=NROD
            GOTO 100
          ENDIF
          WRITE(IOUT,*)'@DEVDGD: GIVEN TOTAL NUMBER OF GROUPS ',DGRP
          WRITE(IOUT,*)'@DEVDGD: CREATED ONLY NUMBER OF GROUPS ',IGRP
          CALL XABORT('@DEVDGD: WRONG NUMBER OF GROUPS.')
        ELSEIF(TEXT.EQ.'GROUP-ID')THEN
          IF(IMPX.GT.0)WRITE(IOUT,1000)IGRP,NROD
          GOTO 10
        ELSE
          CALL XABORT('@DEVDGD: WRONG KEYWORD '//TEXT)
        ENDIF
*----
*  OPTION ROD-ID
*----
      ELSEIF(TEXT.EQ.'ROD-ID')THEN
        NDG=0
        CALL XDISET(RODID,NROD,0)
        KPDEV=LCMDIL(JPDEV,IGRP)
*
   50   CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.EQ.3)THEN
          IF(TEXT.EQ.';')THEN
            IF(IGRP.EQ.DGRP)GOTO 100
            WRITE(IOUT,*)'@DEVDGD: GIVEN TOTAL NUMBER OF GROUPS ',DGRP
            WRITE(IOUT,*)'@DEVDGD: CREATED ONLY NUMBER OF GROUPS ',IGRP
            CALL XABORT('@DEVDGD: WRONG NUMBER OF GROUPS.')
          ELSEIF(TEXT.EQ.'GROUP-ID')THEN
            IF(IMPX.GT.0)WRITE(IOUT,1000)IGRP,NDG
            GOTO 10
          ELSE
            CALL XABORT('@DEVDGD: WRONG KEYWORD '//TEXT)
          ENDIF
*----
*  ROD-ID NUMBERS
*----
        ELSEIF(ITYP.EQ.1)THEN
          ID=NITMA
          IF((ID.GT.NROD).OR.(ID.LE.0))THEN
            WRITE(IOUT,*)'@DEVDGD: FOR THE GROUP #',IGRP
            WRITE(IOUT,*)'@DEVDGD: READ ROD-ID #',ID
            CALL XABORT('@DEVDGD: WRONG ROD-ID NUMBER.')
          ENDIF
          DO I=1,NROD
            IF(ID.EQ.RODID(I))THEN
              WRITE(IOUT,*)'@DEVDGD: FOR THE GROUP #',IGRP
              WRITE(IOUT,*)'@DEVDGD: REPEATED ROD-ID #',ID
              CALL XABORT('@DEVDGD: WRONG ROD-ID NUMBER.')
            ENDIF
          ENDDO
*
          NDG=NDG+1
          IF(NDG.GT.NROD)THEN
            WRITE(IOUT,*)'@DEVDGD: FOR THE GROUP #',IGRP
            WRITE(IOUT,*)'@DEVDGD: WRONG TOTAL NUMBER OF RODS ',NDG
            CALL XABORT('@DEVDGD: INVALID INPUT OF ROD-DEVICES.')
          ENDIF
          RODID(NDG)=ID
          CALL LCMPUT(KPDEV,'GROUP-ID',1,1,IGRP)
          CALL LCMPUT(KPDEV,'NUM-ROD',1,1,NDG)
          CALL LCMPUT(KPDEV,'ROD-ID',NDG,1,RODID)
        ELSE
          CALL XABORT('@DEVDGD: WRONG INPUT DATA.')
        ENDIF
        GOTO 50
      ELSE
        CALL XABORT('@DEVDGD: WRONG KEYWORD '//TEXT)
      ENDIF
  100 IF(IMPX.GT.0)WRITE(IOUT,1000)IGRP,NDG
      IF(IMPX.GT.0)WRITE(IOUT,1002)DGRP
      RETURN
*
 1000 FORMAT(/1X,' => CREATED A GROUP #',I2.2,
     1        4X,'INCLUDES TOTAL NUMBER OF RODS:',I3)
 1001 FORMAT(/1X,'**  CREATING GROUPS FOR ROD-DEVICES  **')
 1002 FORMAT(/1X,39('-')/1X,'TOTAL NUMBER OF GROUPS CREATED: ',I2)
      END
