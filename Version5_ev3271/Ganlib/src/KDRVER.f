*DECK KDRVER
      SUBROUTINE KDRVER(REV,DATE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To extract CVS or SVN version and production date.
*
*Copyright:
* Copyright (C) 2006 Ecole Polytechnique de Montreal
*
*Author(s): A. Hebert
*
*Parameters: output
* REV     revision character identification
* DATE    revision character date 
*
*-----------------------------------------------------------------------
*
      CHARACTER REV*48,DATE*64
*
      REV='Version 5.0.10 ($Revision: 3196 $)'
      DATE='$Date: 2024-03-06 07:25:15 -0500 (Wed, 06 Mar 2024) $'
      IF(REV(22:).EQ.'ion$)') THEN
*        CVS or SVN keyword expansion not performed
         REV='Version 5.0.10'
         DATE='March 6, 2024'
      ENDIF
      RETURN
      END
