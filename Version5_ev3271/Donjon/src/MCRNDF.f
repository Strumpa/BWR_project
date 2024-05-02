*DECK MCRNDF
      SUBROUTINE MCRNDF(IMPX,NBISO,ISO,IBM,HNOMIS,IPLIB,MY1,MY2,YLDS,
     1 IADRY,ISTYP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Store records PYNAM, PYMIX and PYIELD into a Microlib.
*
*Copyright:
* Copyright (C) 2022 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* IMPX    print parameter (equal to zero for no print).
* NBISO  number of particularized isotopes.
* ISO     particularized isotope index.
* IBM     material mixture.
* HNOMIS  array containing the names of the particularized isotopes.
* IPLIB   address of the output microlib LCM object.
* MY1     number of fissile isotopes including macroscopic sets.
* MY2     number of fission fragment.
* YLDS    fission yields.
* IADRY   index in YLDS (<0: fission product; >0: fissile isotope).
*
*Parameters: output
* ISTYP   type of isotope ISO (=1: stable;=2: fissile; =3: fission
*         product).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB
      INTEGER IMPX,NBISO,ISO,IBM,MY1,MY2,ISTYP,IADRY(NBISO)
      DOUBLE PRECISION YLDS(MY1,MY2)
      CHARACTER(LEN=24) HNOMIS(NBISO)
*----
*  LOCAL VARIABLES
*----
      INTEGER I,IY1,IY2,JSO
*----
*  ALLOCATABLE AYYAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IPYMIX
      REAL, ALLOCATABLE, DIMENSION(:) :: PYIELD
      CHARACTER(LEN=24), ALLOCATABLE, DIMENSION(:) :: HPYNAM
*
      IF(IADRY(ISO).GT.0) THEN
*       ISO is a fissile isotope
        ISTYP=2
      ELSE IF(IADRY(ISO).LT.0) THEN
*       ISO is a fission product
        ISTYP=3
        IY2=-IADRY(ISO)
        IF(IY2.GT.MY2) CALL XABORT('MCRNDF: MY2 OVERFLOW.')
        ALLOCATE(HPYNAM(MY1),IPYMIX(MY1),PYIELD(MY1))
        HPYNAM(:MY1)=' '
        IPYMIX(:MY1)=0
        PYIELD(:MY1)=0.0
        IF(IMPX.GT.2) THEN
          WRITE(6,'(25H MCRNDF: fission product=,A24,9H mixture=,I8)')
     1    HNOMIS(ISO),IBM
        ENDIF
        DO JSO=1,NBISO
          IF(IADRY(JSO).GT.0) THEN
            IY1=IADRY(JSO)
            IF(IY1.GT.MY1) CALL XABORT('MCRNDF: MY1 OVERFLOW.')
            HPYNAM(IY1)=HNOMIS(JSO)
            IPYMIX(IY1)=IBM
            PYIELD(IY1)=REAL(YLDS(IY1,IY2))
            IF(IMPX.GT.2) THEN
              WRITE(6,'(9X,16Hfissile isotope(,I4,2H)=,A24,9H mixture=,
     1        I8)') IY1,HPYNAM(IY1),IPYMIX(IY1)
            ENDIF
          ENDIF
        ENDDO
        CALL LCMPTC(IPLIB,'PYNAM',8,MY1,HPYNAM(:8))
        CALL LCMPUT(IPLIB,'PYMIX',MY1,1,IPYMIX)
        CALL LCMPUT(IPLIB,'PYIELD',MY1,2,PYIELD)
        IF(IMPX.GT.2) THEN
          WRITE(6,'(3X,7HPYIELD=,1P,8E12.4/(8X,10E12.4))') (PYIELD(I),
     1    I=1,MY1)
        ENDIF
        DEALLOCATE(PYIELD,IPYMIX,HPYNAM)
      ENDIF
      RETURN
      END
