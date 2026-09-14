*DECK OUT
      SUBROUTINE OUT(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Simple edition module for TRIVAC-3.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1): create or modification type(L_MACROLIB);
*         HENTRY(2): read-only type(L_FLUX) or type(L_FVIEW);
*         HENTRY(3): read-only type(L_TRACK);
*         HENTRY(4): read-only type(L_MACROLIB);
*         HENTRY(5): optional read-only type(L_GEOM) containing the
*         original geometry;
*         HENTRY(6): optional read-only type(L_GEOM) containing the
*         macrogeometry;
*         HENTRY(7): optional read-only type(L_MATEX);
* IENTRY  type of each LCM object or file:
*         =1 LCM memory object; =2 XSM file; =3 sequential binary file;
*         =4 sequential ascii file.
* JENTRY  access of each LCM object or file:
*         =0 the LCM object or file is created;
*         =1 the LCM object or file is open for modifications;
*         =2 the LCM object or file is open in read-only mode.
* KENTRY  LCM object address or file unit number.
*
*Comments:
* The OUT: calling specifications are:
* MACRO2 := OUT: FLUX TRACK MACRO GEOM :: (out\_data) ;
* where
*   MACRO2 : name of the \emph{lcm} object (type L\_MACROLIB) containing the 
*     extended \emph{macrolib}.
*   FLUX   : name of the \emph{lcm} object (type L\_FLUX) containing a solution
*   TRACK  : name of the \emph{lcm} object (type L\_TRACK) containing a 
*     \emph{tracking}.
*   MACRO  : name of the \emph{lcm} object (type L\_MACROLIB) containing the 
*     reference \emph{macrolib}.
*   GEOM   : name of the \emph{lcm} object (type L\_GEOM) containing the 
*     reference \emph{geometry}.
*   out\_data}] : structure containing the data to module OUT:
* 
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER      NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR)  KENTRY(NENTRY)
      CHARACTER    HENTRY(NENTRY)*12
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      CHARACTER TEXT12*12,TITLE*72,HTRACK*12,HSIGN*12,MACGEO*12
      INTEGER IGP(NSTATE)
      TYPE(C_PTR) IPMAC1,IPMAC2,IPFLUX,IPVAL,IPTRK,IPGEO1,IPGEO2,IPMTX
      INTEGER, DIMENSION(:),ALLOCATABLE :: MAT,IDL
      REAL, DIMENSION(:),ALLOCATABLE :: VOL
*----
*  RECOVER EXTENDED OUTPUT MACROLIB
*----
      IF(NENTRY.LE.1) CALL XABORT('OUT: TWO PARAMETERS EXPECTED.')
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2)) CALL XABORT('OUT: LCM '
     1 //'OBJECT EXPECTED AT LHS.')
      IF((JENTRY(1).NE.0).AND.(JENTRY(1).NE.1)) CALL XABORT('OUT: ENTR'
     1 //'Y IN CREATE OR MODIFICATION MODE EXPECTED.')
      IPMAC2=KENTRY(1)
      IF(JENTRY(1).EQ.0) THEN
        HSIGN='L_MACROLIB'
        CALL LCMPTC(IPMAC2,'SIGNATURE',12,HSIGN)
      ELSE
        CALL LCMGTC(IPMAC2,'SIGNATURE',12,HSIGN)
        IF(HSIGN.NE.'L_MACROLIB') THEN
          TEXT12=HENTRY(1)
          CALL XABORT('OUT: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1    '. L_MACROLIB EXPECTED.')
        ENDIF
      ENDIF
*----
*  RECOVER READ-ONLY LCM OBJECTS
*----
      IPTRK=C_NULL_PTR
      IPMAC1=C_NULL_PTR
      IPGEO1=C_NULL_PTR
      IPFLUX=C_NULL_PTR
      IPVAL=C_NULL_PTR
      IPGEO2=C_NULL_PTR
      IPMTX=C_NULL_PTR
      MACGEO=' '
      DO 10 I=2,NENTRY
      IF((JENTRY(I).NE.2).OR.((IENTRY(I).NE.1).AND.(IENTRY(I).NE.2)))
     1 CALL XABORT('OUT: LCM OBJECT IN READ-ONLY MODE EXPECTED AT RHS.')
      CALL LCMGTC(KENTRY(I),'SIGNATURE',12,HSIGN)
      IF(HSIGN.EQ.'L_TRACK') THEN
         IPTRK=KENTRY(I)
      ELSE IF(HSIGN.EQ.'L_MACROLIB') THEN
         IPMAC1=KENTRY(I)
      ELSE IF((HSIGN.EQ.'L_GEOM').AND.(.NOT.C_ASSOCIATED(IPGEO1))) THEN
         IPGEO1=KENTRY(I)
      ELSE IF(HSIGN.EQ.'L_FLUX') THEN
         IPFLUX=KENTRY(I)
      ELSE IF(HSIGN.EQ.'L_FVIEW') THEN
         IPVAL=KENTRY(I)
      ELSE IF((HSIGN.EQ.'L_GEOM').AND.(C_ASSOCIATED(IPVAL))) THEN
         IPGEO2=KENTRY(I)
         MACGEO=HENTRY(I)
      ELSE IF(HSIGN.EQ.'L_MATEX') THEN
         IPMTX=KENTRY(I)
      ENDIF
   10 CONTINUE
*----
*  CHECK OBJECT CONSISTENCY
*----
      IF(.NOT.C_ASSOCIATED(IPTRK)) THEN
        CALL XABORT('OUT: MISSING READ-ONLY TRACKING OBJECT.')
      ELSE IF(.NOT.C_ASSOCIATED(IPMAC1)) THEN
        CALL XABORT('OUT: MISSING READ-ONLY REFERENCE MACROLIB OBJECT.')
      ENDIF
      IF(C_ASSOCIATED(IPFLUX).AND.C_ASSOCIATED(IPVAL)) THEN
        CALL XABORT('OUT: SELECT A SINGLE L_FLUX OR L_FVIEW OBJECT.')
      ELSE IF(C_ASSOCIATED(IPVAL).AND.(.NOT.C_ASSOCIATED(IPGEO2))) THEN
        CALL XABORT('OUT: REQUIRED MACROGEOMETRY WITH L_FVIEW OBJECT.')
      ENDIF
*----
*  RECOVER REFERENCE MACROLIB INFORMATION.
*----
      CALL LCMGET(IPMAC1,'STATE-VECTOR',IGP)
      NGRP=IGP(1)
      NBMIX=IGP(2)
      NL=IGP(3)
      NBFIS=IGP(4)
      NALBP=IGP(8)
*----
*  FIND TYPE OF TRACKING.
*----
      CALL LCMGTC(IPTRK,'TRACK-TYPE',12,HTRACK)
*----
*  RECOVER GENERAL TRACKING INFORMATION.
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',IGP)
      NEL=IGP(1)
      NUN=IGP(2)
      IF(HTRACK.EQ.'BIVAC') THEN
        IELEM=IGP(8)
        ICOL=IGP(9)
        IBFP=0
      ELSE IF(HTRACK.EQ.'TRIVAC') THEN
        IELEM=IGP(9)
        ICOL=IGP(10)
        IBFP=0
      ELSE IF(HTRACK.EQ.'SN') THEN
        IELEM=IGP(8)
        ICOL=0
        IBFP=IGP(31)
      ELSE
        ICOL=0
        IBFP=0
      ENDIF
      MAXNEL=NEL
      CALL LCMLEN(IPTRK,'KEYFLX',LKFL,ITYLCM)
      ALLOCATE(MAT(NEL),VOL(NEL),IDL(LKFL))
      CALL LCMGET(IPTRK,'MATCOD',MAT)
      CALL LCMGET(IPTRK,'VOLUME',VOL)
      CALL LCMGET(IPTRK,'KEYFLX',IDL)
      CALL LCMLEN(IPTRK,'TITLE',LENGT,ITYLCM)
      IF(LENGT.GT.0) THEN
         CALL LCMGTC(IPTRK,'TITLE',72,TITLE)
         CALL LCMPTC(IPMAC2,'TITLE',72,TITLE)
      ELSE
         TITLE='*** NO TITLE PROVIDED ***'
      ENDIF
*----
*  EDITION.
*----
      IF(C_ASSOCIATED(IPFLUX)) THEN
        CALL OUTDRV(IPGEO1,IPMAC1,IPFLUX,IPMAC2,IPMTX,MAXNEL,NBMIX,NL,
     1  NBFIS,NGRP,NEL,NUN,NALBP,HTRACK,IELEM,ICOL,MAT,VOL,IDL,TITLE,
     2  IBFP)
      ELSE IF(C_ASSOCIATED(IPVAL)) THEN
        CALL OUTVOX(IPGEO2,IPMAC1,IPVAL,IPMAC2,NBMIX,NL,NBFIS,NGRP,
     1  NALBP,TITLE,MACGEO)
      ENDIF
*----
*  RELEASE GENERAL TRACKING INFORMATION.
*----
      DEALLOCATE(IDL,VOL,MAT)
      RETURN
      END
