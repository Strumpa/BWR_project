*DECK NSSF
      SUBROUTINE NSSF(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Flux solution for a nodal method (NEM or ANM).
*
*Copyright:
* Copyright (C) 2021 Ecole Polytechnique de Montreal
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
*         HENTRY(1): create type(L_FLUX) nodal flux;
*         HENTRY(2): read-only type(L_TRACK) nodal tracking;
*         HENTRY(3): read-only type(L_MACROLIB) nodal macrolib.
* IENTRY  type of each LCM object or file:
*         =1 LCM memory object; =2 XSM file; =3 sequential binary file;
*         =4 sequential ascii file.
* JENTRY  access of each LCM object or file:
*         =0 the LCM object or file is created;
*         =1 the LCM object or file is open for modifications;
*         =2 the LCM object or file is open in read-only mode.
* KENTRY  LCM object address or file unit number.
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
      TYPE(C_PTR) IPFLX,IPTRK,IPMAC
      CHARACTER HSIGN*12,TEXT4*4,TEXT12*12,HSMG*131,BNDTL*12
      LOGICAL LNODF
      INTEGER ISTATE(NSTATE)
      REAL REALIR
      DOUBLE PRECISION DBLLIR
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ITRIAL
*----
*  PARAMETER VALIDATION.
*----
      IF(NENTRY.NE.3) CALL XABORT('NSSF: 3 PARAMETERS EXPECTED.')
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2)) CALL XABORT('NSSF: LCM'
     1 //' OBJECT EXPECTED AT LHS.')
      DO IEN=2,3
        IF((IENTRY(IEN).NE.1).AND.(IENTRY(IEN).NE.2)) CALL XABORT('NSS'
     1  //'F: LCM OBJECT EXPECTED AT LHS.')
        IF(JENTRY(IEN).NE.2) CALL XABORT('NSSF: ENTRY IN READ-ONLY MOD'
     1  //'E EXPECTED.')
        CALL LCMGTC(KENTRY(IEN),'SIGNATURE',12,1,HSIGN)
        TEXT12=HENTRY(IEN)
        IF(IEN.EQ.2) THEN
          IF(HSIGN.NE.'L_TRACK') THEN
            CALL XABORT('NSSF: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1      '. L_TRACK EXPECTED.')
          ENDIF
          IPTRK=KENTRY(2)
          CALL LCMPTC(KENTRY(1),'LINK.TRACK',12,1,HENTRY(2))
        ELSE IF(IEN.EQ.3) THEN
          IF(HSIGN.NE.'L_MACROLIB') THEN
            CALL XABORT('NSSF: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1      '. L_MACROLIB EXPECTED.')
          ENDIF
          IPMAC=KENTRY(3)
          CALL LCMPTC(KENTRY(1),'LINK.MACRO',12,1,HENTRY(3))
        ENDIF
      ENDDO
      CALL LCMGTC(IPTRK,'TRACK-TYPE',12,1,TEXT12)
      IF(TEXT12.NE.'TRIVAC') CALL XABORT('NSSF: TRIVAC TRACKING EXPECT'
     1 //'ED.')
*----
*  PROCESS MACROLIB AND TRACKING.
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      NEL=ISTATE(1)
      NUN=ISTATE(2)
      NMIX=ISTATE(4)
      NADI=ISTATE(33)
      IGMAX=ISTATE(39)
      ICHX=ISTATE(12)
      IF((ICHX.LT.4).OR.(ICHX.GT.6)) THEN
        CALL XABORT('NSSF: CMFD, NEM OR ANM DISCRETIZATION EXPECTED.')
      ENDIF
      IF(ISTATE(6).EQ.2) THEN
        IDIM=1
      ELSE IF((ISTATE(6).EQ.5).AND.(ICHX.EQ.6)) THEN
        IDIM=2
      ELSE IF((ISTATE(6).EQ.7).AND.(ICHX.EQ.6)) THEN
        IDIM=3
      ELSE
        CALL XABORT('NSSF: 1D SLAB/2D-3D CARTESIAN GEOMETRY EXPECTED.')
      ENDIF
      IF(ISTATE(38).NE.0) CALL XABORT('NSSF: LUMP OPTION FORBIDDEN.')
      CALL LCMGET(IPMAC,'STATE-VECTOR',ISTATE)
      NG=ISTATE(1)
      IF(ISTATE(2).NE.NMIX) THEN
        WRITE(HSMG,'(39HNSSF: INVALID NUMBER OF MIXTURES (GEOM=,I5,
     1  10H MACROLIB=,I5,2H).)') NMIX,ISTATE(2)
        CALL XABORT(HSMG)
      ENDIF
*----
*  CREATE OR RECOVER THE FLUX.
*----
      IPFLX=KENTRY(1)
      IF(JENTRY(1).EQ.0) THEN
        HSIGN='L_FLUX'
        CALL LCMPTC(IPFLX,'SIGNATURE',12,1,HSIGN)
      ELSE IF(JENTRY(1).EQ.1) THEN
        CALL LCMGTC(IPFLX,'SIGNATURE',12,1,HSIGN)
        TEXT12=HENTRY(IEN)
        IF(HSIGN.NE.'L_FLUX') THEN
          CALL XABORT('NSSF: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1    '. L_FLUX EXPECTED.')
        ENDIF
        CALL LCMGET(IPFLX,'STATE-VECTOR',ISTATE)
        IF(ISTATE(1).NE.NG) THEN
          WRITE(HSMG,'(41HNSSF: INVALID NUMBER OF GROUPS (MACROLIB=,I5,
     1    6H FLUX=,I5,2H).)') NG,ISTATE(1)
          CALL XABORT(HSMG)
        ELSE IF(ISTATE(2).NE.NUN) THEN
          WRITE(HSMG,'(43HNSSF: INVALID NUMBER OF UNKNOWNS (TRACKING=,
     1    I10,6H FLUX=,I10,2H).)') NUN,ISTATE(2)
          CALL XABORT(HSMG)
        ENDIF
      ELSE
        CALL XABORT('NSSF: FLUX IN CREATE OR MODIFICATION MODE EXPECTE'
     1  //'D.')
      ENDIF
*---
*  READ DATA
*---
      ALLOCATE(ITRIAL(NMIX,NG))
      IPRINT=1
      ICL1=3
      ICL2=3
      MAXNOD=300
      MAXTHR=0
      MAXOUT=100
      EPSNOD=1.0E-6
      EPSTHR=1.0E-6
      EPSOUT=1.0E-5
      LNODF=.FALSE.
      BB2=0.0
      BNDTL='quadratic'
      NPASS=3
      ITRIAL(:,:)=1
   10 CALL REDGET(ITYPLU,INTLIR,REALIR,TEXT4,DBLLIR)
      IF(ITYPLU.EQ.10) GO TO 100
   20 IF(ITYPLU.NE.3) CALL XABORT('NSSF: READ ERROR - CHARACTER VARIAB'
     > //'LE EXPECTED')
      IF(TEXT4.EQ.';') THEN
        GO TO 100
      ELSE IF(TEXT4.EQ.'EDIT') THEN
        CALL REDGET(ITYPLU,IPRINT,REALIR,TEXT4,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('NSSF: INTEGER DATA EXPECTED(1).')
      ELSE IF((TEXT4.EQ.'VAR1').OR.(TEXT4.EQ.'ACCE')) THEN
        CALL REDGET(ITYPLU,ICL1,REALIR,TEXT4,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('NSSF: INTEGER DATA EXPECTED(2).')
        CALL REDGET(ITYPLU,ICL2,REALIR,TEXT4,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('NSSF: INTEGER DATA EXPECTED(3).')
      ELSE IF(TEXT4=='ADI') THEN
        CALL REDGET(ITYPLU,NADI,REALIR,TEXT4,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('NSSF: INTEGER DATA EXPECTED(5).')
      ELSE IF(TEXT4=='NUPD') THEN
        ! maximum number of nodal correction iterations
   30   CALL REDGET(ITYPLU,INTLIR,REALIR,TEXT4,DBLLIR)
        IF(ITYPLU.EQ.1) THEN
          MAXNOD=INTLIR
          CALL REDGET(ITYPLU,INTLIR,REALIR,TEXT4,DBLLIR)
          IF(ITYPLU.EQ.1) THEN
            NPASS=INTLIR
            GO TO 30
          ENDIF
        ELSE IF(ITYPLU.EQ.2) THEN
          EPSNOD=REALIR
        ELSE
          GO TO 20
        ENDIF
        GO TO 30
      ELSE IF(TEXT4=='EXTE') THEN
        ! maximum number and convergence criterion of Keff iterations
   40   CALL REDGET(ITYPLU,INTLIR,REALIR,TEXT4,DBLLIR)
        IF(ITYPLU.EQ.1) THEN
          MAXOUT=INTLIR
        ELSE IF(ITYPLU.EQ.2) THEN
          EPSOUT=REALIR
        ELSE
          GO TO 20
        ENDIF
        GO TO 40
      ELSE IF(TEXT4=='THER') THEN
        ! maximum number and convergence criterion of thermal iterations
   50   CALL REDGET(ITYPLU,INTLIR,REALIR,TEXT4,DBLLIR)
        IF(ITYPLU.EQ.1) THEN
          MAXTHR=INTLIR
        ELSE IF(ITYPLU.EQ.2) THEN
          EPSTHR=REALIR
        ELSE
          GO TO 20
        ENDIF
        GO TO 50
      ELSE IF(TEXT4.EQ.'NODF') THEN
        LNODF=.TRUE.
      ELSE IF(TEXT4=='LEAK') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,BNDTL,DBLLIR)
        IF(ITYPLU/=3) CALL XABORT('NSSF: READ ERROR - CHARACTER VARIAB'
     >  //'LE EXPECTED')
        IF((BNDTL.NE.'flat').AND.(BNDTL.NE.'quadratic')) THEN
         CALL XABORT('NSSF: flat OR quadratic KEYWORD EXPECTED')
        ENDIF
      ELSE IF(TEXT4.EQ.'BUCK') THEN
        CALL REDGET(ITYPLU,INTLIR,BB2,TEXT4,DBLLIR)
        IF(ITYPLU.NE.2) CALL XABORT('NSSF: READ ERROR - REAL VARIABLE '
     >  //'EXPECTED')
      ELSE
        CALL XABORT('NSSF: ILLEGAL KEYWORD '//TEXT4)
      ENDIF
      GO TO 10
  100 IF(IGMAX.GT.NG) CALL XABORT('NSSF: IGMAX>NG.')
      IF(IPRINT.GT.0) THEN
        WRITE(6,'(/47H NSSF: number of transverse current iterations=,
     >  I3)') NPASS
      ENDIF
      IF(IGMAX.GT.0) ITRIAL(:NMIX,IGMAX:NG)=2
      CALL NSSDRV(IPTRK,IPMAC,IPFLX,ICHX,IDIM,NUN,NG,NEL,NMIX,ITRIAL,
     > ICL1,ICL2,NADI,EPSNOD,MAXNOD,EPSTHR,MAXTHR,EPSOUT,MAXOUT,LNODF,
     > BNDTL,NPASS,BB2,IPRINT)
      DEALLOCATE(ITRIAL)
      RETURN
      END
