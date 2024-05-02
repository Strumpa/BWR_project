*DECK MPO
      SUBROUTINE MPO(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Creation and construction of a MPO database object.
*
*Copyright:
* Copyright (C) 2022 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): A. Hebert
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1) MPO database object;
*         HENTRY(I) I>1 read-only type(L_BURNUP, L_LIBRARY, L_EDIT
*         or MPO file).
* IENTRY  type of each LCM object or file:
*         =1 LCM memory object; =2 XSM file; =3 sequential binary file;
*         =4 sequential ascii file; =6 HDF5 file.
* JENTRY  access of each LCM object or file:
*         =0 the LCM object or file is created;
*         =1 the LCM object or file is open for modifications;
*         =2 the LCM object or file is open in read-only mode.
* KENTRY  LCM object address or file unit number.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE hdf5_wrap
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER      NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR)  KENTRY(NENTRY)
      CHARACTER    HENTRY(NENTRY)*12
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40,MAXPAR=50,MAXISO=800,NKEYS=8,NREAK=10,
     1 MAXMAC=2,MAXREA=50,MAXCAD=153)
      TYPE(C_PTR) IPMPO,IPLB1,IPLB2,IPDEPL,IPEDIT,IPRHS
      CHARACTER TEXT24*24,TEXT8*8,TEXT12*12,TEXT20*20,HSMPO*132,
     1 HSMG*131,HEDIT*12,CDIRO*12,RECNAM*72,HSIGN*12,KEYWRD(NKEYS)*4,
     2 NOMISO(MAXISO)*20,NOMEVO(MAXISO)*12,REANAM(NREAK)*20,
     3 NOMREA(MAXREA)*20,REV*48,DATE*64
      DOUBLE PRECISION DFLOTT
      LOGICAL LINIT,LWARN,LGNEW(MAXPAR)
      INTEGER ISTATE(NSTATE),TYPISO(MAXISO),MUPLET(MAXPAR),DIMSR(5),
     1 RANK,TYPE
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NVALUE,PARADR,PARADL,
     1 REACTION,ISOTOPE,ADDRISO,DIMS_MPO
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: OUPUTID,OUPUTID2
      REAL, ALLOCATABLE, DIMENSION(:) :: TIMES,VOL,ENERG
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: PARFMT
      CHARACTER(LEN=24), ALLOCATABLE, DIMENSION(:) :: PARTYP,PARKEY,
     1 PARCAD,PARTYL,PARKEL,PARCAL
*----
*  DATA STATEMENTS
*----
      DATA KEYWRD/'NOML','PARA','LOCA','ISOT','MACR','REAC','NAME',
     1            ';   '/
      DATA REANAM/'Total               ','Absorption          ',
     1            'Diffusion           ','Fission             ',
     2            'FissionSpectrum     ','Nexcess             ',
     3            'NuFission           ','Scattering          ',
     4            'CaptureEnergyCapture','FissionEnergyFission'/
*----
*  PARAMETER VALIDATION.
*----
      LINIT=.FALSE.
      IF(NENTRY.EQ.0) CALL XABORT('MPO: PARAMETERS EXPECTED.')
      IF((IENTRY(1).EQ.6).AND.(JENTRY(1).EQ.0)) THEN
         IPMPO=KENTRY(1)
         LINIT=.TRUE.
         CALL KDRVER(REV,DATE)
         HSMPO='DRAGON5 generated file'
         CALL hdf5_create_group(IPMPO,"info")
         CALL hdf5_write_data(IPMPO,"/info/APOLLO3_VERSION",
     1   TRIM(HSMPO))
         IVERS=1
         CALL hdf5_write_data(IPMPO,"/info/MPO_VERSION",IVERS)
         HSMPO='MPO LIBRARY'//REV//DATE
         CALL hdf5_write_data(IPMPO,"/info/MPO_CREATION_INFO",
     1   TRIM(HSMPO))
      ELSE IF(IENTRY(1).EQ.6) THEN
         IPMPO=KENTRY(1)
         CALL hdf5_info(IPMPO,"info/MPO_VERSION",RANK,TYPE,NBYTE,DIMSR)
         IF(TYPE.EQ.99) THEN
            TEXT12=HENTRY(1)
            CALL XABORT('MPO: HDF FILE '//TEXT12//' CANNOT BE READ.')
         ENDIF
         LINIT=.FALSE.
      ELSE
         CALL XABORT('MPO: MPO HDF5 OBJECT EXPECTED.')
      ENDIF
      TYPISO(:MAXISO)=0
      IPLB1=C_NULL_PTR
      IPLB2=C_NULL_PTR
      IPDEPL=C_NULL_PTR
      IPEDIT=C_NULL_PTR
      IPRHS=C_NULL_PTR
      DO 10 I=2,NENTRY
      IF(JENTRY(I).NE.2) CALL XABORT('MPO: READ-ONLY RHS EXPECTED.')
      IF(IENTRY(I).LE.2) THEN
         CALL LCMGTC(KENTRY(I),'SIGNATURE',12,1,HSIGN)
         IF(HSIGN.EQ.'L_LIBRARY') THEN
            IF(.NOT.C_ASSOCIATED(IPLB1)) THEN
               IPLB1=KENTRY(I)
            ELSE
               IF(.NOT.C_ASSOCIATED(IPLB2)) IPLB2=KENTRY(I)
            ENDIF
         ELSE IF(HSIGN.EQ.'L_BURNUP') THEN
            IPDEPL=KENTRY(I)
         ELSE IF(HSIGN.EQ.'L_EDIT') THEN
            IPEDIT=KENTRY(I)
         ENDIF
      ELSE IF(IENTRY(I).EQ.6) THEN
         IPRHS=KENTRY(I)
      ELSE
         CALL XABORT('MPO: LCM OR HDF5 OBJECTS EXPECTED AT RHS.')
      ENDIF
   10 CONTINUE
*----
*  READ THE INPUT DATA.
*----
*     DEFAULT OPTIONS:
      ALLOCATE(PARADR(MAXPAR+1),PARCAD(MAXCAD),NVALUE(MAXPAR))
      ALLOCATE(PARKEL(MAXPAR),PARTYL(MAXPAR),PARADL(MAXPAR+1),
     1 PARCAL(MAXCAD))
      IMPX=1
      IF(LINIT) THEN
         ALLOCATE(PARKEY(MAXPAR),PARFMT(MAXPAR),PARTYP(MAXPAR))
         NPAR=0
         NPCHR=0
         NLOC=0
         NPCHL=0
         NISO=0
         NMIL=0
         PARADR(1)=0
         PARADL(1)=0
         DO 15 IKEY=1,NREAK
            NOMREA(IKEY)=REANAM(IKEY)
   15    CONTINUE
         NREA=NREAK
      ELSE
         GO TO 300
      ENDIF
   20 CALL REDGET(INDIC,NITMA,FLOTT,TEXT8,DFLOTT)
      IF(INDIC.NE.3) CALL XABORT('MPO: CHARACTER DATA EXPECTED(1).')
   30 IF(TEXT8.EQ.'EDIT') THEN
*        READ THE PRINT INDEX.
         CALL REDGET(INDIC,IMPX,FLOTT,TEXT24,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('MPO: INTEGER DATA EXPECTED(1).')
      ELSE IF(TEXT8.EQ.'COMM') THEN
         CALL REDGET(INDIC,NITMA,FLOTT,HSMPO,DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('MPO: COMMENTS EXPECTED.')
         CALL hdf5_delete(IPMPO,"info/MPO_CREATION_INFO")
         CALL hdf5_write_data(IPMPO,"info/MPO_CREATION_INFO",
     1   TRIM(HSMPO))
      ELSE IF(TEXT8.EQ.'PARA') THEN
         NPAR=NPAR+1
         IF(NPAR.GT.MAXPAR) CALL XABORT('MPO: TOO MANY PARAMETERS.')
         PARKEY(NPAR)=' '
         CALL REDGET(INDIC,NITMA,FLOTT,PARKEY(NPAR),DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('MPO: CHARACTER DATA EXPECTED'
     1   //'(4).')
         DO 50 I=1,NPAR-1
         IF(PARKEY(NPAR).EQ.PARKEY(I)) CALL XABORT('MPO: PARKEY '//
     1   PARKEY(NPAR)//' ALREADY DEFINED(1).')
   50    CONTINUE
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT24,DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('MPO: CHARACTER DATA EXPECTED'
     1   //'(4).')
         PARTYP(NPAR)=TEXT24
         IF(TEXT24.EQ.'TEMP') THEN
           IF(NPCHR+2.GT.MAXCAD) CALL XABORT('MPO: MAXCAD OVERFLOW(1).')
           NPCHR=NPCHR+1
           CALL REDGET(INDIC,NITMA,FLOTT,PARCAD(NPCHR),DFLOTT)
           IF(INDIC.NE.3) CALL XABORT('MPO: CHARACTER DATA EXPECTED'
     1     //'(5).')
           NPCHR=NPCHR+1
           CALL REDGET(INDIC,NITMA,FLOTT,TEXT24,DFLOTT)
           IF(INDIC.NE.1) CALL XABORT('MPO: INTEGER DATA EXPECTE'//
     1     'D(2).')
           WRITE(PARCAD(NPCHR),'(3HMIX,I9.9)') NITMA
           PARFMT(NPAR)='FLOAT'
           PARTYP(NPAR)='TEMPERATURE'
         ELSE IF(TEXT24.EQ.'CONC') THEN
           IF(NPCHR+3.GT.MAXCAD) CALL XABORT('MPO: MAXCAD OVERFLOW(2).')
           NPCHR=NPCHR+1
           CALL REDGET(INDIC,NITMA,FLOTT,PARCAD(NPCHR),DFLOTT)
           IF(INDIC.NE.3) CALL XABORT('MPO: CHARACTER DATA EXPECTED'
     1     //'(6).')
           NPCHR=NPCHR+1
           CALL REDGET(INDIC,NITMA,FLOTT,PARCAD(NPCHR),DFLOTT)
           IF(INDIC.NE.3) CALL XABORT('MPO: CHARACTER DATA EXPECTED'
     1     //'(7).')
           NPCHR=NPCHR+1
           CALL REDGET(INDIC,NITMA,FLOTT,TEXT24,DFLOTT)
           IF(INDIC.NE.1) CALL XABORT('MPO: INTEGER DATA EXPECTE'//
     1     'D(3).')
           WRITE(PARCAD(NPCHR),'(3HMIX,I9.9)') NITMA
           PARFMT(NPAR)='FLOAT'
           PARTYP(NPAR)='CONCENTRATION_MATERIAL'
         ELSE IF(TEXT24.EQ.'IRRA') THEN
           PARFMT(NPAR)='FLOAT'
           PARTYP(NPAR)='BURNUP'
         ELSE IF(TEXT24.EQ.'FLUX') THEN
           PARFMT(NPAR)='FLOAT'
         ELSE IF(TEXT24.EQ.'FLUB') THEN
           PARFMT(NPAR)='FLOAT'
         ELSE IF(TEXT24.EQ.'PUIS') THEN
           PARFMT(NPAR)='FLOAT'
         ELSE IF(TEXT24.EQ.'TIME') THEN
           PARFMT(NPAR)='FLOAT'
           PARTYP(NPAR)='TIME'
         ELSE IF(TEXT24.EQ.'VALU') THEN
           CALL REDGET(INDIC,NITMA,FLOTT,TEXT8,DFLOTT)
           IF(INDIC.NE.3) CALL XABORT('MPO: CHARACTER DATA EXPECTED'
     1     //'(8).')
           IF(TEXT8.EQ.'REAL')THEN
             PARFMT(NPAR)='FLOAT'
           ELSEIF(TEXT8.EQ.'CHAR')THEN
             PARFMT(NPAR)='STRING'
           ELSEIF(TEXT8.EQ.'INTE')THEN
             PARFMT(NPAR)='INTEGER'
           ELSE
             CALL XABORT('MPO: INVALID KEYWORD='//TEXT24//'(1).')
           ENDIF
         ELSE
           CALL XABORT('MPO: INVALID KEYWORD='//TEXT24//'(2).')
         ENDIF
         NVALUE(NPAR)=0
         PARADR(NPAR+1)=NPCHR
      ELSE IF(TEXT8.EQ.'LOCA') THEN
         NLOC=NLOC+1
         IF(NLOC.GT.MAXPAR) CALL XABORT('MPO: TOO MANY LOCAL VAR'//
     1   'IABLES.')
         CALL REDGET(INDIC,NITMA,FLOTT,PARKEL(NLOC),DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('MPO: CHARACTER DATA EXPECTED'
     1   //'(10).')
         DO 70 I=1,NLOC-1
         IF(PARKEL(NLOC).EQ.PARKEL(I)) CALL XABORT('MPO: PARKEY '//
     1   PARKEL(NLOC)//' ALREADY DEFINED(2).')
   70    CONTINUE
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT24,DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('MPO: CHARACTER DATA EXPECTED'
     1   //'(11).')
         IF(TEXT24.EQ.'CONC') THEN
           IF(NPCHL+1.GT.MAXCAD) CALL XABORT('MPO: MAXCAD OVERFLOW(3).')
           NPCHL=NPCHL+1
           CALL REDGET(INDIC,NITMA,FLOTT,PARCAL(NPCHL),DFLOTT)
           IF(INDIC.NE.3) CALL XABORT('MPO: CHARACTER DATA EXPECTED'
     1     //'(12).')
         ELSE IF((TEXT24.NE.'IRRA').AND.(TEXT24.NE.'FLUG').AND.
     1           (TEXT24.NE.'FLUB').AND.(TEXT24.NE.'PUIS').AND.
     2           (TEXT24.NE.'MASL').AND.(TEXT24.NE.'FLUX').AND.
     3           (TEXT24.NE.'EQUI').AND.(TEXT24.NE.'TEMP')) THEN
           CALL XABORT('MPO: INVALID KEYWORD='//TEXT24//'(3).')
         ENDIF
         PARTYL(NLOC)=TEXT24
         PARADL(NLOC+1)=NPCHL
      ELSE IF(TEXT8.EQ.'ISOT') THEN
   80    CALL REDGET(INDIC,NITMA,FLOTT,TEXT8,DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('MPO: CHARACTER DATA EXPECTED'
     1   //'(13).')
         DO 90 IKEY=1,NKEYS
         IF(TEXT8.EQ.KEYWRD(IKEY)) GO TO 30
   90    CONTINUE
         IF(TEXT8.EQ.'TOUT') THEN
            CALL COMISO(-1,MAXISO,IPLB1,NISO,NOMISO,NOMEVO,TYPISO)
            GO TO 20
         ELSE IF(TEXT8.EQ.'FISS') THEN
            CALL COMISO(-2,MAXISO,IPLB1,NISO,NOMISO,NOMEVO,TYPISO)
         ELSE IF(TEXT8.EQ.'PF') THEN
            CALL COMISO(-3,MAXISO,IPLB1,NISO,NOMISO,NOMEVO,TYPISO)
         ELSE IF(TEXT8.EQ.'MILI') THEN
            CALL REDGET(INDIC,NITMA,FLOTT,TEXT24,DFLOTT)
            IF(INDIC.NE.1) CALL XABORT('MPO: INTEGER DATA EXPECTE'//
     1      'D(4).')
            CALL COMISO(NITMA,MAXISO,IPLB1,NISO,NOMISO,NOMEVO,TYPISO)
         ELSE
            DO 100 IKEY=1,NKEYS
            IF(TEXT8.EQ.KEYWRD(IKEY)) GO TO 30
  100       CONTINUE
            NISO=NISO+1
            IF(NISO.GT.MAXISO) CALL XABORT('MPO: TOO MANY ISOTOPES.')
            NOMISO(NISO)=TEXT8
            TYPISO(NISO)=0
         ENDIF
         GO TO 80
      ELSE IF(TEXT8.EQ.'REAC') THEN
  110    CALL REDGET(INDIC,NITMA,FLOTT,TEXT20,DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('MPO: CHARACTER DATA EXPECTED'
     1   //'(16).')
         IF(TEXT20.EQ.';') GO TO 160
         DO 120 IKEY=1,NKEYS
         IF(TEXT20.EQ.KEYWRD(IKEY)) GO TO 30
  120    CONTINUE
         DO 130 IKEY=1,NREAK
         IF(TEXT20.EQ.REANAM(IKEY)) GO TO 110
  130    CONTINUE
         NREA=NREA+1
         IF(NREA.GT.MAXREA) CALL XABORT('MPO: TOO MANY REACTIONS.')
         NOMREA(NREA)=TEXT20
         GO TO 110
      ELSE IF(TEXT8.EQ.';') THEN
         GO TO 160
      ELSE
         CALL XABORT('MPO: INVALID KEYWORD='//TEXT8//'(4).')
      ENDIF
      GO TO 20
*
* ADD MACROSCOPIC RESIDUAL TO ISOTOPES.
  160 NISO=NISO+1
      IF(NISO.GT.MAXISO) CALL XABORT('MPO: TOO MANY ISOTOPES.')
      NOMISO(NISO)='TotalResidual_mix'
      TYPISO(NISO)=0
      MUPLET(:NPAR)=-99
*
* CREATE MPO GROUPS.
      CALL hdf5_create_group(IPMPO,"/contents")
      CALL hdf5_create_group(IPMPO,"/contents/isotopes")
      CALL hdf5_create_group(IPMPO,"/contents/mixtures")
      CALL hdf5_create_group(IPMPO,"/contents/reactions")
      CALL hdf5_create_group(IPMPO,"/parameters")
      CALL hdf5_create_group(IPMPO,"/parameters/info")
      CALL hdf5_create_group(IPMPO,"/parameters/tree")
      CALL hdf5_create_group(IPMPO,"/parameters/values")
      CALL hdf5_create_group(IPMPO,"/energymesh")
      CALL hdf5_create_group(IPMPO,"/geometry")
      ICAL=0
      CALL hdf5_write_data(IPMPO,"/parameters/tree/NSTATEPOINT",ICAL)
*
* PRINT THE TITLE.
      IF(IMPX.GT.0) THEN
         CALL hdf5_read_data(IPMPO,"info/MPO_CREATION_INFO",HSMPO)
         CALL hdf5_read_data(IPMPO,"info/MPO_VERSION",IVERS)
         WRITE(6,400) HSMPO,IVERS
      ENDIF
*
* ADD THE TIME PARAMETER.
      DO 170 I=1,NPAR
      IF((PARTYP(I).EQ.'BURNUP').OR.(PARTYP(I).EQ.'FLUB')) GO TO 180
  170 CONTINUE
      GO TO 220
  180 DO 210 I=1,NPAR
      IF(PARTYP(I).EQ.'TIME') GO TO 220
  210 CONTINUE
      NPAR=NPAR+1
      IF(NPAR.GT.MAXPAR) CALL XABORT('MPO: TOO MANY PARAMETERS.')
      PARKEY(NPAR)='Time'
      PARTYP(NPAR)='TIME'
      PARFMT(NPAR)='FLOAT'
      PARADR(NPAR+1)=NPCHR
      NVALUE(NPAR)=0
*----
*  STORE THE MPO INITIALIZATION INFORMATION.
*----
  220 IF(NISO.GT.0) THEN
         CALL COMISO(0,MAXISO,IPLB1,NISO-1,NOMISO,NOMEVO,TYPISO)
         CALL hdf5_write_data(IPMPO,"/contents/isotopes/ISOTOPENAME",
     1   NOMISO(:NISO))
      ENDIF
      IF(NREA.GT.0) THEN
         CALL hdf5_write_data(IPMPO,"/contents/reactions/REACTIONAME",
     1   NOMREA(:NREA))
      ENDIF
*
      IF(NPAR.GT.0) THEN
         CALL hdf5_write_data(IPMPO,"/parameters/info/PARAMNAME",
     1   PARKEY(:NPAR))
         CALL hdf5_write_data(IPMPO,"/parameters/info/PARAMTYPE",
     1   PARTYP(:NPAR))
         CALL hdf5_write_data(IPMPO,"/parameters/info/PARAMFORM",
     1   PARFMT(:NPAR))
         IF(NPCHR.GT.0) THEN
           CALL hdf5_write_data(IPMPO,"/parameters/info/PARAMINFO",
     1     PARCAD(:NPCHR))
         ENDIF
         CALL hdf5_write_data(IPMPO,"/parameters/info/PARAMINFOADR",
     1   PARADR(:NPAR+1))
         CALL hdf5_write_data(IPMPO,"/parameters/info/NVALUE",
     1   NVALUE(:NPAR))
      ENDIF
*
      IF(NLOC.GT.0) THEN
         CALL hdf5_write_data(IPMPO,"/local_values/LOCVALNAME",
     1   PARKEL(:NLOC))
         CALL hdf5_write_data(IPMPO,"/local_values/LOCVALTYPE",
     1   PARTYL(:NLOC))
         CALL hdf5_write_data(IPMPO,"/local_values/LOCVALINFOADR",
     1   PARADL(:NLOC+1))
         CALL hdf5_write_data(IPMPO,"/local_values/NLOCVALINFO",NPCHL)
         IF(NPCHL.GT.0) THEN
           CALL hdf5_write_data(IPMPO,"/local_values/LOCVALINFO",
     1     PARCAL(:NPCHL))
         ENDIF
      ENDIF
      IF(NLOC.GT.0) DEALLOCATE(PARCAL,PARADL,PARTYL,PARKEL)
      IF(NPAR.GT.0) DEALLOCATE(NVALUE,PARCAD,PARADR,PARTYP,PARFMT,
     1 PARKEY)
      RETURN
* END OF MPO INITIALIZATION. **********************************
*----
*  INPUT AN ELEMENTARY CALCULATION. *******************************
*----
  300 CALL hdf5_read_data(IPMPO,"/info/MPO_VERSION",IVERS)
      IF(IVERS.NE.1) CALL XABORT('MPO: INVALID VERSION OF MPO SPECIF'//
     1 'ICATION.')
      NPAR=0
      IF(hdf5_group_exists(IPMPO,"/parameters/info/")) THEN
        CALL hdf5_info(IPMPO,"/parameters/info/NVALUE",RANK,TYPE,NBYTE,
     1  DIMSR)
        IF(RANK.GT.0) NPAR=DIMSR(1)
      ENDIF
      IF(NPAR.GT.0) THEN
         CALL hdf5_read_data(IPMPO,"/parameters/info/PARAMNAME",PARKEY)
         CALL hdf5_read_data(IPMPO,"/parameters/info/PARAMFORM",PARFMT)
         CALL hdf5_read_data(IPMPO,"/parameters/info/PARAMTYPE",PARTYP)
      ENDIF
*
      ITIM=0
      LWARN=.FALSE.
      HEDIT='output_0'
  310 CALL REDGET(INDIC,NITMA,FLOTT,TEXT24,DFLOTT)
      IF(INDIC.EQ.10) GO TO 350
      IF(INDIC.NE.3) CALL XABORT('MPO: CHARACTER DATA EXPECTED(18).')
      IF(TEXT24.EQ.'EDIT') THEN
*       READ THE PRINT INDEX.
        CALL REDGET(INDIC,IMPX,FLOTT,TEXT24,DFLOTT)
        IF(INDIC.NE.1) CALL XABORT('MPO: INTEGER DATA EXPECTED(5).')
      ELSE IF(TEXT24.EQ.'STEP') THEN
*       SET THE NAME OF THE OUTPUT SET.
        CALL REDGET(INDIC,NITMA,FLOTT,HEDIT,DFLOTT)
        IF(INDIC.NE.3) CALL XABORT('MPO: DIR-NAME EXPECTED.')
        IF(HEDIT(:7).NE.'output_') CALL XABORT('MPO: output_ EXPECTED.')
        IF(IMPX.GT.0) WRITE(6,'(/28H MPO: ACCESS A GROUP NAMED '',A,
     1  31H'' TO STORE THE MPO INFORMATION.)') TRIM(HEDIT)
      ELSE IF(TEXT24.EQ.'SET') THEN
        CALL REDGET(INDIC,NITMA,XT,TEXT24,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('MPO: REAL DATA EXPECTED(1).')
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT24,DFLOTT)
        IF(INDIC.NE.3) CALL XABORT('MPO: CHARACTER DATA EXPECTED'
     1  //'(19).')
        IF(TEXT24.EQ.'S') THEN
          XT=XT*1.0E-8
        ELSE IF(TEXT24.EQ.'DAY') THEN
          XT=XT*8.64E-4
        ELSE IF(TEXT24.EQ.'YEAR') THEN
          XT=XT*3.1536E-1
        ELSE
          CALL XABORT('MPO: S, DAY OR YEAR EXPECTED.')
        ENDIF
        IF(.NOT.C_ASSOCIATED(IPDEPL)) CALL XABORT('MPO: DEPLETION OBJ'
     1  //'ECT EXPECTED AT RHS.')
        CALL LCMLEN(IPDEPL,'DEPL-TIMES',NTIM,ITYLCM)
        IF(NTIM.EQ.0) CALL XABORT('MPO: NO DEPLETION TIME STEPS.')
        ALLOCATE(TIMES(NTIM))
        CALL LCMGET(IPDEPL,'DEPL-TIMES',TIMES)
        DO 320 I=1,NTIM
        IF(ABS(TIMES(I)-XT).LE.1.0E-4*XT) ITIM=I
  320   CONTINUE
        IF(ITIM.EQ.0) THEN
          WRITE(HSMG,'(39HMPO: UNABLE TO FIND A DEPLETION DIRECTO,
     1    12HRY AT TIME =,1P,E12.4,5H DAY.)') XT/8.64E-4
          CALL XABORT(HSMG)
        ENDIF
        DEALLOCATE(TIMES)
        IF(IMPX.GT.0) THEN
          WRITE(TEXT12,'(8HDEPL-DAT,I4.4)') ITIM
          WRITE(6,430) XT,XT/8.64E-4,TEXT12
        ENDIF
      ELSE IF(TEXT24.EQ.';') THEN
        GO TO 350
      ELSE IF(TEXT24.EQ.'WARNING-ONLY') THEN
        LWARN=.TRUE.
      ELSE
        DO 330 IKEY=1,NPAR
        IF(TEXT24.EQ.PARKEY(IKEY)) THEN
          IPAR=IKEY
          GO TO 340
        ENDIF
  330   CONTINUE
        CALL XABORT('MPO: INVALID KEYWORD='//TEXT24//'(5).')
  340   IF(PARTYP(IPAR).NE.'VALU') CALL XABORT('MPO: '//TEXT24//
     1  ' IS NOT OF VALU TYPE.')
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT24,DFLOTT)
        IF(PARFMT(IPAR).EQ.'INTEGER') THEN
          IF(INDIC.NE.1) CALL XABORT('MPO: INTEGER DATA EXPECTE'//
     1    'D(7).')
          IF(IMPX.GT.0) WRITE(6,450) TRIM(PARKEY(IPAR)),NITMA
        ELSE IF(PARFMT(IPAR).EQ.'FLOAT') THEN
          IF(INDIC.NE.2) CALL XABORT('MPO: REAL DATA EXPECTED(2).')
          IF(IMPX.GT.0) WRITE(6,440) TRIM(PARKEY(IPAR)),FLOTT
        ELSE IF(PARFMT(IPAR).EQ.'STRING') THEN
          IF(INDIC.NE.3) CALL XABORT('MPO: CHARACTER DATA EXPEC'//
     1    'TED(20).')
          IF(IMPX.GT.0) WRITE(6,460) TRIM(PARKEY(IPAR)),TEXT24
        ENDIF
        CALL MPOPAV(IPMPO,HEDIT,IPAR,NPAR,PARFMT(IPAR),FLOTT,NITMA,
     1  TEXT24,MUPLET(IPAR),LGNEW(IPAR))
      ENDIF
      GO TO 310
*----
*  RECOVER AN ELEMENTARY CALCULATION FROM EDITION.
*----
  350 NCALS=0
      IF(hdf5_group_exists(IPMPO,"/parameters/tree")) THEN
        CALL hdf5_read_data(IPMPO,"/parameters/tree/NSTATEPOINT",NCALS)
      ENDIF
      READ(HEDIT,'(7X,I2)') ID
      IF(C_ASSOCIATED(IPRHS)) GO TO 370 ! concatenation
      IF(hdf5_group_exists(IPMPO,"/output/")) THEN
        CALL hdf5_read_data(IPMPO,"/output/NOUTPUT",NOUTPUT)
        CALL hdf5_read_data(IPMPO,"/output/OUPUTID",OUPUTID)
        CALL hdf5_read_data(IPMPO,"/energymesh/NENERGYMESH",NENERG)
        CALL hdf5_read_data(IPMPO,"/geometry/NGEOMETRY",NGEOME)
        DO I=1,NGEOME
          DO J=1,NENERG
            IF(OUPUTID(J,I).EQ.ID) GO TO 360
          ENDDO
        ENDDO
        ALLOCATE(OUPUTID2(NENERG+1,NGEOME+1))
        OUPUTID2(:NENERG+1,:NGEOME+1)=0
        OUPUTID2(:NENERG,:NGEOME)=OUPUTID(:NENERG,:NGEOME)
        CALL hdf5_delete(IPMPO,"/energymesh/NENERGYMESH")
        CALL hdf5_delete(IPMPO,"/geometry/NGEOMETRY")
        CALL hdf5_delete(IPMPO,"/output/OUPUTID")
        DEALLOCATE(OUPUTID)
      ELSE
        CALL hdf5_create_group(IPMPO,"/output")
        ALLOCATE(OUPUTID2(1,1))
        OUPUTID2(1,1)=ID
        NENERG=0
        NGEOME=0
        NOUTPUT=0
      ENDIF
      NENERG=NENERG+1
      NGEOME=NGEOME+1
      CALL hdf5_write_data(IPMPO,"/energymesh/NENERGYMESH",NENERG)
      CALL hdf5_write_data(IPMPO,"/geometry/NGEOMETRY",NGEOME)
      OUPUTID2(NENERG,NGEOME)=ID
      CALL hdf5_write_data(IPMPO,"/output/OUPUTID",OUPUTID2)
      DEALLOCATE(OUPUTID2)
      NOUTPUT=NOUTPUT+1
      CALL hdf5_write_data(IPMPO,"/output/NOUTPUT",NOUTPUT)
*----
*  RECOVER ENERGY GROUP AND VOLUME INFORMATION.
*----
      CALL LCMGTC(IPEDIT,'LAST-EDIT',12,1,CDIRO)
      CALL LCMSIX(IPEDIT,CDIRO,1)
        CALL LCMSIX(IPEDIT,'MACROLIB',1)
          CALL LCMGET(IPEDIT,'STATE-VECTOR',ISTATE)
          NG=ISTATE(1)
          NMIL=ISTATE(2)
          WRITE(RECNAM,'(23H/energymesh/energymesh_,I0,1H/)') NENERG-1
          CALL hdf5_create_group(IPMPO,TRIM(RECNAM))
          CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"NG",NG)
          ALLOCATE(ENERG(NG+1))
          CALL LCMLEN(IPEDIT,'ENERGY',ILONG,ITYLCM)
          IF(ILONG.NE.NG+1) CALL XABORT('MPO: BAD VALUE OF NG.')
          CALL LCMGET(IPEDIT,'ENERGY',ENERG)
          WRITE(RECNAM,'(23H/energymesh/energymesh_,I0,1H/)') NENERG-1
          DO IGR=1,NG+1
            ENERG(IGR)=ENERG(IGR)*1.0E-6
          ENDDO
          CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"ENERGY",ENERG)
          DEALLOCATE(ENERG)
          WRITE(RECNAM,'(19H/geometry/geometry_,I0,1H/)') NGEOME-1
          CALL hdf5_create_group(IPMPO,TRIM(RECNAM))
          CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"NZONE",NMIL)
          ALLOCATE(VOL(NMIL))
          CALL LCMGET(IPEDIT,'VOLUME',VOL)
          WRITE(RECNAM,'(18Hgeometry/geometry_,I0,1H/)') NGEOME-1
          CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"ZONEVOLUME",VOL)
          DEALLOCATE(VOL)
        CALL LCMSIX(IPEDIT,' ',2)
      CALL LCMSIX(IPEDIT,' ',2)
*----
*  CREATE DUMMY DATASETS REACTION, ADDRISO AND ISOTOPE
*----
      CALL hdf5_get_shape(IPMPO,"/contents/isotopes/ISOTOPENAME",
     1 DIMS_MPO)
      NISO=DIMS_MPO(1)
      DEALLOCATE(DIMS_MPO)
      CALL hdf5_get_shape(IPMPO,"/contents/reactions/REACTIONAME",
     1 DIMS_MPO)
      NREA=DIMS_MPO(1)
      DEALLOCATE(DIMS_MPO)
      ALLOCATE(REACTION(NREA),ADDRISO(2),ISOTOPE(NISO))
      ADDRISO(1)=0
      ADDRISO(2)=NISO
      DO I=1,NREA
        REACTION(I)=I-1
      ENDDO
      DO I=1,NISO
        ISOTOPE(I)=I-1
      ENDDO
      WRITE(RECNAM,'(8H/output/,A)') TRIM(HEDIT)
      CALL hdf5_create_group(IPMPO,TRIM(RECNAM))
      WRITE(RECNAM,'(8H/output/,A,6H/info/)') TRIM(HEDIT)
      CALL hdf5_create_group(IPMPO,TRIM(RECNAM))
      CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"REACTION",REACTION)
      CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"ISOTOPE",ISOTOPE)
      CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"ADDRISO",ADDRISO)
      CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"NISO",NISO)
      CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"NREA",NREA)
      DEALLOCATE(ISOTOPE,ADDRISO,REACTION)
*----
*  RECOVER CALCULATION.
*----
  360 IF(IMPX.GT.0) WRITE(6,420) NCALS+1
      IF(ITIM.GT.0) THEN
        WRITE(TEXT12,'(8HDEPL-DAT,I4.4)') ITIM
        CALL LCMSIX(IPDEPL,TEXT12,1)
      ENDIF
*     -------------------------------------------
      CALL MPOCAL(IMPX,IPMPO,IPDEPL,IPEDIT,HEDIT)
*     -------------------------------------------
      IF(ITIM.GT.0) CALL LCMSIX(IPDEPL,' ',2)
      NG=0
      CALL MPOTOC(IPMPO,HEDIT,IMPX,NREA,NBISO,NMIL,NPAR,NLOC,NISOF,
     1 NISOP,NISOS,NCALS,NG,NSURFD,NALBP,NPRC)
*----
*  RECOVER REMAINING GLOBAL PARAMETER AND LOCAL VALUES.
*----
      NCALS=NCALS+1
      CALL hdf5_delete(IPMPO,"/parameters/tree/NSTATEPOINT")
      CALL hdf5_write_data(IPMPO,"/parameters/tree/NSTATEPOINT",NCALS)
      CALL MPOGEP(IPMPO,IPDEPL,IPLB1,IPLB2,IPEDIT,HEDIT,IMPX,ITIM,NPAR,
     1 NLOC,MUPLET,LGNEW,NMIL,NG,NCALS)
      IF(NPAR.GT.0) DEALLOCATE(PARTYP,PARFMT,PARKEY)
      RETURN
*----
*  MPO CONCATENATION.
*----
  370 NG=0
      CALL MPOTOC(IPRHS,HEDIT,IMPX,NREA,NBISO,NMIL,NPARR,NLOC,NISOF,
     1 NISOP,NISOS,NCALR,NG,NSURFD,NALBP,NPRC)
      IF(IMPX.GT.0) WRITE(6,470) NCALS+1,NCALS+NCALR
*     ------------------------------------------------------
      CALL MPOCAT(IPMPO,IPRHS,NPAR,MUPLET,LGNEW,LWARN,HEDIT)
*     ------------------------------------------------------
      IF(IMPX.GT.0) THEN
        CALL MPOTOC(IPMPO,HEDIT,IMPX,NREA,NBISO,NMIL,NPAR,NLOC,NISOF,
     1  NISOP,NISOS,NCALS,NG,NSURFD,NALBP,NPRC)
      ENDIF
      RETURN
*
  400 FORMAT(/6H MPO: ,A/6X,8HVERSION=,I3)
  420 FORMAT(/1X,43(1H*)/34H * MPO: ELEMENTARY CALCULATION NB.,I8,
     1 2H */1X,43(1H*))
  430 FORMAT(/41H MPO: RECOVER INFORMATION RELATED TO TIME,1P,E12.4,
     1 8H E+8 S (,E12.4,32H DAY) FROM LCM DIRECTORY NAMED ',A12,2H'.)
  440 FORMAT(28H MPO: SET GLOBAL PARAMETER ',A,3H' =,1P,E12.4)
  450 FORMAT(28H MPO: SET GLOBAL PARAMETER ',A,3H' =,I10)
  460 FORMAT(28H MPO: SET GLOBAL PARAMETER ',A,5H' = ',A,1H')
  470 FORMAT(/1X,55(1H*)/35H * MPO: ELEMENTARY CALCULATIONS NB.,I8,
     1 3H TO,I8,2H */1X,55(1H*))
      END
