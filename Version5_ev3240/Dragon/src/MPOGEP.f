*DECK MPOGEP
      SUBROUTINE MPOGEP(IPMPO,IPDEPL,IPLB1,IPLB2,IPEDIT,HEDIT,IMPX,
     1 ITIM,NPAR,NLOC,MUPLET,LGNEW,NMIL,NG,NCALAR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To recover remaining global parameters and local values. Update the
* parameter tree for a new elementary calculation.
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
* IPMPO   pointer to the MPO file.
* IPDEPL  pointer to the burnup object.
* IPLB1   pointer to the first microlib object.
* IPLB2   pointer to the second (optional) microlib object.
* IPEDIT  pointer to the edition object.
* HEDIT   name of output group for a (multigroup mesh, output geometry)
*         couple (generally equal to 'output_0').
* IMPX    print parameter.
* ITIM    index of the current burnup step.
* NPAR    number of global parameters.
* NLOC    number of local parameters.
* MUPLET  tuple of indices associated to each global parameter of the
*         elementary calculation.
* LGNEW   parameter modification flag (.TRUE. only if the I-th global
*         parameter has changed in the new elementary calculation).
* NMIL    number of mixtures in the MPO file
* NG      number of energy groups.
* NCALAR  index of the new elementary calculation.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE hdf5_wrap
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMPO,IPDEPL,IPLB1,IPLB2,IPEDIT
      INTEGER IMPX,ITIM,NPAR,NLOC,MUPLET(NPAR),NMIL,NG,NCALAR
      LOGICAL LGNEW(NPAR)
      CHARACTER(LEN=12) HEDIT
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) IPLB3
      PARAMETER (MAXPAR=50,NSTATE=40)
      INTEGER ISTATE(NSTATE)
      CHARACTER TEXT8*8,TEXT12*12,NAMLCM*12,NAMMY*12,HSMG*131,RECNAM*80
      LOGICAL EMPTY,LCM
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: PARADR,PARADL,LOCADR,
     1 DIMS_MPO
      REAL, ALLOCATABLE, DIMENSION(:) :: RVALO
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: PARFMT
      CHARACTER(LEN=24), ALLOCATABLE, DIMENSION(:) :: PARTYP,PARKEY,
     1 PARCAD,PARTYL,PARKEL,PARCAL
*----
*  VALIDATE NPAR
*----
      IF(NPAR.EQ.0) GO TO 45
      CALL hdf5_get_shape(IPMPO,"/parameters/info/PARAMNAME",DIMS_MPO)
      IF(NPAR.NE.DIMS_MPO(1)) CALL XABORT('MPOGEP: INVALID NPAR.')
      DEALLOCATE(DIMS_MPO)
*----
*  RECOVER INFORMATION FROM THE /parameters GROUP.
*----
      CALL hdf5_read_data(IPMPO,"/parameters/info/PARAMNAME",PARKEY)
      CALL hdf5_read_data(IPMPO,"/parameters/info/PARAMFORM",PARFMT)
      CALL hdf5_read_data(IPMPO,"/parameters/info/PARAMTYPE",PARTYP)
      CALL hdf5_read_data(IPMPO,"/parameters/info/PARAMINFOADR",PARADR)
      NPCHR=PARADR(NPAR+1)
      IF(NPCHR.GT.0) THEN
        CALL hdf5_read_data(IPMPO,"/parameters/info/PARAMINFO",PARCAD)
      ENDIF
*----
*  RECOVER REMAINING GLOBAL PARAMETERS.
*----
      DO 10 IPAR=1,NPAR
      IF(PARTYP(IPAR).EQ.'VALU') THEN
         GO TO 10
      ELSE IF((PARTYP(IPAR).EQ.'BURNUP').OR.(PARTYP(IPAR).EQ.'TIME').OR.
     1        (PARTYP(IPAR).EQ.'PUIS').OR.(PARTYP(IPAR).EQ.'FLUB').OR.
     2        (PARTYP(IPAR).EQ.'FLUX').OR.(PARTYP(IPAR).EQ.'MASL')) THEN
*
*        RECOVER GLOBAL PARAMETER VALUES FROM THE DEPLETION OBJECT.
         IF(.NOT.C_ASSOCIATED(IPDEPL)) CALL XABORT('MPOGEP: NO DEPLETI'
     1   //'ON OBJECT AVAILABLE AMONG THE RHS LCM OBJECTS.')
         CALL LCMGET(IPDEPL,'STATE-VECTOR',ISTATE)
         NBURN=ISTATE(3)
         NBISO=ISTATE(4)
         NREAC=ISTATE(6)
         NVAR=ISTATE(7)
         NBMIX=ISTATE(8)
         CALL COMGEM(IPDEPL,ITIM,PARTYP(IPAR),0,NBURN,NBMIX,NBISO,
     1   NREAC,NVAR,VALPAR)
      ELSE IF((PARTYP(IPAR).EQ.'TEMP').OR.(PARTYP(IPAR).EQ.'CONC'))
     1   THEN
*
*        RECOVER GLOBAL PARAMETER VALUES FROM A MICROLIB OBJECT.
         IF(.NOT.C_ASSOCIATED(IPLB1)) CALL XABORT('MPOGEP: MICROLIB EX'
     1   //'PECTED AT RHS.')
         IF(NPCHR.EQ.0) CALL XABORT('MPOGEP: MISSING PARAMINFO.')
         TEXT8=' '
         TEXT12=' '
         IMILI=0
         IPCHR=PARADR(IPAR)+1
         IF(PARTYP(IPAR).EQ.'CONC') THEN
           TEXT8=PARCAD(IPCHR)(:8)
           IPCHR=IPCHR+1
         ENDIF
         TEXT12=PARCAD(IPCHR)(:8)
         IPCHR=IPCHR+1
         READ(PARCAD(IPCHR),'(3X,I9)') IMILI
         CALL LCMGET(IPLB1,'STATE-VECTOR',ISTATE)
         MAXNBI=ISTATE(2)
         IF(C_ASSOCIATED(IPLB2)) THEN
            CALL LCMGET(IPLB2,'STATE-VECTOR',ISTATE)
            MAXNBI=MAX(MAXNBI,ISTATE(2))
         ENDIF
         CALL COMBIB(IPLB1,IPLB2,PARTYP(IPAR),IMILI,TEXT12,TEXT8,MAXNBI,
     1   VALPAR)
      ELSE
         CALL XABORT('MPOGEP: '//PARTYP(IPAR)//' IS AN UNKNOWN PARAM'//
     1   'ETER TYPE.')
      ENDIF
      IF(IMPX.GT.0) WRITE(6,100) PARKEY(IPAR),VALPAR
*
      CALL MPOPAV(IPMPO,HEDIT,IPAR,NPAR,PARFMT(IPAR),VALPAR,NITMA,
     1 TEXT12,MUPLET(IPAR),LGNEW(IPAR))
   10 CONTINUE
      IF(IMPX.GT.2) THEN
         WRITE(6,110) (MUPLET(I),I=1,NPAR)
         WRITE(6,'(/)')
      ENDIF
      DO 15 I=1,NPAR
      IF(MUPLET(I).EQ.-99) THEN
         WRITE(HSMG,'(33HMPOGEP: UNDEFINED MUPLET ELEMENT=,I6)') I
         CALL XABORT(HSMG)
      ENDIF
   15 CONTINUE
      WRITE(RECNAM,'(8H/output/,A,9H/statept_,I0)') TRIM(HEDIT),NCALAR-1
      CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"/PARAMVALUEORD",MUPLET)
      IF(NPCHR.GT.0) DEALLOCATE(PARCAD)
      DEALLOCATE(PARADR,PARTYP,PARFMT,PARKEY)
*----
*  RECOVER INFORMATION FROM THE 'varlocdescri' GROUP.
*----
   45 IF(NLOC.EQ.0) RETURN
      ALLOCATE(LOCADR(NLOC+1))
      CALL hdf5_read_data(IPMPO,"/local_values/LOCVALNAME",PARKEL)
      CALL hdf5_read_data(IPMPO,"/local_values/LOCVALTYPE",PARTYL)
      CALL hdf5_read_data(IPMPO,"/local_values/LOCVALINFOADR",PARADL)
      CALL hdf5_read_data(IPMPO,"/local_values/NLOCVALINFO",NPCHL)
      IF(NPCHL.GT.0) THEN
        CALL hdf5_read_data(IPMPO,"/local_values/LOCVALINFO",PARCAL)
      ENDIF
*
      CALL LCMGTC(IPEDIT,'LAST-EDIT',12,1,TEXT12)
*----
*  INITIALIZE LOCADR AND ALLOCATE RVALOC.
*----
      IADR=0
      LOCADR(1)=0
      DO 50 IPAR=1,NLOC
      IF((PARTYL(IPAR).EQ.'EQUI').OR.(PARTYL(IPAR).EQ.'VITE')) THEN
         IADR=IADR+NG
      ELSE IF(PARTYL(IPAR).EQ.'COUR') THEN
         IADR=IADR+2*NG
      ELSE
         IADR=IADR+1
      ENDIF
      LOCADR(IPAR+1)=IADR
   50 CONTINUE
      NVLC=LOCADR(NLOC+1)
      ALLOCATE(RVALO(NVLC*NMIL))
*----
*  RECOVER LOCAL VARIABLES.
*----
      DO 70 IPAR=1,NLOC
      IF((PARTYL(IPAR).EQ.'BURNUP').OR.(PARTYL(IPAR).EQ.'TIME').OR.
     1   (PARTYL(IPAR).EQ.'PUIS').OR.(PARTYL(IPAR).EQ.'FLUG').OR.
     2   (PARTYL(IPAR).EQ.'FLUB').OR.(PARTYL(IPAR).EQ.'FLUX').OR.
     3   (PARTYL(IPAR).EQ.'MASL')) THEN
*
*        RECOVER LOCAL VARIABLES FROM THE DEPLETION OBJECT.
         IF(.NOT.C_ASSOCIATED(IPDEPL)) CALL XABORT('MPOGEP: NO DEPLET'
     1   //'ION OBJECT AVAILABLE AMONG THE RHS LCM OBJECTS.')
         CALL LCMGET(IPDEPL,'STATE-VECTOR',ISTATE)
         NBURN=ISTATE(3)
         NBISO=ISTATE(4)
         NREAC=ISTATE(6)
         NVAR=ISTATE(7)
         NBMIX=ISTATE(8)
         CALL LCMGET(IPEDIT,'STATE-VECTOR',ISTATE)
         NREG=ISTATE(17)
         CALL COMGEN(IPDEPL,IPEDIT,NREG,NMIL,ITIM,PARTYL(IPAR),NBURN,
     1   NBMIX,NBISO,NREAC,NVAR,LOCADR(IPAR),NVLC,RVALO)
      ELSE IF((PARTYL(IPAR).EQ.'TEMP').OR.(PARTYL(IPAR).EQ.'CONC'))
     1   THEN
*
*        RECOVER LOCAL VARIABLES FROM THE MICROLIB IN EDIT OBJECT.
         TEXT8=' '
         IF(PARTYL(IPAR).EQ.'CONC') THEN
           IF(NPCHL.EQ.0) CALL XABORT('MPOGEP: MISSING LOCVALINFO.')
           IPCHL=PARADL(IPAR)+1
           TEXT8=PARCAL(IPCHL)(:8)
         ENDIF
         CALL LCMSIX(IPEDIT,TEXT12,1)
         CALL LCMGET(IPEDIT,'STATE-VECTOR',ISTATE)
         MAXNBI=ISTATE(2)
         CALL LCMINF(IPEDIT,NAMLCM,NAMMY,EMPTY,ILONG,LCM)
         IPLB3=C_NULL_PTR
         DO 60 IBM=1,NMIL
         CALL COMBIB(IPEDIT,IPLB3,PARTYL(IPAR),IBM,NAMLCM,TEXT8,MAXNBI,
     1   VALPAR)
         RVALO((IBM-1)*NVLC+LOCADR(IPAR))=VALPAR
   60    CONTINUE
         CALL LCMSIX(IPEDIT,' ',2)
      ELSE IF(PARTYL(IPAR).EQ.'EQUI') THEN
*        RECOVER A SET OF SPH EQUIVALENCE FACTORS.
         CALL SAPSPH(IPEDIT,NG,NMIL,LOCADR(IPAR),NVLC,RVALO)
      ELSE
         CALL XABORT('MPOGEP: '//PARTYL(IPAR)//' IS AN UNKNOWN LOCAL'//
     1   ' VARIABLE TYPE.')
      ENDIF
      IF(IMPX.GT.1) WRITE(6,120) PARKEY(IPAR),
     1            (RVALO((IBM-1)*NVLC+LOCADR(IPAR)),IBM=1,NMIL)
   70 CONTINUE
      DO 80 IBM=1,NMIL
      WRITE(RECNAM,'(8H/output/,A,9H/statept_,I0,6H/zone_,I0,1H/)')
     1 TRIM(HEDIT),NCALAR-1,IBM-1
      CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"LOCALVALUE",
     1 RVALO((IBM-1)*NVLC+1:IBM*NVLC))
      CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"LOCALVALADDR",LOCADR)
   80 CONTINUE
      DEALLOCATE(RVALO)
      IF(NPCHL.GT.0) DEALLOCATE(PARCAL)
      DEALLOCATE(PARADL,PARTYL,PARKEL,LOCADR)
      RETURN
*
  100 FORMAT(31H MPOGEP: SET GLOBAL PARAMETER ',A,3H' =,1P,E12.4)
  110 FORMAT(/16H MPOGEP: MUPLET=,10I6:/(16X,10I6))
  120 FORMAT(29H MPOGEP: SET LOCAL VARIABLE ',A,3H' =,1P,5E12.4/(36X,
     1 5E12.4))
      END
