*DECK MPOCAL
      SUBROUTINE MPOCAL(IMPX,IPMPO,IPDEPL,IPEDIT,HEDIT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Store the results of an elementary calculation in the MPO file.
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
* IMPX    print parameter.
* IPMPO   pointer to the MPO file.
* IPDEPL  pointer to the burnup object (L_BURNUP signature).
* IPEDIT  pointer to the edition object (L_EDIT signature).
* HEDIT   name of output group for a (multigroup mesh, output geometry)
*         couple (generally equal to 'output_0').
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE hdf5_wrap
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMPO,IPDEPL,IPEDIT,IPSPH
      INTEGER IMPX
      CHARACTER(LEN=12) HEDIT
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      INTEGER IPAR(NSTATE),RANK,TYPE,NBYTE,DIMSR(5)
      REAL BIRRAD(2)
      CHARACTER CDIRO*12,HSMG*131,RECNAM*80
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: DIMS_MPO
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: OUPUTID
      REAL, ALLOCATABLE, DIMENSION(:) :: VOLMIL
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FLXMIL
*----
*  RECOVER MICROLIB AND MACROLIB INFORMATION
*----
      CALL LCMGTC(IPEDIT,'LAST-EDIT',12,1,CDIRO)
      CALL LCMSIX(IPEDIT,CDIRO,1)
      CALL LCMLEN(IPEDIT,'STATE-VECTOR',ILONG,ITYLCM)
      IF(ILONG.NE.0) THEN
         CALL LCMGET(IPEDIT,'STATE-VECTOR',IPAR)
         NBISO=IPAR(2)
         NED=IPAR(13)
         NPRC=IPAR(19)
         NDFI=IPAR(20)
      ELSE
         NBISO=0
         NDFI=0
      ENDIF
      CALL LCMSIX(IPEDIT,'MACROLIB',1)
        CALL LCMGET(IPEDIT,'STATE-VECTOR',IPAR)
        NG=IPAR(1)
        NMIL=IPAR(2)
        NL=IPAR(3)
        IF(IPAR(4).GT.1) CALL XABORT('MPOCAL: CANNOT PROCESS MULTIPLE '
     1  //'FISSION SPECTRA.')
        NED=IPAR(5)
        ITRANC=IPAR(6)
        NPRC=IPAR(7)
        NALBP=IPAR(8)
        ILEAK=IPAR(9)
        IDF=IPAR(12)
        CALL LCMLEN(IPEDIT,'SPH',ILEN,ITYLCM)
        IF(ILEN.NE.0) THEN
          IPSPH=LCMGID(IPEDIT,'SPH')
          CALL LCMGET(IPSPH,'STATE-VECTOR',IPAR)
          IMC=IPAR(6)
        ELSE
          IMC=0
        ENDIF
      CALL LCMSIX(IPEDIT,' ',2)
*----
*  RECOVER ENERGY ID_G AND ID_E
*----
      CALL hdf5_read_data(IPMPO,"/output/NOUTPUT",NOUTPUT)
      ID_G=-1
      ID_E=-1
      CALL hdf5_read_data(IPMPO,"/output/OUPUTID",OUPUTID)
      CALL hdf5_read_data(IPMPO,"/energymesh/NENERGYMESH",NENERG)
      CALL hdf5_read_data(IPMPO,"/geometry/NGEOMETRY",NGEOME)
      READ(HEDIT,'(7X,I2)') ID
      DO I=1,NGEOME
        DO J=1,NENERG
          IF(OUPUTID(J,I).EQ.ID) THEN
            ID_G=I-1
            ID_E=J-1
            GO TO 10
          ENDIF
        ENDDO
      ENDDO
      CALL XABORT('MPOCAL: no ID found in /output/OUPUTID.')
   10 CALL hdf5_read_data(IPMPO,"/parameters/tree/NSTATEPOINT",NCALS)
      ICAL=NCALS+1
*----
*  RECOVER THE FLUX NORMALIZATION FACTOR.
*----
      IF(C_ASSOCIATED(IPDEPL)) THEN
         CALL LCMGET(IPDEPL,'BURNUP-IRRAD',BIRRAD)
         BURN=BIRRAD(1)
         CALL LCMLEN(IPDEPL,'FLUX-NORM',ILONG,ITYLCM)
         IF(ILONG.EQ.0) THEN
            WRITE(HSMG,'(40HMPOCAL: THE ''FLUX-NORM'' RECORD IS NOT SE,
     1      20HT FOR BURNUP STEP AT,E12.5,14H MW-DAY/TONNE.)') BURN
            CALL XABORT(HSMG)
         ENDIF
         CALL LCMGET(IPDEPL,'FLUX-NORM',FNORM)
         IF(IMPX.GT.0) WRITE(6,100) FNORM,BURN
      ELSE
         FNORM=1.0
         IF(IMPX.GT.0) WRITE(6,110)
      ENDIF
*----
*  RECOVER THE NUMBER OF ADRX (NADRX) AND TRANSPROFILE (NADRI) SETS.
*----
      WRITE(RECNAM,'(8H/output/,A,6H/info/)') TRIM(HEDIT)
      NADRX=0
      CALL hdf5_info(IPMPO,TRIM(RECNAM)//"ADDRXS",RANK,TYPE,NBYTE,DIMSR)
      IF(TYPE.NE.99) THEN
        NADRX=DIMSR(3)
      ENDIF
      NADRI=0
      CALL hdf5_info(IPMPO,TRIM(RECNAM)//"TRANSPROFILE",RANK,TYPE,NBYTE,
     1 DIMSR)
      IF(TYPE.NE.99) NADRI=DIMSR(1)/(2*NG+1)
*----
*  RECOVER THE CROSS SECTIONS.
*----
      CALL hdf5_get_shape(IPMPO,"/contents/isotopes/ISOTOPENAME",
     1 DIMS_MPO)
      NISO=DIMS_MPO(1)
      DEALLOCATE(DIMS_MPO)
      CALL hdf5_get_shape(IPMPO,"/contents/reactions/REACTIONAME",
     1 DIMS_MPO)
      NREA=DIMS_MPO(1)
      DEALLOCATE(DIMS_MPO)
      MAXRDA=(NREA*NG+NL*NG+NL*NG*NG)*NISO
      MAXIDA=(2*NG+1)*(NADRI+NISO*NMIL)
      ALLOCATE(VOLMIL(NMIL),FLXMIL(NMIL,NG))
      CALL MPOCA2(IPMPO,IPEDIT,HEDIT,NREA,NISO,NADRX,NED,NPRC,ILEAK,
     1 NG,NMIL,NL,ITRANC,NALBP,IMC,NBISO,ICAL,MAXRDA,MAXIDA,FNORM,IMPX,
     2 NISOTS,NISFS,NISPS,VOLMIL,FLXMIL)
*----
*  RECOVER DISCONTINUITY FACTOR INFORMATION.
*----
      IF((IDF.EQ.2).OR.(IDF.EQ.3).OR.(NALBP.GT.0)) THEN
        CALL MPOIDF(IPMPO,IPEDIT,HEDIT,NG,NMIL,ICAL,IDF,NALBP,FNORM,
     1  VOLMIL,FLXMIL)
      ENDIF
      DEALLOCATE(FLXMIL,VOLMIL)
*----
*  RECOVER THE FISSION YIELDS.
*----
      IF(NISFS*NISPS.GT.0) THEN
         CALL MPOGEY(IPMPO,IPEDIT,HEDIT,NISO,NG,NMIL,NBISO,ICAL,NDFI,
     1   NISFS,NISPS)
      ENDIF
*
      CALL LCMSIX(IPEDIT,' ',2)
      RETURN
*
  100 FORMAT(45H MPOCAL: NORMALIZE THE FLUX WITH THE FACTOR =,1P,E12.5,
     1 26H TAKEN FROM BURNUP STEP AT,E12.5,14H MW-DAY/TONNE.)
  110 FORMAT(36H MPOCAL: THE FLUX IS NOT NORMALIZED.)
      END
