*DECK MCRTRP
      SUBROUTINE MCRTRP(IPMPO,LCUB2,IMPX,NPAR,NCAL,MUPLET,MUTYPE,PARTYP,
     1 VALR,VARVAL,MUBASE,TERP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the TERP interpolation/derivation/integration factors using
* table-of-content information of the MPO file.
*
*Copyright:
* Copyright (C) 2022 Ecole Polytechnique de Montreal
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* IPMPO   address of the multidimensional MPO file.
* LCUB2   interpolation type for each parameter (=.TRUE.: cubic Ceschino
*         interpolation; =.FALSE: linear Lagrange interpolation).
* IMPX    print parameter (equal to zero for no print).
* NPAR    number of global parameters.
* NCAL    number of elementary calculations in the MPO file.
* MUPLET  tuple used to identify an elementary calculation.
* MUTYPE  type of interpolation (=1: interpolation; =2: delta-sigma).
* PARTYP  parameter types.
* VALR    real values of the interpolated point.
* VARVAL  exit burnup used if MUTYPE(IPAR(ID))=3.
* MUBASE  muplet database.
*
*Parameters: output
* TERP    interpolation factors.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE hdf5_wrap
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER, PARAMETER::MAXPAR=50
      TYPE(C_PTR) IPMPO
      INTEGER IMPX,NPAR,NCAL,MUPLET(NPAR),MUTYPE(NPAR),MUBASE(NPAR,NCAL)
      REAL VALR(2*MAXPAR,2),VARVAL,TERP(NCAL)
      LOGICAL LCUB2(NPAR)
      CHARACTER(LEN=24) PARTYP(NPAR)
*----
*  LOCAL VARIABLES
*----
      INTEGER, PARAMETER::IOUT=6
      INTEGER, PARAMETER::MAXDIM=10
      INTEGER, PARAMETER::MAXVAL=200
      INTEGER IPAR(MAXDIM),NVAL(MAXDIM),IDDIV(MAXDIM)
      REAL BURN0, BURN1, DENOM, TERTMP
      INTEGER I, ICAL, ID, IDTMP, IDTOT, JD, NDELTA, NDIM, NID, NTOT,
     1 MCRCAL, IBURN, ITIME
      REAL T1D(MAXVAL,MAXDIM),WORK(MAXVAL)
      CHARACTER HSMG*131,RECNAM*80
      LOGICAL LCUBIC,LSINGL
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NVALUE,MUPLE2
      REAL, ALLOCATABLE, DIMENSION(:) :: TERPA,VREAL
      CHARACTER(LEN=80), ALLOCATABLE, DIMENSION(:) :: PARNAM
*----
*  RECOVER TREE INFORMATION
*----
      IBURN=0
      ITIME=0
      IF(NPAR.GT.0) THEN
        CALL hdf5_read_data(IPMPO,"/parameters/info/NVALUE",NVALUE)
        DO I=1,NPAR
          IF(PARTYP(I).EQ.'BURNUP') IBURN=I
          IF(PARTYP(I).EQ.'TIME') ITIME=I
        ENDDO
      ENDIF
*----
*  COMPUTE TERP FACTORS
*----
      TERP(:NCAL)=0.0
      IPAR(:MAXDIM)=0
      NDIM=0
      NDELTA=0
      DO 10 I=1,NPAR
        IF(MUPLET(I).EQ.-1) THEN
          NDIM=NDIM+1
          IF(MUTYPE(I).NE.1) NDELTA=NDELTA+1
          IF(NDIM.GT.MAXDIM) THEN
            WRITE(HSMG,'(7HMCRTRP:,I4,29H-DIMENSIONAL INTERPOLATION NO,
     1      14HT IMPLEMENTED.)') NDIM
            CALL XABORT(HSMG)
          ENDIF
          IPAR(NDIM)=I
        ELSE IF((MUPLET(I).EQ.0).AND.(NVALUE(I).EQ.1)) THEN
          MUPLET(I)=1
        ENDIF
   10 CONTINUE
      IF(IMPX.GT.2) THEN
        WRITE(IOUT,'(16H MCRTRP: MUPLET=,10I4/(16X,10I4))')
     1  (MUPLET(I),I=1,NPAR)
        WRITE(IOUT,'(8H MCRTRP:,I4,31H-DIMENSIONAL INTERPOLATION IN M,
     1  3HPO.)') NDIM
      ENDIF
      ALLOCATE(MUPLE2(NPAR))
      IF(NDIM.EQ.0) THEN
        MUPLE2(:NPAR)=MUPLET(:NPAR)
        IF((MUPLET(IBURN).NE.0).AND.(MUPLET(ITIME).EQ.0)) THEN
          MUPLE2(ITIME)=MUPLE2(IBURN)
        ELSE IF((MUPLET(IBURN).EQ.0).AND.(MUPLET(ITIME).NE.0)) THEN
          MUPLE2(IBURN)=MUPLE2(ITIME)
        ENDIF
        ICAL=0
        IF(NPAR.GT.0) ICAL=MCRCAL(NPAR,NCAL,MUPLE2,MUBASE)
        IF(ICAL.GT.NCAL) CALL XABORT('MCRTRP: TERP OVERFLOW(1).')
        IF(ICAL.EQ.0) GO TO 200
        IF(ICAL.EQ.-1) GO TO 210
        TERP(ICAL)=1.0
      ELSE
        NTOT=1
        IDDIV(:MAXDIM)=1
        DO 70 ID=1,NDIM
        IF(IPAR(ID).LE.NPAR) THEN
          WRITE(RECNAM,'(25H/parameters/values/PARAM_,I0)') IPAR(ID)-1
          NID=NVALUE(IPAR(ID))
        ELSE
          CALL XABORT('MCRTRP: PARAMETER INDEX OVERFLOW.')
        ENDIF
        NTOT=NTOT*NID
        DO 15 IDTMP=1,NDIM-ID
        IDDIV(IDTMP)=IDDIV(IDTMP)*NID
   15   CONTINUE
        CALL hdf5_read_data(IPMPO,RECNAM,VREAL)
        BURN0=VALR(IPAR(ID),1)
        BURN1=VALR(IPAR(ID),2)
        LSINGL=(BURN0.EQ.BURN1)
        LCUBIC=LCUB2(IPAR(ID))
        IF((MUTYPE(IPAR(ID)).EQ.1).AND.LSINGL) THEN
          CALL ALTERP(LCUBIC,NID,VREAL,BURN0,.FALSE.,T1D(1,ID))
        ELSE IF(MUTYPE(IPAR(ID)).EQ.1) THEN
          IF(BURN0.GE.BURN1) CALL XABORT('MCRTRP: INVALID BURNUP'
     1     //' LIMITS(1).')
          CALL ALTERI(LCUBIC,NID,VREAL,BURN0,BURN1,T1D(1,ID))
          DO 20 I=1,NID
          T1D(I,ID)=T1D(I,ID)/(BURN1-BURN0)
   20     CONTINUE
        ELSE IF((MUTYPE(IPAR(ID)).EQ.2).AND.(.NOT.LSINGL)) THEN
          CALL ALTERP(LCUBIC,NID,VREAL,BURN0,.FALSE.,WORK)
          CALL ALTERP(LCUBIC,NID,VREAL,BURN1,.FALSE.,T1D(1,ID))
          DO 30 I=1,NID
          T1D(I,ID)=T1D(I,ID)-WORK(I)
   30     CONTINUE
        ELSE IF((MUTYPE(IPAR(ID)).EQ.2).AND.(LSINGL)) THEN
          T1D(:NID,ID)=0.0
        ELSE IF(MUTYPE(IPAR(ID)).EQ.3) THEN
*          DERIVATIVE WITH RESPECT TO A SINGLE EXIT BURNUP. USE
*          EQ.(3.3) OF RICHARD CHAMBON'S THESIS.
          IF(BURN0.GE.BURN1) CALL XABORT('MCRTRP: INVALID BURNUP'
     1    //' LIMITS(2).')
          CALL hdf5_read_data(IPMPO,"/paramdescrip/PARNAM",PARNAM)
          IF(PARNAM(IPAR(ID)).NE.'Burnup') THEN
            CALL XABORT('MCRTRP: Burnup EXPECTED.')
          ENDIF
          DEALLOCATE(PARNAM)
          ALLOCATE(TERPA(NID))
          CALL ALTERI(LCUBIC,NID,VREAL,BURN0,BURN1,TERPA)
          DO 40 I=1,NID
          T1D(I,ID)=-TERPA(I)
   40     CONTINUE
          CALL ALTERP(LCUBIC,NID,VREAL,BURN0,.FALSE.,TERPA)
          DO 50 I=1,NID
          T1D(I,ID)=T1D(I,ID)-TERPA(I)*BURN0
   50     CONTINUE
          CALL ALTERP(LCUBIC,NID,VREAL,BURN1,.FALSE.,TERPA)
          DENOM=VARVAL*(BURN1-BURN0)
          DO 60 I=1,NID
          T1D(I,ID)=(T1D(I,ID)+TERPA(I)*BURN1)/DENOM
   60     CONTINUE
          DEALLOCATE(TERPA)
        ELSE
          CALL XABORT('MCRTRP: INVALID OPTION.')
        ENDIF
        DEALLOCATE(VREAL)
        NVAL(ID)=NID
   70   CONTINUE

* Example: NDIM=3, NVALUE=(3,2,2)
* IDTOT 1  2  3  4  5  6  7  8  9 10 11 12
* ID(1) 1  2  3  1  2  3  1  2  3  1  2  3
* ID(2) 1  1  1  2  2  2  1  1  1  2  2  2
* ID(3) 1  1  1  1  1  1  2  2  2  2  2  2
* (NTOT=12, IDDIV=(6,3,1))
        DO 100 IDTOT=1,NTOT               ! Ex.: IDTOT       = 9
          TERTMP=1.0
          IDTMP=IDTOT
          DO 80 JD=1,NDIM                 ! Ex.: JD          = 1,2,3
            ID=(IDTMP-1)/IDDIV(JD)+1      ! Ex.: ID(NDIM...1)= 2,1,3
            IDTMP=IDTMP-(ID-1)*IDDIV(JD)  ! Ex.: IDTMP       = 3,3,1
            MUPLET(IPAR(NDIM-JD+1))=ID
            TERTMP=TERTMP*T1D(ID,NDIM-JD+1)
   80     CONTINUE
          MUPLE2(:NPAR)=MUPLET(:NPAR)
          IF((MUPLET(IBURN).NE.0).AND.(MUPLET(ITIME).EQ.0)) THEN
            MUPLE2(ITIME)=MUPLE2(IBURN)
          ELSE IF((MUPLET(IBURN).EQ.0).AND.(MUPLET(ITIME).NE.0)) THEN
            MUPLE2(IBURN)=MUPLE2(ITIME)
          ENDIF
          ICAL=MCRCAL(NPAR,NCAL,MUPLE2,MUBASE)
          IF(ICAL.GT.NCAL) CALL XABORT('MCRTRP: TERP OVERFLOW(2).')
          IF(ICAL.EQ.0) GO TO 200
          IF(ICAL.EQ.-1) GO TO 210
          TERP(ICAL)=TERP(ICAL)+TERTMP
  100   CONTINUE
      ENDIF
      IF(IMPX.GT.3) THEN
        WRITE(IOUT,'(25H MCRTRP: TERP PARAMETERS:/(1X,1P,10E12.4))')
     1  (TERP(I),I=1,NCAL)
      ENDIF
      DEALLOCATE(MUPLE2)
      IF(NPAR.GT.0) DEALLOCATE(NVALUE)
      RETURN
*----
*  MISSING ELEMENTARY CALCULATION EXCEPTION.
*----
  200 WRITE(IOUT,'(16H MCRTRP: MUPLET=,10I4/(16X,10I4))')
     1 (MUPLET(I),I=1,NPAR)
      CALL XABORT('MCRTRP: MISSING ELEMENTARY CALCULATION.')
  210 WRITE(IOUT,'(16H MCRTRP: MUPLET=,10I4/(16X,10I4))')
     1 (MUPLET(I),I=1,NPAR)
      WRITE(IOUT,'(9X,7HNVALUE=,10I4/(16X,10I4))') (NVALUE(I),I=1,NPAR)
      CALL XABORT('MCRTRP: DEGENERATE ELEMENTARY CALCULATION.')
      END
