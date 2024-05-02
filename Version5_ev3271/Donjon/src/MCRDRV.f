*DECK MCRDRV
      SUBROUTINE MCRDRV(IPMPO,LCUBIC,NMIX,IMPX,NMIL,NCAL,NBISO,NPAR,
     1 HEDIT,ITER,MAXNIS,MIXC,TERP,NISO,LISO,HISO,CONC,ITODO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute TERP factors for MPO file interpolation. Use user-defined
* global parameters.
*
*Copyright:
* Copyright (C) 2022 Ecole Polytechnique de Montreal
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* IPMPO   address of the MPO file.
* LCUBIC  =.TRUE.: cubic Ceschino interpolation; =.FALSE: linear
*         Lagrange interpolation.
* NMIX    maximum number of material mixtures in the microlib.
* IMPX    print parameter (equal to zero for no print).
* NMIL    number of material mixtures in the MPO file.
* NCAL    number of elementary calculations in the MPO file.
* NBISO   number of particularized isotopes in the MPO file.
* NPAR    number of parameters
* HEDIT   name of output group for a (multigroup mesh, output geometry)
*         couple (generally equal to 'output_0').
*
*Parameters: output
* ITER    completion flag (=0: all over; =1: use another MPO file;
*         =2 use another L_MAP + MPO file).
* MAXNIS  maximum value of NISO(I) in user data.
* MIXC    mixture index in the MPO file corresponding to each microlib
*         mixture.
* TERP    interpolation factors.
* NISO    number of user-selected isotopes.
* LISO    type of treatment (=.true.: ALL; =.false.: ONLY).
* HISO    name of the user-selected isotopes.
* CONC    user-defined number density of the user-selected isotopes. A
*         value of -99.99 is set to indicate that the compo value is
*         used.
* ITODO   non-depletion mask (=1 to force a user-selected isotope to be
*         non-depleting)
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE hdf5_wrap
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMPO
      INTEGER NMIX,IMPX,NMIL,NCAL,NBISO,NPAR,ITER,MAXNIS,MIXC(NMIX),
     1 NISO(NMIX),ITODO(NMIX,NBISO)
      REAL TERP(NCAL,NMIX),CONC(NMIX,NBISO)
      LOGICAL LCUBIC,LISO(NMIX)
      CHARACTER(LEN=8) HISO(NMIX,NBISO)
      CHARACTER(LEN=12) HEDIT
*----
*  LOCAL VARIABLES
*----
      INTEGER, PARAMETER::IOUT=6
      INTEGER, PARAMETER::MAXLIN=132
      INTEGER, PARAMETER::MAXPAR=50
      INTEGER, PARAMETER::MAXVAL=200
      REAL, PARAMETER::REPS=1.0E-4
      INTEGER I, J, IBM, IBMOLD, ICAL, INDIC, IPAR, ITYPE, JBM, NITMA
      REAL SUM, FLOTT
      CHARACTER TEXT72*72,HSMG*131,TEXT132*132,VALH(MAXPAR)*12,
     1 RECNAM*80,HCUBIC*12
      INTEGER VALI(MAXPAR)
      INTEGER RANK,TYPE,NBYTE,DIMSR(5)
      DOUBLE PRECISION DFLOTT
      REAL VALR(2*MAXPAR,2)
      LOGICAL LCUB2(MAXPAR)
*----
*  ALLOCATABLE ARRAYS
*----
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: LDELTA
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MUPLET,MUTYPE,NVALUE,VINTE
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: MUBASE
      REAL, ALLOCATABLE, DIMENSION(:) :: VREAL
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: PARFMT
      CHARACTER(LEN=12), ALLOCATABLE, DIMENSION(:) :: VCHAR
      CHARACTER(LEN=24), ALLOCATABLE, DIMENSION(:) :: PARTYP,PARKEY
      CHARACTER(LEN=132), ALLOCATABLE, DIMENSION(:) :: TEXT132V1
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(MUPLET(NPAR),MUTYPE(NPAR),LDELTA(NMIX),MUBASE(NPAR,NCAL))
*----
*  RECOVER INFORMATION FOR THE MPO FILE.
*----
      CALL hdf5_info(IPMPO,"/info/MPO_CREATION_INFO",RANK,TYPE,NBYTE,
     1 DIMSR)
      IF(RANK.GT.MAXLIN) CALL XABORT('MCRDRV: MAXLIN OVERFLOW.')
      IF(NPAR.GT.MAXPAR) CALL XABORT('MCRDRV: MAXPAR OVERFLOW.')
      IF((RANK.EQ.1).AND.(DIMSR(1).EQ.1)) THEN
        CALL hdf5_read_data(IPMPO,"/info/MPO_CREATION_INFO",TEXT132)
        IF(IMPX.GT.0) WRITE(IOUT,'(1X,A)') TEXT132
      ELSE IF(RANK.EQ.1) THEN
        CALL hdf5_read_data(IPMPO,"/info/MPO_CREATION_INFO",TEXT132V1)
        IF(IMPX.GT.0) THEN
          DO I=1,DIMSR(1)
            WRITE(IOUT,'(1X,A)') TEXT132V1(I)
          ENDDO
        ENDIF
        DEALLOCATE(TEXT132V1)
      ENDIF
      IF(NPAR.GT.0) THEN
        CALL hdf5_read_data(IPMPO,"/parameters/info/PARAMTYPE",PARTYP)
        CALL hdf5_read_data(IPMPO,"/parameters/info/PARAMNAME",PARKEY)
        CALL hdf5_read_data(IPMPO,"/parameters/info/PARAMFORM",PARFMT)
        IF(IMPX.GT.1) THEN
          WRITE(IOUT,*) 'NPAR=',NPAR,SIZE(PARKEY,1)
          DO I=1,NPAR
            WRITE(IOUT,*)'PARKEY(',I,')=',PARKEY(I),' PARFMT=',PARFMT(I)
          ENDDO
        ENDIF
      ENDIF
      TERP(:NCAL,:NMIX)=0.0
      MIXC(:NMIX)=0
*----
*  SCAN THE MPO FILE INFORMATION TO RECOVER THE MUPLET DATABASE
*----
      IF(IMPX.GT.5) THEN
        WRITE(IOUT,'(24H MCRDRV: MUPLET DATABASE/12H CALCULATION,5X,
     1  10HMUPLET....)')
      ENDIF
      DO ICAL=1,NCAL
        WRITE(RECNAM,'(8H/output/,A,9H/statept_,I0)') TRIM(HEDIT),ICAL-1
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"/PARAMVALUEORD",VINTE)
        IF(SIZE(VINTE).NE.NPAR) THEN
          WRITE(HSMG,'(43HMCRDRV: INCONSISTENT PARAMVALUEORD LENGTH (,
     1    I5,3H VS,I5,2H).)') SIZE(VINTE),NPAR
          CALL XABORT(HSMG)
        ENDIF
        DO IPAR=1,NPAR
          MUBASE(IPAR,ICAL)=VINTE(IPAR)+1
        ENDDO
        IF(IMPX.GT.5) THEN
          WRITE(IOUT,'(I8,6X,20I4/(14X,20I4))') ICAL,
     1    MUBASE(:,ICAL)
        ENDIF
        DEALLOCATE(VINTE)
      ENDDO
*----
*  READ (INTERP_DATA) AND SET VALI, VALR AND VALH PARAMETERS
*  CORRESPONDING TO THE INTERPOLATION POINT. FILL MUPLET FOR
*  PARAMETERS SET WITHOUT INTERPOLATION.
*----
      IBM=0
      MAXNIS=0
      NISO(:NMIX)=0
      LISO(:NMIX)=.TRUE.
      LDELTA(:NMIX)=.FALSE.
      ITODO(:NMIX,:NBISO)=0
   10 CALL REDGET(INDIC,NITMA,FLOTT,TEXT72,DFLOTT)
      IF(INDIC.NE.3) CALL XABORT('MCRDRV: CHARACTER DATA EXPECTED.')
   20 IF(TEXT72.EQ.'MIX') THEN
         MUPLET(:NPAR)=0
         MUTYPE(:NPAR)=0
         VALI(:NPAR)=0
         VALR(:NPAR,1)=0.0
         VALR(:NPAR,2)=0.0
         DO 30 I=1,NPAR
         VALH(I)=' '
   30    CONTINUE
         LCUB2(:NPAR)=LCUBIC
         CALL REDGET(INDIC,IBM,FLOTT,TEXT72,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('MCRDRV: INTEGER DATA EXPECTED.')
         IF(IBM.GT.NMIX) THEN
            WRITE(HSMG,'(27HMCRDRV: NMIX OVERFLOW (IBM=,I8,6H NMIX=,I8,
     1      2H).)') IBM,NMIX
            CALL XABORT(HSMG)
         ENDIF
         IBMOLD=1
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT72,DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('MCRDRV: CHARACTER DATA EXPECTED.')
         IF(TEXT72.EQ.'FROM') THEN
            CALL REDGET(INDIC,IBMOLD,FLOTT,TEXT72,DFLOTT)
            IF(INDIC.NE.1) CALL XABORT('MCRDRV: INTEGER DATA EXPECTED.')
            IF(IBMOLD.GT.NMIL) CALL XABORT('MCRDRV: MPO MIX OVERFLOW'
     1      //'(1).')
            MIXC(IBM)=IBMOLD
            GO TO 10
         ELSE IF(TEXT72.EQ.'USE') THEN
            IF(IBM.GT.NMIL) CALL XABORT('MCRDRV: MPO MIX OVERFLOW(2).')
            MIXC(IBM)=IBM
            GO TO 10
         ENDIF
         MIXC(IBM)=IBMOLD
         GO TO 20
      ELSE IF(TEXT72.EQ.'MICRO') THEN
         IF(IBM.EQ.0) CALL XABORT('MCRDRV: MIX NOT SET (1).')
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT72,DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('MCRDRV: CHARACTER DATA EXPECTED.')
         IF(TEXT72.EQ.'ALL') THEN
            LISO(IBM)=.TRUE.
         ELSE IF(TEXT72.EQ.'ONLY') THEN
            LISO(IBM)=.FALSE.
         ENDIF
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT72,DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('MCRDRV: CHARACTER DATA EXPECTED.')
   40    IF(TEXT72.EQ.'ENDMIX') THEN
           GO TO 20
         ELSE IF(TEXT72.EQ.'NOEV') THEN
           IF(NISO(IBM).EQ.0) CALL XABORT('MCRDRV: MISPLACED NOEV.')
           ITODO(IBM,NISO(IBM))=1
         ELSE
           NISO(IBM)=NISO(IBM)+1
           IF(NISO(IBM).GT.NBISO) CALL XABORT('MCRDRV: NBISO OVERFLOW.')
           MAXNIS=MAX(MAXNIS,NISO(IBM))
           HISO(IBM,NISO(IBM))=TEXT72(:8)
           CALL REDGET(INDIC,NITMA,FLOTT,TEXT72,DFLOTT)
           IF(INDIC.EQ.2) THEN
             CONC(IBM,NISO(IBM))=FLOTT
           ELSE IF((INDIC.EQ.3).AND.(TEXT72.EQ.'*')) THEN
             CONC(IBM,NISO(IBM))=-99.99
           ELSE
             CALL XABORT('MCRDRV: INVALID HISO DATA.')
           ENDIF
         ENDIF
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT72,DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('MCRDRV: CHARACTER DATA EXPECTED.')
         GO TO 40
      ELSE IF((TEXT72.EQ.'SET').OR.(TEXT72.EQ.'DELTA')) THEN
         IF(IBM.EQ.0) CALL XABORT('MCRDRV: MIX NOT SET (2).')
         ITYPE=0
         IF(TEXT72.EQ.'SET') THEN
            ITYPE=1
         ELSE IF(TEXT72.EQ.'DELTA') THEN
            ITYPE=2
            LDELTA(IBM)=.TRUE.
         ENDIF
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT72,DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('MCRDRV: CHARACTER DATA EXPECTED.')
         IF((TEXT72.EQ.'LINEAR').OR.(TEXT72.EQ.'CUBIC')) THEN
            HCUBIC=TEXT72(:12)
            CALL REDGET(INDIC,NITMA,FLOTT,TEXT72,DFLOTT)
         ELSE
            HCUBIC=' '
         ENDIF
         IF(INDIC.NE.3) CALL XABORT('MCRDRV: CHARACTER DATA EXPECTED.')
         IPAR=0
         DO 50 I=1,NPAR
         IF(TEXT72.EQ.PARKEY(I)) THEN
            IPAR=I
            GO TO 60
         ENDIF
   50    CONTINUE
         CALL XABORT('MCRDRV: PARAMETER '//TEXT72//' NOT FOUND.')
   60    IF(HCUBIC.EQ.'LINEAR') THEN
            LCUB2(IPAR)=.FALSE.
         ELSE IF(HCUBIC.EQ.'CUBIC') THEN
            LCUB2(IPAR)=.TRUE.
         ENDIF
         CALL hdf5_read_data(IPMPO,"/parameters/info/NVALUE",NVALUE)
         IF(NVALUE(IPAR).GT.MAXVAL) CALL XABORT('MCRDRV: MAXVAL OVERFL'
     1   //'OW.')
         WRITE(RECNAM,'(25H/parameters/values/PARAM_,I0)') IPAR-1
         CALL hdf5_info(IPMPO,RECNAM,RANK,TYPE,NBYTE,DIMSR)
         IF(TYPE.EQ.99) THEN
            WRITE(HSMG,'(25HMCRDRV: GLOBAL PARAMETER ,A,9H NOT SET.)')
     1      TRIM(PARKEY(IPAR))
            CALL XABORT(HSMG)
         ENDIF
         IF(PARFMT(IPAR).EQ.'INTEGER') THEN
            IF(ITYPE.NE.1) CALL XABORT('MCRDRV: SET MANDATORY WITH INT'
     1      //'EGER PARAMETERS.')
            CALL REDGET(INDIC,VALI(IPAR),FLOTT,TEXT72,DFLOTT)
            IF(INDIC.NE.1) CALL XABORT('MCRDRV: INTEGER DATA EXPECTED.')
            CALL hdf5_read_data(IPMPO,RECNAM,VINTE)
            DO 70 J=1,NVALUE(IPAR)
            IF(VALI(IPAR).EQ.VINTE(J)) THEN
               MUPLET(IPAR)=J
               MUTYPE(IPAR)=ITYPE
               DEALLOCATE(NVALUE,VINTE)
               GO TO 10
            ENDIF
   70       CONTINUE
            WRITE(HSMG,'(26HMCRDRV: INTEGER PARAMETER ,A,11H WITH VALUE,
     1      I5,32H NOT FOUND IN MPO FILE DATABASE.)') 
     2      TRIM(PARKEY(IPAR)),VALI(IPAR)
            CALL XABORT(HSMG)
         ELSE IF(PARFMT(IPAR).EQ.'FLOAT') THEN
            CALL REDGET(INDIC,NITMA,VALR(IPAR,1),TEXT72,DFLOTT)
            IF(INDIC.NE.2) CALL XABORT('MCRDRV: REAL DATA EXPECTED.')
            VALR(IPAR,2)=VALR(IPAR,1)
            CALL REDGET(INDIC,NITMA,FLOTT,TEXT72,DFLOTT)
            IF(INDIC.EQ.2) THEN
               VALR(IPAR,2)=FLOTT
               CALL REDGET(INDIC,NITMA,FLOTT,TEXT72,DFLOTT)
            ENDIF
            CALL hdf5_read_data(IPMPO,RECNAM,VREAL)
            IF(VALR(IPAR,1).EQ.VALR(IPAR,2)) THEN
               DO 80 J=1,NVALUE(IPAR)
               IF(ABS(VALR(IPAR,1)-VREAL(J)).LE.REPS*ABS(VREAL(J))) THEN
                  MUPLET(IPAR)=J
                  IF(ITYPE.NE.1) MUPLET(IPAR)=-1
                  MUTYPE(IPAR)=ITYPE
                  DEALLOCATE(NVALUE,VREAL)
                  GO TO 20
               ENDIF
   80          CONTINUE
            ENDIF
            IF(VALR(IPAR,1).LT.VREAL(1)) THEN
               WRITE(HSMG,'(23HMCRDRV: REAL PARAMETER ,A,10H WITH VALU,
     1         1HE,1P,E12.4,25H IS OUTSIDE THE DOMAIN (<,E12.4,1H))')
     2         TRIM(PARKEY(IPAR)),VALR(IPAR,1),VREAL(1)
               CALL XABORT(HSMG)
            ELSE IF(VALR(IPAR,2).GT.VREAL(NVALUE(IPAR))) THEN
               WRITE(HSMG,'(23HMCRDRV: REAL PARAMETER ,A,10H WITH VALU,
     1         1HE,1P,E12.4,25H IS OUTSIDE THE DOMAIN (>,E12.4,1H))')
     2         TRIM(PARKEY(IPAR)),VALR(IPAR,1),VREAL(NVALUE(IPAR))
               CALL XABORT(HSMG)
            ELSE IF(VALR(IPAR,1).GT.VALR(IPAR,2)) THEN
               WRITE(HSMG,'(23HMCRDRV: REAL PARAMETER ,A,9H IS DEFIN,
     1         7HED WITH,1P,E12.4,2H >,E12.4,1H.)') TRIM(PARKEY(IPAR)),
     2         VALR(IPAR,1),VALR(IPAR,2)
               CALL XABORT(HSMG)
            ENDIF
            MUPLET(IPAR)=-1
            MUTYPE(IPAR)=ITYPE
            DEALLOCATE(NVALUE,VREAL)
            GO TO 20
         ELSE IF(PARFMT(IPAR).EQ.'STRING') THEN
            IF(ITYPE.NE.1) CALL XABORT('MCRDRV: SET MANDATORY WITH STR'
     1      //'ING PARAMETERS.')
            CALL REDGET(INDIC,NITMA,FLOTT,VALH(IPAR),DFLOTT)
            IF(INDIC.NE.3) CALL XABORT('MCRDRV: STRING DATA EXPECTED.')
            CALL hdf5_read_data(IPMPO,RECNAM,VCHAR)
            DO 90 J=1,NVALUE(IPAR)
            IF(VALH(IPAR).EQ.VCHAR(J)) THEN
               MUPLET(IPAR)=J
               MUTYPE(IPAR)=ITYPE
               DEALLOCATE(NVALUE,VCHAR)
               GO TO 10
            ENDIF
   90       CONTINUE
            WRITE(HSMG,'(25HMCRDRV: STRING PARAMETER ,A,10H WITH VALU,
     1      2HE ,A12,32H NOT FOUND IN MPO FILE DATABASE.)')
     2      TRIM(PARKEY(IPAR)), VALH(IPAR)
            CALL XABORT(HSMG)
         ELSE
            CALL XABORT('MCRDRV: INVALID FORMAT='//PARFMT(IPAR))
         ENDIF
      ELSE IF(TEXT72.EQ.'ENDMIX') THEN
*----
*  COMPUTE THE TERP FACTORS USING TABLE-OF-CONTENT INFORMATION.
*----
         IF(IMPX.GT.0) THEN
           DO IPAR=1,NPAR
             IF(PARFMT(IPAR).EQ.'FLOAT')THEN
               IF(LCUB2(IPAR)) THEN
                 WRITE(IOUT,'(26H MCRDRV: GLOBAL PARAMETER:,A,7H ->CUBI,
     1           16HC INTERPOLATION.)') TRIM(PARKEY(IPAR))
               ELSE
                 WRITE(IOUT,'(26H MCRDRV: GLOBAL PARAMETER:,A,7H ->LINE,
     1           17HAR INTERPOLATION.)') TRIM(PARKEY(IPAR))
               ENDIF
             ENDIF
           ENDDO
         ENDIF
         IF(IBMOLD.GT.NMIL)CALL XABORT('MCRDRV: MPO MIX OVERFLOW(3).')
         IF(IBM.GT.NMIX)CALL XABORT('MCRDRV: MIX OVERFLOW (MICROLIB).')
         IF(NCAL.EQ.1) THEN
           TERP(1,IBM)=1.0
         ELSE
           CALL MCRTRP(IPMPO,LCUB2,IMPX,NPAR,NCAL,MUPLET,MUTYPE,PARTYP,
     1     VALR,0.0,MUBASE,TERP(1,IBM))
         ENDIF
         IBM=0
      ELSE IF((TEXT72.EQ.'MPO').OR.(TEXT72.EQ.'TABLE').OR.
     1   (TEXT72.EQ.'CHAIN').OR.(TEXT72.EQ.';')) THEN
*----
*  CHECK TERP FACTORS AND RETURN
*----
         IF(TEXT72.EQ.';') ITER=0
         IF(TEXT72.EQ.'MPO') ITER=1
         IF(TEXT72.EQ.'TABLE') ITER=2
         IF(TEXT72.EQ.'CHAIN') ITER=3
         DO 150 IBM=1,NMIX
         IBMOLD=MIXC(IBM)
         IF(IBMOLD.EQ.0) GO TO 150
         IF(NISO(IBM).GT.MAXNIS) CALL XABORT('MCRDRV: MAXNIS OVERFLOW.')
         IF(LDELTA(IBM)) THEN
            SUM=0.0
         ELSE
            SUM=1.0
         ENDIF
         DO 140 ICAL=1,NCAL
         SUM=SUM-TERP(ICAL,IBM)
  140    CONTINUE
         IF(ABS(SUM).GT.1.0E-4) THEN
            WRITE(HSMG,'(43HMCRDRV: INVALID INTERPOLATION FACTORS IN MI,
     1      5HXTURE,I4,1H.)') IBM
            CALL XABORT(HSMG)
         ENDIF
  150    CONTINUE
         GO TO 160
      ELSE
         CALL XABORT('MCRDRV: '//TEXT72//' IS AN INVALID KEYWORD.')
      ENDIF
      GO TO 10
*----
*  PRINT INTERPOLATION (TERP) FACTORS
*----
  160 IF(IMPX.GT.2) THEN
        WRITE(IOUT,'(/30H MCRDRV: INTERPOLATION FACTORS)')
        DO ICAL=1,NCAL
          DO IBM=1,NMIX
            IF(TERP(ICAL,IBM).NE.0.0) THEN
              WRITE(IOUT,170) ICAL,(TERP(ICAL,JBM),JBM=1,NMIX)
              EXIT
            ENDIF
          ENDDO
        ENDDO
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      IF(NPAR.GT.0) DEALLOCATE(PARFMT,PARKEY,PARTYP)
      DEALLOCATE(MUBASE,LDELTA,MUTYPE,MUPLET)
      RETURN
  170 FORMAT(6H CALC=,I8,6H TERP=,1P,8E13.5/(20X,8E13.5))
      END
