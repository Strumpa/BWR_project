*DECK CVRDRV
      SUBROUTINE CVRDRV(IPMAP,NCH,NB,NFUEL,NPARM,NX,NY,NZ,NVOID,IVOID,
     1 IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read the input data required for the voiding simulations.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* D. Sekki
*
*Parameters: input
* IPMAP  pointer to the perturbed fuel-map.
* NCH    number of reactor channels.
* NB     number of fuel bundles per channel.
* NFUEL  number of fuel types.
* NPARM  total number of recorded parameters.
* NX     number of elements along x-axis in fuel map.
* NY     number of elements along y-axis in fuel map.
* NZ     number of elements along z-axis in fuel map.
* IMPX   printing index (=0 for no print).
*
*Parameters: output
* IVOID  index associated with the core-voiding pattern.
* NVOID  total number of voided channels.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAP
      INTEGER NCH,NB,NFUEL,NPARM,NX,NY,NZ,NVOID,IVOID,IMPX
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOUT=6)
      CHARACTER TEXT*12,PNAME*12,CHANX*4,CHANY*4
      INTEGER INAME(3)
      DOUBLE PRECISION DFLOT
      REAL PVALUE(NCH,NB)
      LOGICAL LCOOL
      TYPE(C_PTR) JPMAP,KPMAP
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NXX,NYY
*----
*  FUEL-TYPE INDICES
*----
      JPMAP=LCMGID(IPMAP,'FUEL')
      DO IFUEL=1,NFUEL
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.3)CALL XABORT('@CVRDRV: CHARACTER DATA EXPECTED(1).')
        IF(TEXT.NE.'MIX-FUEL')CALL XABORT('@CVRDRV: KEYWORD MIX-FUEL'
     1   //' EXPECTED.')
*       UNPERTURBED-CELL FUEL MIXTURE NUMBER
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@CVRDRV: INTEGER DATA EXPECTED(1).')
        KPMAP=LCMGIL(JPMAP,IFUEL)
        CALL LCMGET(KPMAP,'MIX',IMIX)
        IF(IMIX.NE.NITMA)THEN
          WRITE(IOUT,*)'@CVRDRV: RECORDED FUEL MIXTURE NUMBER ',IMIX
          WRITE(IOUT,*)'@CVRDRV: READ FROM INPUT THE MIXTURE ',NITMA
          CALL XABORT('@CVRDRV: WRONG INPUT ORDER OF FUEL MIXTURES.')
        ENDIF
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.3)CALL XABORT('@CVRDRV: CHARACTER DATA EXPECTED(2).')
        IF(TEXT.NE.'MIX-VOID')CALL XABORT('@CVRDRV: KEYWORD MIX-VOI'
     1   //'D EXPECTED.')
*       PERTURBED-CELL FUEL MIXTURE NUMBER
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@CVRDRV: INTEGER DATA EXPECTED(2).')
        IF(NITMA.LE.0)CALL XABORT('@CVRDRV: MIX-VOID NUMBER MUST BE'
     1   //' POSITIVE AND GREATER THAN ZERO.')
        CALL LCMPUT(KPMAP,'MIX-VOID',1,1,NITMA)
      ENDDO
*----
*  COOLANT DENSITIES
*----
      LCOOL=.FALSE.
      CALL XDRSET(PVALUE,NCH*NB,0.)
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.3)CALL XABORT('@CVRDRV: CHARACTER DATA EXPECTED(3).')
      IF(TEXT.NE.'DENS-COOL')GOTO 20
      IF(NPARM.EQ.0)CALL XABORT('@CVRDRV: NO DEFINED PARAMETERS IN T'
     1 //'HE FUEL-MAP NPARM=0')
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.3)CALL XABORT('@CVRDRV: CHARACTER DATA FOR PARAMETE'
     1 //'R PNAME EXPECTED.')
      JPMAP=LCMGID(IPMAP,'PARAM')
      DO IPAR=1,NPARM
        KPMAP=LCMGIL(JPMAP,IPAR)
        CALL LCMGET(KPMAP,'P-NAME',INAME)
        WRITE(PNAME,'(3A4)') (INAME(I),I=1,3)
        IF(PNAME.EQ.TEXT)THEN
          CALL LCMGET(KPMAP,'P-VALUE',PVALUE)
          GOTO 10
        ENDIF
      ENDDO
      CALL XABORT('@CVRDRV: UNABLE TO FIND PARAMETER WITH PNAME '//TEXT)
*
   10 CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(TEXT.NE.'SET')CALL XABORT('@CVRDRV: KEYWORD SET EXPECTED.')
      CALL REDGET(ITYP,NITMA,VCOOL,TEXT,DFLOT)
      IF(ITYP.NE.2)CALL XABORT('@CVRDRV: REAL DATA FOR THE COOLANT DEN'
     1 //'SITY EXPECTED.')
      IF(VCOOL.LT.0.)CALL XABORT('@CVRDRV: INVALID VALUE FOR THE COOLA'
     1 //'NT DENSITY <0.')
      LCOOL=.TRUE.
*----
*  CORE-VOIDING PATTERN
*----
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
   20 IF(TEXT.NE.'VOID-PATTERN')CALL XABORT('@CVRDRV: KEYWORD VOID-'
     1 //'PATTERN EXPECTED.')
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.3)CALL XABORT('@CVRDRV: CHARACTER DATA EXPECTED(4).')
      IF(TEXT.EQ.'FULL')THEN
        IVOID=1
        NVOID=NCH
      ELSEIF(TEXT.EQ.'HALF')THEN
        IVOID=2
        NVOID=NCH/2
      ELSEIF(TEXT.EQ.'QUARTER')THEN
        IVOID=3
        NVOID=NCH/4
      ELSEIF(TEXT.EQ.'CHECKER')THEN
        IVOID=4
        NVOID=NCH
      ELSEIF(TEXT.EQ.'CHECKER-1/2')THEN
        IVOID=5
        NVOID=NCH/2
      ELSEIF(TEXT.EQ.'CHECKER-1/4')THEN
        IVOID=6
        NVOID=NCH/4
      ELSEIF(TEXT.EQ.'CHAN-VOID')THEN
*----
*  USER-DEFINED PATTERN
*----
        IVOID=7
*       TOTAL NUMBER OF VOIDED CHANNELS
        CALL REDGET(ITYP,NVOID,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@CVRDRV: INTEGER TOTAL NUMBER OF V'
     1   //'OIDED CHANNELS EXPECTED.')
        IF((NVOID.LT.1).OR.(NVOID.GT.NCH))CALL XABORT('@CVRDRV: TH'
     1   //'E NUMBER OF VOIDED CHANNELS MUST BE > 0 AND < NCH')
        ALLOCATE(NXX(NVOID),NYY(NVOID))
        DO I=1,NVOID
*         VOIDED-CHANNEL YNAME
          CALL REDGET(ITYP,NITMA,FLOT,CHANY,DFLOT)
          IF(ITYP.NE.3)CALL XABORT('@CVRDRV: CHANNEL YNAME EXPECTED.')
          READ(CHANY,'(A4)') NYY(I)
*         VOIDED-CHANNEL XNAME
          CALL REDGET(ITYP,NITMA,FLOT,CHANX,DFLOT)
          IF(ITYP.NE.3)CALL XABORT('@CVRDRV: CHANNEL XNAME EXPECTED.')
          READ(CHANX,'(A4)') NXX(I)
        ENDDO
        CALL CVRUSR(IPMAP,NCH,NB,NFUEL,NX,NY,NZ,NVOID,NXX,NYY,NPARM,
     1  PNAME,PVALUE,VCOOL,LCOOL,IMPX)
        DEALLOCATE(NXX,NYY)
      ELSE
        CALL XABORT('@CVRDRV: WRONG KEYWORD '//TEXT)
      ENDIF
*
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(TEXT.NE.';')CALL XABORT('@CVRDRV: FINAL ; EXPECTED.')
*     SPECIFIED CORE-VOIDING PATTERN
      IF(IVOID.LT.7) CALL CVRCOR(IPMAP,NCH,NB,NFUEL,NX,NY,NZ,IVOID,
     1 NVOID,NPARM,PNAME,PVALUE,VCOOL,LCOOL,IMPX)
      RETURN
      END