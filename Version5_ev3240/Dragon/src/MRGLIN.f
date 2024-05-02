*DECK MRGLIN
      SUBROUTINE MRGLIN(IPRINT,IFTRKO,NSOUTO,NVOUTO,IFTRKN,
     >                  IMERGE,NDIM,IFMT,MXSUB,MXSEG)
*
*----------
*
*Purpose:
* Merge volume surface information on track file.
*
*Copyright:
* Copyright (C) 1997 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* G. Marleau
*
*Parameters: input
* IPRINT  print level.
* IFTRKO  old tracking file.
* NSOUTO  old number of surfaces.
* NVOUTO  old number of regions.
* IFTRKN  new tracking file.
* IFMT    file format
* IMERGE  merged position.
* IFMT    track format: =0 short; =1 long.
* MXSUB  maximum number of subtracks in a track.
* MXSEG   maximum number of segments.
*
*----------
*
      IMPLICIT         NONE
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='MRGLIN')
*----
*  ROUTINE PARAMETERS
*----
      INTEGER          IPRINT,IFTRKO,NSOUTO,NVOUTO,IFTRKN,
     >                 NDIM,IFMT,MXSUB,MXSEG
      INTEGER          IMERGE(-NSOUTO:NVOUTO)
*----
*  LOCAL VARIABLES
*----
      INTEGER          ITRAK,NLINEO,NLINEN,ILINE,
     >                 ISEG,IVSO,NSUB,IADD(4),IRA,ISU
      DOUBLE PRECISION WEIGHT
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NRSEG,IANGL
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: PATH
      DOUBLE PRECISION , ALLOCATABLE, DIMENSION(:,:) :: DADD
*----
*  LOOP OVER TRACKS
*----
      ALLOCATE(NRSEG(MXSEG),PATH(MXSEG))
      ALLOCATE(IANGL(MXSUB),DADD(NDIM,MXSUB))
      ITRAK=0
 1000 CONTINUE
        IF(IFMT.EQ.1) THEN
          READ (IFTRKO,END=1010) NSUB,NLINEO,WEIGHT,
     >                           (IANGL(IRA),IRA=1,NSUB),
     >                           (NRSEG(ILINE),ILINE=1,NLINEO),
     >                           (PATH(ILINE),ILINE=1,NLINEO),
     >                           (IADD(IRA),IRA=1,4),
     >                          ((DADD(IRA,ISU),IRA=1,NDIM),ISU=1,NSUB)
       ELSE
          READ (IFTRKO,END=1010) NSUB,NLINEO,WEIGHT,
     >                           (IANGL(IRA),IRA=1,NSUB),
     >                           (NRSEG(ILINE),ILINE=1,NLINEO),
     >                           (PATH(ILINE),ILINE=1,NLINEO)
       ENDIF
*----
*  SCAN NRSEG AND RESET TO NEW VOLUME AND SURFACE NUMBER
*----
        ITRAK=ITRAK+1
        IF(IPRINT.GE.1000) THEN
          WRITE(IOUT,6000) ITRAK,NLINEO,WEIGHT,IANGL
          WRITE(IOUT,6010)
     >      (NRSEG(ILINE),PATH(ILINE),ILINE=1,NLINEO)
        ENDIF
        DO 100 ILINE=1,NLINEO
          DO 110 IVSO=-NSOUTO,NVOUTO
            IF(NRSEG(ILINE) .EQ. IVSO ) THEN
              NRSEG(ILINE) = IMERGE(IVSO)
              GO TO 115
            ENDIF
 110      CONTINUE
 115      CONTINUE
 100    CONTINUE
*----
*  COMPRESS REGION OF SUCCESSIVE IDENTICAL REGION
*  EXCEPT FOR SURFACES
*----
        NLINEN=1
        ISEG=NRSEG(NLINEN)
        DO 120 ILINE=2,NLINEO
          IF(NRSEG(ILINE) .EQ. ISEG .AND.
     >       ISEG .GT. 0                  ) THEN
            PATH(NLINEN)=PATH(NLINEN)+PATH(ILINE)
          ELSE
            NLINEN=NLINEN+1
            NRSEG(NLINEN)=NRSEG(ILINE)
            PATH(NLINEN)=PATH(ILINE)
            ISEG=NRSEG(NLINEN)
          ENDIF
 120    CONTINUE
        IF(IFMT.EQ.1) THEN
          WRITE(IFTRKN) NSUB,NLINEN,WEIGHT,
     >                  (IANGL(IRA),IRA=1,NSUB),
     >                  (NRSEG(ILINE),ILINE=1,NLINEN),
     >                  (PATH(ILINE),ILINE=1,NLINEN),
     >                  (IADD(IRA),IRA=1,4),
     >                 ((DADD(IRA,ISU),IRA=1,NDIM),ISU=1,NSUB)
        ELSE
          WRITE(IFTRKN) NSUB,NLINEN,WEIGHT,
     >                  (IANGL(IRA),IRA=1,NSUB),
     >                  (NRSEG(ILINE),ILINE=1,NLINEN),
     >                  (PATH(ILINE),ILINE=1,NLINEN)
        ENDIF
        IF(IPRINT.GE.1000) THEN
          WRITE(IOUT,6001) ITRAK,NLINEN,WEIGHT,IANGL
          WRITE(IOUT,6010)
     >      (NRSEG(ILINE),PATH(ILINE),ILINE=1,NLINEN)
        ENDIF
      GO TO 1000
 1010 CONTINUE
      DEALLOCATE(DADD,IANGL)
      DEALLOCATE(PATH,NRSEG)
*----
*  FORMAT
*----
 6000 FORMAT(' INITIAL LINE = ',I10/
     >       ' PARAMETERS = ',I10,F15.7,10I10)
 6001 FORMAT(' FINAL LINE = ',I10/
     >       ' PARAMETERS = ',I10,F15.7,10I10)
 6010 FORMAT(1P,5(I10,E15.7))
      RETURN
      END
