*DECK S2MGET
      SUBROUTINE S2MGET(IFIN,HNAME,IDX,LRMS,NGRP,XS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To read a record in the Matlab-formatted SERPENT output file.
*
*Copyright:
* Copyright (C) 2014 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert and C. Bienvenue
*
*Parameters: input
* IFIN    SERPENT unit number.
* HNAME   record name.
* IDX     instance index on the record (generally equal to the burnup
*         step).
* LRMS    standard deviation flag (.TRUE. if the standard deviation is
*         present).
* NGRP    number of values expected in the record.
*
*Parameters: output
* XS      record recovered from SERPENT file.
*
*-----------------------------------------------------------------------
*
      USE ISO_FORTRAN_ENV, ONLY : IOSTAT_EOR
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IFIN,IDX,NGRP
      LOGICAL LRMS
      CHARACTER*(*) HNAME
      REAL XS(NGRP)
*----
*  LOCAL VARIABLES
*----
      CHARACTER(LEN=:), ALLOCATABLE :: HLINE
      CHARACTER PREFIX*26,HSMG*131
      INTEGER N,IOS,IB1,IB2,IB3,IB4,ICOL,NVAL,NEXP,IGR
      REAL, ALLOCATABLE, DIMENSION(:) :: VALUES
*
      IF(LEN(HNAME).GT.26) CALL XABORT('S2MGET: PREFIX OVERFLOW.')
      PREFIX=HNAME
*----
*  LOCATE THE IDX-TH INSTANCE OF THE RECORD
*----
      REWIND(IFIN)
      N=0
      DO
        CALL S2MRDL(IOS)
        IF(IOS.NE.0) GO TO 10
        IF(INDEX(HLINE,PREFIX).GT.0) THEN
          N=N+1
          IF(N.EQ.IDX) GO TO 20
        ENDIF
      ENDDO
   10 WRITE(HSMG,'(''S2MGET: UNABLE TO FIND '',I5,'' INSTANCE(S) OF ''
     > ,A,'' RECORD ('',I5,'' FOUND).'')') IDX,TRIM(PREFIX),N
      CALL XABORT(HSMG)
*----
*  RECOVER THE DECLARED NUMBER OF ITEMS FROM THE '(idx, [1:   N])' HEADER
*----
   20 IB1=INDEX(HLINE,'[')
      IB2=INDEX(HLINE,']')
      IF((IB1.EQ.0).OR.(IB2.LE.IB1)) THEN
        WRITE(HSMG,'(''S2MGET: MALFORMED HEADER ON THE '',A,
     >  '' RECORD.'')') TRIM(PREFIX)
        CALL XABORT(HSMG)
      ENDIF
      ICOL=INDEX(HLINE(IB1:IB2),':')
      IF(ICOL.EQ.0) THEN
        WRITE(HSMG,'(''S2MGET: MALFORMED HEADER ON THE '',A,
     >  '' RECORD.'')') TRIM(PREFIX)
        CALL XABORT(HSMG)
      ENDIF
      READ(HLINE(IB1+ICOL:IB2-1),*,IOSTAT=IOS) NVAL
      IF(IOS.NE.0) THEN
        WRITE(HSMG,'(''S2MGET: UNABLE TO RECOVER THE NUMBER OF ITEMS''
     >  ,'' ON THE '',A,'' RECORD.'')') TRIM(PREFIX)
        CALL XABORT(HSMG)
      ENDIF
      IF(LRMS) THEN
        NEXP=2*NGRP
      ELSE
        NEXP=NGRP
      ENDIF
      IF(NVAL.NE.NEXP) THEN
        WRITE(HSMG,'(''S2MGET: '',A,'' RECORD: '',I8,'' ITEM(S) FOUND''
     >  ,'' WHILE '',I8,'' WERE EXPECTED.'')') TRIM(PREFIX),NVAL,NEXP
        CALL XABORT(HSMG)
      ENDIF
*----
*  RECOVER THE VALUES FROM THE DATA BRACKET
*----
      IB3=IB2+INDEX(HLINE(IB2+1:),'[')
      IB4=IB2+INDEX(HLINE(IB2+1:),']')
      IF((IB3.LE.IB2).OR.(IB4.LE.IB3)) THEN
        WRITE(HSMG,'(''S2MGET: NO DATA BRACKET ON THE '',A,
     >  '' RECORD.'')') TRIM(PREFIX)
        CALL XABORT(HSMG)
      ENDIF
      ALLOCATE(VALUES(NVAL))
      READ(HLINE(IB3+1:IB4-1),*,IOSTAT=IOS) VALUES
      IF(IOS.NE.0) THEN
        WRITE(HSMG,'(''S2MGET: UNABLE TO DECODE THE '',A,
     >  '' RECORD.'')') TRIM(PREFIX)
        CALL XABORT(HSMG)
      ENDIF
      IF(LRMS) THEN
        DO IGR=1,NGRP
          XS(IGR)=VALUES(2*IGR-1)
        ENDDO
      ELSE
        DO IGR=1,NGRP
          XS(IGR)=VALUES(IGR)
        ENDDO
      ENDIF
      DEALLOCATE(VALUES)
      IF(ALLOCATED(HLINE)) DEALLOCATE(HLINE)
      RETURN
*
      CONTAINS
*
      SUBROUTINE S2MRDL(JOS)
*----
*  Recover a complete record of arbitrary length from unit IFIN into the
*  host-associated allocatable string HLINE. JOS is set to a non-zero
*  value at end of file.
*----
      INTEGER JOS,NCHAR,NOLD
      CHARACTER CHUNK*1024
      CHARACTER(LEN=:), ALLOCATABLE :: HTMP
*
      IF(ALLOCATED(HLINE)) DEALLOCATE(HLINE)
      ALLOCATE(CHARACTER(LEN=0) :: HLINE)
      DO
        NCHAR=0
        READ(IFIN,'(A)',ADVANCE='NO',SIZE=NCHAR,IOSTAT=JOS) CHUNK
        IF(NCHAR.GT.0) THEN
          NOLD=LEN(HLINE)
          ALLOCATE(CHARACTER(LEN=NOLD+NCHAR) :: HTMP)
          IF(NOLD.GT.0) HTMP(:NOLD)=HLINE(:NOLD)
          HTMP(NOLD+1:)=CHUNK(:NCHAR)
          CALL MOVE_ALLOC(HTMP,HLINE)
        ENDIF
        IF(JOS.NE.0) EXIT
      ENDDO
      IF(JOS.EQ.IOSTAT_EOR) JOS=0
      END SUBROUTINE S2MRDL
*
      END
