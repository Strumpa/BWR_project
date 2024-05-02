*DECK SNADPT
      SUBROUTINE SNADPT(IELEM,NM,NMX,M,MX,TB,W,ISFIX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the weighting parameters for SN adaptive flux calculation,
* based on flux variation in the cell.
*
*Copyright:
* Copyright (C) 2020 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): C. Bienvenue
*
*Parameters: input
* IELEM   measure of order of the closure relation.
* NM      number of flux moments.
* NMX     number of incoming flux moments.
* M       moments of the flux.
* MX      incoming flux moments.
* TB      ratio (space=1, energy=ratio of boundary stopping powers).
*
*Parameters: input and output
* W       weighting parameters of the closure relation.
* ISFIX   flag indicating if moments sould be recalculated.
*
*-----------------------------------------------------------------------
*
*----
* SUBROUTINE ARGUMENTS
*----
      INTEGER IELEM,NM,NMX
      REAL W(IELEM+1),TB
      DOUBLE PRECISION M(NM),MX(NMX)
      LOGICAL ISFIX
*----
* LOCAL VARIABLES
*----
      REAL P,U1,F1
      PARAMETER(EPS=1.0E-8,B=2.0)
*----
* CONSTANT ORDER ADAPTIVE CALCULATIONS 
*----
      IF(IELEM.EQ.1) THEN

      ! EXTRACT P
      P=-W(1)/TB

      ! CASE P = 1
      IF(ABS(P-1).LE.EPS) THEN
      IF(M(1).NE.0) THEN
        U1=REAL((MX(1)-M(1))/M(1)) ! FLUX VARIATION IN CELL
        F1=2.0*B*ABS(U1)
        IF(F1.LE.1) THEN
          P=1.0
          ISFIX=.TRUE.
        ELSE
          P=1.0/F1
          ISFIX=.FALSE.
        ENDIF
      ELSE
        P=0.0
        ISFIX=.FALSE.
      ENDIF

      ! CASE 0 <= P < 1
      ELSE
        ISFIX=.TRUE.
      ENDIF

      ! COMPUTE WEIGHTING FACTORS
      W(1)=-P*TB
      W(2)=1+P*TB
*----
* LINEAR ORDER ADAPTIVE CALCULATIONS
*----
      ELSEIF(IELEM.EQ.2) THEN
        CALL XABORT('SNADPT: LINEAR ORDER ADAPTIVE CALCULATIONS NOT'
     1   //'IMPLEMENTED YET')
      ELSE
        CALL XABORT('SNADPT: QUADRATIC AND HIGHER ORDER ADAPTIVE'
     1   //'CALCULATIONS NOT IMPLEMENTED YET')
      ENDIF

      RETURN
      END
