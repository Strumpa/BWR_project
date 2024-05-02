*DECK AUTPRD
      SUBROUTINE AUTPRD(NGRP,LBIN,NFS,SIGT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Cross section or source spreading.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NGRP    number of macro energy groups.
* LBIN    number of fine energy groups.
* NFS     number of fine energy groups in each coarse energy group.
* SIGT    cross section or source before spreading.
*
*Parameters: output
* SIGT    cross section or source after spreading.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NGRP,LBIN,NFS(NGRP)
      REAL SIGT(LBIN)
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: GAR
*
      ALLOCATE(GAR(NGRP))
      GAR(:NGRP)=SIGT(:NGRP)
      CALL XDRSET(SIGT,LBIN,0.0)
      IPO=LBIN
      DO J=NGRP,1,-1
        ND=ABS(NFS(J))
        SS=GAR(J)
        DO L=1,ND
          K=IPO-L+1
          SIGT(K)=SS
        ENDDO
      IPO=IPO-ND
      ENDDO
      DEALLOCATE(GAR)
      IF(IPO.NE.0) CALL XABORT('AUTPRD: SPREAD FAILURE.')
      RETURN
      END
