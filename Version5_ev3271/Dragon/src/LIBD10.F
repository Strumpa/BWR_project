*DECK LIBD10
#if defined(HDF5_LIB)
      SUBROUTINE LIBD10 (MAXDIL,NAMFIL,HNISOR,NDIL,DILUT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Find the dilutions corresponding to a resonant isotope within a
* library in Apolib-3 format.
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
* MAXDIL  maximum number of dilutions.
* NAMFIL  name of the APOLIB-3 file.
* HNISOR  library name of the resonant isotope.
*
*Parameters: output
* NDIL    number of finite dilutions.
* DILUT   dilutions.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE hdf5_wrap
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER MAXDIL,NDIL
      CHARACTER HNISOR*12,NAMFIL*(*)
      REAL DILUT(MAXDIL)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) IPAP
      CHARACTER RECNAM*80
      REAL, ALLOCATABLE, DIMENSION(:) :: BGXS
*----
*  OPEN AND PROBE THE APOLIB-3 FILE.
*----
      CALL hdf5_open_file(NAMFIL, IPAP, .TRUE.)
      WRITE(RECNAM,'(9HIsotopes/,A,15H/HomoRates/BgXS)') TRIM(HNISOR)
      CALL hdf5_read_data(IPAP, RECNAM, BGXS)
      NDIL=SIZE(BGXS)
      IF(NDIL.GT.MAXDIL) CALL XABORT('LIBD10: MAXDIL IS TOO SMALL.')
      DILUT(:NDIL)=BGXS(:)
      NDIL=NDIL-1
      DEALLOCATE(BGXS)
      CALL hdf5_close_file(IPAP)
#else
      SUBROUTINE LIBD10
      CALL XABORT('LIBD10: THE HDF5 API IS NOT AVAILABLE.')
#endif /* defined(HDF5_LIB) */
      RETURN
      END
