*DECK LIBNFI
      SUBROUTINE LIBNFI(IPLIB,NGRO,NBISO,NBMIX,NDEL,NESP,IPISO,MIX,
     1 MAXNFI,NFISS0,NFISSI,LSAME)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the maximum number of fissile isotopes in a mixture.
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
* IPLIB   pointer to the lattice microscopic cross section library
*         (L_LIBRARY signature).
* NGRO    number of energy groups.
* NBISO   number of isotopes present in the calculation domain.
* NBMIX   number of mixtures present in the calculation domain.
* NDEL    number of delayed precursor groups.
* NESP    number of energy-dependent fission spectra.
* IPISO   pointer array towards microlib isotopes.
* MIX     mixture number of each isotope (can be zero for void).
* MAXNFI  second dimension of array INDFIS.
*
*Parameters: output
* NFISS0  number of fissile isotopes in a mixture before adding new
*         fissile isotopes.
* NFISSI  number of fissile isotopes in a mixture after adding new
*         fissile isotopes.
* LSAME   fission spectrum mask (=.true. if all the isotopes have the
*         same fission spectrum and the same precursor group decay
*         constants.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB,IPISO(NBISO)
      INTEGER NGRO,NBISO,NBMIX,NDEL,NESP,MIX(NBISO),MAXNFI,NFISS0,NFISSI
      LOGICAL LSAME
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPLIB
      INTEGER NSTATE
      PARAMETER (NSTATE=40)
      CHARACTER HSMG*131
      INTEGER IDATA(NSTATE),ISOT,IBM,IFIS,IGR,ILONG,ITYLCM,IWFIS,KFIS,
     1 LENGT1,LENGT2,LENGTZ
      LOGICAL LFISS
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IWRK
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: INDFIS
      REAL, ALLOCATABLE, DIMENSION(:) :: CHI1,CHI2,LAM1,LAM2
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(INDFIS(NBMIX,MAXNFI),CHI1(NGRO),CHI2(NGRO),LAM1(NDEL),
     1 LAM2(NDEL))
*
      NFISS0=0
      NFISSI=0
      LSAME=(NESP.EQ.1)
      CALL LCMLEN(IPLIB,'MACROLIB',ILONG,ITYLCM)
      IF(ILONG.EQ.-1) THEN
         CALL LCMSIX(IPLIB,'MACROLIB',1)
         CALL LCMLEN(IPLIB,'SIGNATURE',ILONG,ITYLCM)
         IF(ILONG.EQ.0) GO TO 40
         CALL LCMGET(IPLIB,'STATE-VECTOR',IDATA)
         IF(IDATA(1).NE.NGRO) THEN
            WRITE(HSMG,'(38HLIBNFI: EXISTING MACROLIB HAVE NGROUP=,I4,
     1      26H NEW MACROLIB HAVE NGROUP=,I4,1H.)') IDATA(1),NGRO
            CALL XABORT(HSMG)
         ELSE IF(IDATA(2).GT.NBMIX) THEN
            WRITE(HSMG,'(37HLIBNFI: EXISTING MACROLIB HAVE NBMIX=,I4,
     1      25H NEW MACROLIB HAVE NBMIX=,I4,1H.)') IDATA(2),NBMIX
            CALL XABORT(HSMG)
         ELSE IF(IDATA(4).GT.MAXNFI*NESP) THEN
            WRITE(HSMG,'(38HLIBNFI: EXISTING MACROLIB HAVE NFISSI=,I4,
     1      13H GREATER THAN,I5,1H.)') IDATA(4)/NESP,MAXNFI
         ENDIF
         NFISSI=IDATA(4)/NESP
         NFISS0=NFISSI
         LSAME=(NESP.EQ.1).AND.(NFISSI.LE.1)
         IF(NFISSI.GT.0) THEN
            CALL LCMLEN(IPLIB,'FISSIONINDEX',ILONG,ITYLCM)
            IF(ILONG.EQ.0) THEN
*              THE NAMES ARE NOT DEFINED.
               INDFIS(:NBMIX,:NFISSI)=0
            ELSE IF(ILONG.EQ.NFISSI*NBMIX) THEN
               CALL LCMGET(IPLIB,'FISSIONINDEX',INDFIS)
            ELSE IF(ILONG.LT.NFISSI*NBMIX) THEN
*              REORDER THE 'FISSIONINDEX' MATRIX.
               ALLOCATE(IWRK(ILONG))
               CALL LCMGET(IPLIB,'FISSIONINDEX',IWRK)
               DO 30 IFIS=1,NFISSI
               DO 20 IBM=1,IDATA(2)
               INDFIS(IBM,IFIS)=IWRK((IFIS-1)*IDATA(2)+IBM)
   20          CONTINUE
               INDFIS(IDATA(2)+1:NBMIX,IFIS)=0
   30          CONTINUE
               DEALLOCATE(IWRK)
            ELSE
               WRITE(HSMG,'(40HLIBNFI: NUMBER OF FISSIONINDEX MIXTURES=,
     1         I6,30H; NUMBER OF MACROLIB MIXTURES=,I6,1H.)')
     2         ILONG/NFISSI,NBMIX
               CALL XABORT(HSMG)
            ENDIF
         ENDIF
   40    CALL LCMSIX(IPLIB,' ',2)
      ENDIF
      DO 120 ISOT=1,NBISO
      IBM=MIX(ISOT)
      IF(IBM.GT.0) THEN
         JPLIB=IPISO(ISOT)
         IF(C_ASSOCIATED(JPLIB)) THEN
            CALL LCMLEN(JPLIB,'NUSIGF',ILONG,ITYLCM)
            IF(NESP.EQ.1) THEN
               CALL LCMLEN(JPLIB,'CHI',LENGTZ,ITYLCM)
            ELSE
               CALL LCMLEN(JPLIB,'CHI--01',LENGTZ,ITYLCM)
            ENDIF
            IF((ILONG.GT.0).AND.(LENGTZ.GT.0)) THEN
               IF(NESP.EQ.1) THEN
                  CALL LCMGET(JPLIB,'CHI',CHI1)
               ELSE
                  CALL LCMGET(JPLIB,'CHI--01',CHI1)
               ENDIF
               LFISS=.FALSE.
               DO 50 IGR=1,NGRO
               LFISS=LFISS.OR.(CHI1(IGR).GT.0.0)
   50          CONTINUE
               IF(.NOT.LFISS) GO TO 120
               IF(LSAME) THEN
                  CALL LCMLEN(JPLIB,'LAMBDA-D',LENGT1,ITYLCM)
                  IF((LENGT1.EQ.NDEL).AND.(NDEL.GT.0)) THEN
                     CALL LCMGET(JPLIB,'LAMBDA-D',LAM1)
                  ENDIF
               ENDIF
               DO 60 IFIS=1,NFISSI
               IWFIS=INDFIS(IBM,IFIS)
               IF((IWFIS.EQ.ISOT).OR.(IWFIS.EQ.0)) THEN
                  KFIS=IFIS
                  GO TO 110
               ENDIF
   60          CONTINUE
               IF(LSAME) THEN
                  DO 90 IFIS=1,NFISSI
                  IWFIS=INDFIS(IBM,IFIS)
                  JPLIB=IPISO(IWFIS)
                  IF(NESP.EQ.1) THEN
                    CALL LCMGET(JPLIB,'CHI',CHI2)
                  ELSE
                    CALL LCMGET(JPLIB,'CHI--01',CHI2)
                  ENDIF
                  DO 70 IGR=1,NGRO
                  IF(ABS(CHI1(IGR)-CHI2(IGR)).GT.1.0E-3) THEN
                     LSAME=.FALSE.
                     GO TO 100
                  ENDIF
   70             CONTINUE
                  CALL LCMLEN(JPLIB,'LAMBDA-D',LENGT2,ITYLCM)
                  IF((LENGT1.EQ.NDEL).AND.(LENGT2.EQ.NDEL)
     1                               .AND.(NDEL.GT.0)) THEN
                     CALL LCMGET(JPLIB,'LAMBDA-D',LAM2)
                     DO 80 IGR=1,NDEL
                     IF(LAM1(IGR).NE.LAM2(IGR)) THEN
                        LSAME=.FALSE.
                        GO TO 100
                     ENDIF
   80                CONTINUE
                  ENDIF
   90             CONTINUE
               ENDIF
  100          NFISSI=NFISSI+1
               IF(NFISSI.GT.MAXNFI) CALL XABORT('LIBNFI: INDFIS OVERFL'
     1         //'OW.')
               KFIS=NFISSI
               INDFIS(:NBMIX,KFIS)=0
  110          INDFIS(IBM,KFIS)=ISOT
            ENDIF
         ENDIF
      ENDIF
  120 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(LAM2,LAM1,CHI2,CHI1,INDFIS)
      RETURN
      END
