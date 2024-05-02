*DECK USSONE
      SUBROUTINE USSONE(IPLI0,IPTRK,IPLIB,IFTRAK,CDOOR,IMPX,IGRMIN,
     1 IGRMAX,NIRES,NBNRS,IREX,NGRP,NBMIX,NREG,NUN,NBISO,NL,NED,NDEL,
     2 ISONAM,IHSUF,HCAL,DEN,MIX,IAPT,MAT,VOL,KEYFLX,LEAKSW,ITRANC,
     3 IPHASE,TITR,KSPH,ICORR,ISUBG,MAXST)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform a resonance self-shielding calculation named HCAL and build
* a corresponding internal library.
*
*Copyright:
* Copyright (C) 2003 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPLI0   pointer to the internal microscopic cross section library
*         builded by the self-shielding module (L_LIBRARY signature).
* IPTRK   pointer to the tracking. (L_TRACK signature).
* IPLIB   pointer to the internal microscopic cross section library
*         with subgroups (L_LIBRARY signature).
* IFTRAK  unit number of the sequential binary tracking file.
* CDOOR   name of the geometry/solution operator.
* IMPX    print flag (equal to zero for no print).
* IGRMIN  first group where the self-shielding is applied.
* IGRMAX  most thermal group where the self-shielding is applied.
* NIRES   number of correlated resonant isotopes in fuel regions.
* NBNRS   number of correlated fuel regions. Note that NBNRS=max(IREX).
* IREX    fuel region index assigned to each mixture. Equal to zero
*         in non-resonant mixtures or in mixtures not used.
* NGRP    number of energy groups.
* NBMIX   number of mixtures in the internal library.
* NREG    number of regions.
* NUN     number of unknowns per energy group.
* NBISO   number of isotopes specifications in the internal library.
* NL      number of Legendre orders required in the calculation
*         (NL=1 or higher).
* NED     number of extra vector edits.
* NDEL    number of delayed neutron precursor groups.
* ISONAM  alias name of isotopes.
* IHSUF   suffix name of isotopes.
* HCAL    name of the self-shielding calculation.
* DEN     density of each isotope.
* MIX     mix number of each isotope (can be zero).
* IAPT    resonant isotope index associated with isotope I. Mixed
*         moderator if IAPT(I)=NIRES+1. Out-of-fuel isotope if
*         IAPT(I)=0.
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* KEYFLX  pointers of fluxes in unknown vector.
* LEAKSW  leakage flag (LEAKSW=.TRUE. if neutron leakage through
*         external boundary is present).
* ITRANC  type of transport correction.
* IPHASE  type of flux solution (=1 use a native flux solution door;
*         =2 use collision probabilities).
* TITR    title.
* KSPH    SPH equivalence flag (=0 no SPH correction; =1 SPH correction
*         in the fuel).
* ICORR   mutual resonance shielding flag (=1 to suppress the model
*         in cases it is required in LIB operator).
* ISUBG   type of self-shielding model (=1 use physical probability
*         tables; =3 use original Ribon method; =4 use Ribon extended
*         method; =6 use resonance spectrum expansion method).
* MAXST   maximum number of fixed point iterations for the ST scattering
*         source.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLI0,IPTRK,IPLIB
      INTEGER IFTRAK,IMPX,IGRMIN,IGRMAX,NIRES,NBNRS,IREX(NBMIX),
     1 NGRP,NBMIX,NREG,NUN,NBISO,NL,NED,NDEL,ISONAM(3,NBISO),
     2 IHSUF(NBISO),MIX(NBISO),IAPT(NBISO),MAT(NREG),KEYFLX(NREG),
     3 ITRANC,IPHASE,KSPH,ICORR,ISUBG,MAXST
      REAL DEN(NBISO),VOL(NREG)
      LOGICAL LEAKSW
      CHARACTER CDOOR*12,HCAL*12,TITR*72
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISOBIS
      REAL, ALLOCATABLE, DIMENSION(:) :: SIGAR,UNGAR,DELTAU
      REAL, ALLOCATABLE, DIMENSION(:,:) :: GOLD
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SPH,PHGAR,STGAR,SFGAR,SWGAR
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: SSGAR,SAGAR,SDGAR
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: S0GAR
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: MASKI
      LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: MASKG
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(ISOBIS(3,NBISO))
      ALLOCATE(SPH(NBNRS,NIRES,NGRP),PHGAR(NBNRS,NIRES,NGRP),
     1 STGAR(NBNRS,NIRES,NGRP),SFGAR(NBNRS,NIRES,NGRP),
     2 SSGAR(NBNRS,NIRES,NL,NGRP),S0GAR(NBNRS,NIRES,NL,NGRP,NGRP),
     3 SAGAR(NBNRS,NIRES,NED,NGRP),SDGAR(NBNRS,NIRES,NDEL,NGRP),
     4 SWGAR(NBNRS,NIRES,NGRP),DELTAU(NGRP))
      ALLOCATE(MASKI(NBISO),MASKG(NGRP,NIRES))
*----
*  FIND THE NEW ISOTOPE NAMES IN IPLI0.
*----
      CALL LCMLEN(IPLI0,'ISOTOPESUSED',ILONG,ITYLCM)
      IF(ILONG.NE.0) THEN
         CALL LCMGET(IPLI0,'ISOTOPESUSED',ISOBIS)
      ELSE
         CALL LCMGET(IPLIB,'ISOTOPESUSED',ISOBIS)
      ENDIF
      DO 10 ISO=1,NBISO
      IF((IAPT(ISO).GT.0).AND.(IAPT(ISO).LE.NIRES)) THEN
         ISOBIS(3,ISO)=IHSUF(ISO)
      ENDIF
   10 CONTINUE
      CALL LCMPUT(IPLI0,'ISOTOPESUSED',3*NBISO,3,ISOBIS)
*
      ALLOCATE(SIGAR(4*NBMIX*(NIRES+1)*NGRP),UNGAR(NREG*NIRES*NGRP),
     1 GOLD(NIRES,NGRP))
*----
*  COMPUTE THE NEUTRON FLUX.
*----
      IF(ISUBG.EQ.6) THEN
        ! resonance spectrum expansion method
        CALL USSRSE(IPTRK,IPLIB,IPLI0,IFTRAK,NREG,NUN,NBMIX,NBISO,
     1  NIRES,NL,NED,NDEL,ISONAM,ISOBIS,HCAL,MAT,VOL,KEYFLX,CDOOR,
     2  LEAKSW,IMPX,DEN,MIX,IAPT,IPHASE,NGRP,IGRMIN,IGRMAX,NBNRS,
     3  IREX,TITR,ICORR,MAXST,GOLD,UNGAR,PHGAR,STGAR,SFGAR,SSGAR,
     4  S0GAR,SAGAR,SDGAR,MASKG,SIGAR)
      ELSE
        ! subgroup method
        CALL USSFLU(IPTRK,IPLIB,IPLI0,IFTRAK,NREG,NUN,NBMIX,NBISO,
     1  NIRES,NL,NED,NDEL,ISONAM,ISOBIS,HCAL,MAT,VOL,KEYFLX,CDOOR,
     2  LEAKSW,IMPX,DEN,MIX,IAPT,IPHASE,NGRP,IGRMIN,IGRMAX,NBNRS,
     3  IREX,TITR,ICORR,ISUBG,MAXST,GOLD,UNGAR,PHGAR,STGAR,SFGAR,
     4  SSGAR,S0GAR,SAGAR,SDGAR,SWGAR,MASKG,SIGAR)
      ENDIF
*----
*  COMPUTE THE SPH FACTORS.
*----
      SPH(:NBNRS,:NIRES,:NGRP)=1.0
      IF(KSPH.EQ.1) THEN
         CALL LCMGET(IPLI0,'DELTAU',DELTAU)
         CALL USSSPH(IPLI0,IPTRK,IFTRAK,NREG,NUN,NBMIX,NBISO,NIRES,NL,
     1   NED,NDEL,ISONAM,HCAL,MAT,VOL,KEYFLX,CDOOR,LEAKSW,IMPX,DEN,MIX,
     2   IAPT,ITRANC,IPHASE,NGRP,MASKG,NBNRS,IREX,TITR,ISUBG,SIGAR,
     3   GOLD,UNGAR,PHGAR,STGAR,SFGAR,SSGAR,S0GAR,SAGAR,SDGAR,SWGAR,
     4   DELTAU,SPH)
      ENDIF
*
      DEALLOCATE(GOLD,UNGAR,SIGAR)
*----
*  CREATE THE SELF-SHIELDED INTERNAL LIBRARY USING A SIMPLE
*  TRANSCRIPTION OF THE SELF-SHIELDED CROSS SECTIONS.
*----
      CALL KDRCPU(TK1)
*     SIMPLE TRANSCRIPTION OF THE SELF-SHIELDED CROSS SECTIONS.
      DO 100 ISO=1,NBISO
      MASKI(ISO)=(IAPT(ISO).GT.0).AND.(IAPT(ISO).LE.NIRES)
  100 CONTINUE
      DO 120 ISO=1,NBISO
      IF(MASKI(ISO)) THEN
         DO 110 JSO=ISO+1,NBISO
         IF((ISOBIS(1,ISO).EQ.ISOBIS(1,JSO)).AND.
     1      (ISOBIS(2,ISO).EQ.ISOBIS(2,JSO)).AND.
     2      (ISOBIS(3,ISO).EQ.ISOBIS(3,JSO))) MASKI(JSO)=.FALSE.
  110    CONTINUE
      ENDIF
  120 CONTINUE
      CALL USSIN1(IPLI0,IPLIB,NGRP,NBMIX,NBISO,NIRES,NBNRS,NL,NED,NDEL,
     1 IREX,IMPX,ISONAM,ISOBIS,MIX,IAPT,MASKI,SPH,PHGAR,STGAR,SFGAR,
     2 SSGAR,S0GAR,SAGAR,SDGAR)
      CALL KDRCPU(TK2)
      IF(IMPX.GT.1) WRITE(6,'(/36H USSONE: CPU TIME SPENT TO BUILD THE,
     1 33H SELF-SHIELDED INTERNAL LIBRARY =,F8.1,8H SECOND.)') TK2-TK1
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(MASKG,MASKI)
      DEALLOCATE(DELTAU,SWGAR,SDGAR,SAGAR,S0GAR,SSGAR,SFGAR,STGAR,
     1 PHGAR,SPH)
      DEALLOCATE(ISOBIS)
      RETURN
      END
