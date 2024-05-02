*DECK USSIT3
      SUBROUTINE USSIT3(MAXNOR,NGRP,MASKG,IRES,IPLI0,IPTRK,IFTRAK,CDOOR,
     1 IMPX,NBMIX,NREG,NUN,IPHASE,MAXST,MAT,VOL,KEYFLX,LEAKSW,IREX,
     2 SIGGAR,TITR,ICORR,NIRES,NBNRS,CONR,GOLD,IPPT1,IPPT2,VOLMER,
     3 DELTAU,UNGAR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the snapshot weights as required by the resonance spectrum
* expansion (RSE) method:
* a) assume a single resonant isotope;
* b) use the standard solution doors of Dragon.
*
*Copyright:
* Copyright (C) 2023 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* MAXNOR  maximum number of base points.
* NGRP    number of energy group.
* MASKG   energy group mask pointing on self-shielded groups.
* IRES    index of the resonant isotope.
* IPLI0   pointer to the internal microscopic cross section library
*         builded by the self-shielding module.
* IPTRK   pointer to the tracking (L_TRACK signature).
* IFTRAK  file unit number used to store the tracks.
* CDOOR   name of the geometry/solution operator.
* IMPX    print flag (equal to zero for no print).
* NBMIX   number of mixtures in the internal library.
* NREG    number of regions.
* NUN     number of unknowns in the flux or source vector in one
*         energy group and one band.
* IPHASE  type of flux solution (=1 use a native flux solution door;
*         =2 use collision probabilities).
* MAXST   maximum number of fixed point iterations for the ST scattering
*         source.
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* KEYFLX  pointers of fluxes in unknown vector.
* LEAKSW  leakage switch (LEAKSW=.TRUE. if neutron leakage through
*         external boundary is present).
* IREX    fuel region index assigned to each mixture. Equal to zero
*         in non-resonant mixtures or in mixtures not used.
* SIGGAR  macroscopic x-s of the non-resonant isotopes in each mixture:
*         (*,*,*,1) total; (*,*,*,2) transport correction; 
*         (*,*,*,3) P0 scattering.
* TITR    title.
* ICORR   mutual resonance shielding flag (=1 to suppress the model
*         in cases it is required in LIB operator).
* NIRES   exact number of correlated resonant isotopes.
* NBNRS   number of correlated fuel regions.
* CONR    number density of the resonant isotopes.
* GOLD    type of self-shielding model (=1.0 physical probability
*         tables; =-1001.0 resonance spectrum expansion method).
* IPPT1   pointer to LCM directory of each resonant isotope.
* IPPT2   information related to each resonant isotope:
*         IPPT2(:,1) index of a resonant region (used with infinite
*         dilution case);
*         IPPT2(:,2:4) alias name of resonant isotope.
* VOLMER  volumes of the resonant regions.
* DELTAU  lethargy widths of coarse groups.
*
*Parameters: output
* UNGAR   averaged fluxes per volume.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLI0,IPTRK,IPPT1(NIRES)
      INTEGER MAXNOR,NGRP,IRES,IFTRAK,IMPX,NBMIX,NREG,NUN,IPHASE,
     1 MAXST,MAT(NREG),KEYFLX(NREG),IREX(NBMIX),ICORR,NIRES,NBNRS,
     2 IPPT2(NIRES,4)
      REAL VOL(NREG),SIGGAR(NBMIX,0:NIRES,NGRP,3),CONR(NBNRS,NIRES),
     1 GOLD(NIRES,NGRP),VOLMER(0:NBNRS),UNGAR(NREG,NIRES,NGRP),
     2 DELTAU(NGRP)
      CHARACTER CDOOR*12,TITR*72
      LOGICAL LEAKSW,MASKG(NGRP)
*----
*  LOCAL VARIABLES
*----
      REAL ERR1,ERR2
      DOUBLE PRECISION T1,DXI
      CHARACTER CBDPNM*12,TEXT12*12
      LOGICAL LEXAC,REBFLG
      TYPE(C_PTR) IPLIB,JPLI0,JPLIB1,KPLIB,IPSYS,KPSYS,IOFSET,
     1 IPMACR,IPSOU
*----
*  ALLOCATABLE ARRAYS
*----
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) ::  JPLIB2
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NJJ,NPSYS,MRANK
      REAL, ALLOCATABLE, DIMENSION(:) :: SIGTXS,SIGS0X
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FUN,SUN,XFLUX
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: XFLUX2
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: SIGTAD
      TYPE VECTOR_ARRAY
        DOUBLE PRECISION, POINTER, DIMENSION(:) :: VECTOR
      END TYPE VECTOR_ARRAY
      TYPE MATRIX_ARRAY
        DOUBLE PRECISION, POINTER, DIMENSION(:,:) :: MATRIX
      END TYPE MATRIX_ARRAY
      TYPE(VECTOR_ARRAY), ALLOCATABLE, DIMENSION(:) :: SIGT_V,XI_V,
     1 GAMMA_V
      TYPE(MATRIX_ARRAY), ALLOCATABLE, DIMENSION(:,:) :: SIGT_M
      TYPE(MATRIX_ARRAY), ALLOCATABLE, DIMENSION(:,:,:) :: SCAT_M
      TYPE MATRIX_ARRAY_SP
        REAL, POINTER, DIMENSION(:,:) :: MATRIX
      END TYPE MATRIX_ARRAY_SP
      TYPE(MATRIX_ARRAY_SP), ALLOCATABLE, DIMENSION(:) :: PSI_M
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(JPLIB2(NIRES))
      ALLOCATE(NJJ(NGRP),NPSYS(MAXNOR*NGRP),MRANK(NGRP))
      ALLOCATE(SIGT_V(NGRP),SIGT_M(NGRP,NIRES),SCAT_M(NGRP,NGRP,NIRES),
     1 XI_V(NGRP),GAMMA_V(NGRP),PSI_M(NGRP),SIGTAD(NIRES,NGRP))
      SIGTAD(:NIRES,:NGRP)=0.0D0
*----
*  FIND THE NUMBER OF COMPONENTS REQUIRED AND ALLOCATE THE LIST OF
*  ASSEMBLY MATRICES.
*----
      IPLIB=IPPT1(IRES)
      CALL LCMLEN(IPLIB,'NOR',ILONG,ITYLCM)
      IF(ILONG.NE.NGRP) THEN
        CALL LCMLIB(IPLIB)
        CALL XABORT('USSIT3: RANK ARRAY MISSING.')
      ENDIF
      CALL LCMGET(IPLIB,'NOR',MRANK)
      CALL LCMGET(IPLIB,'NJJS00',NJJ)
      JPLIB1=LCMGID(IPLIB,'GROUP-RSE')
      NASM=0
      DO IG=1,NGRP
        IF(MASKG(IG).AND.(GOLD(IRES,IG).EQ.-1001.)) THEN
          NASM=NASM+MRANK(IG)
        ENDIF
      ENDDO
      IF(NASM.EQ.0) GO TO 50
      DO JRES=1,NIRES
        DO JG=1,NGRP
          DO IG=1,NGRP
            NULLIFY(SCAT_M(IG,JG,JRES)%MATRIX)
          ENDDO
        ENDDO
      ENDDO
*----
*  CREATE A SPECIFIC DIRECTORY FOR IRES-TH RESONANT ISOTOPE.
*----
      WRITE(CBDPNM,'(3HCOR,I4.4,1H/,I4.4)') IRES,NIRES
      CALL LCMSIX(IPLI0,CBDPNM,1)
      JPLI0=LCMGID(IPLI0,'NWT0-PT')
      IPSYS=LCMLID(IPLI0,'ASSEMB-RSE',NASM)
      CALL LCMSIX(IPLI0,' ',2)
*----
*  RECOVER RSE INFORMATION FROM MICROLIB (PART 1)
*----
      DO JRES=1,NIRES
        WRITE(TEXT12,'(3A4)') (IPPT2(JRES,I),I=2,4)
        JPLIB2(JRES)=LCMGID(IPLIB,TEXT12) ! SCAT_M information
      ENDDO
*----
*  LOOP OVER ENERGY GROUPS FOR THE ASSEMBLY CALCULATION.
*----
      IPOS=1
      DO IG=1,NGRP
        IF(.NOT.MASKG(IG).OR.(GOLD(IRES,IG).NE.-1001.)) THEN
          IPOS=IPOS+NJJ(IG)
          CYCLE
        ENDIF
        IF(IMPX.GE.9) WRITE(6,'(22H USSIT3: energy group=,I8)') IG
        MI=MRANK(IG)
*----
*  RECOVER RSE INFORMATION FROM MICROLIB
*----
        KPLIB=LCMGIL(JPLIB1,IG)
        CALL LCMLEN(KPLIB,'SIGT_V',ILONG,ITYLCM)
        IF(ILONG.GT.MAXNOR) CALL XABORT('USSIT3: MAXNOR OVERFLOW.')
        CALL LCMGPD(KPLIB,'SIGT_V',IOFSET)
        CALL C_F_POINTER(IOFSET,SIGT_V(IG)%VECTOR,(/MI/))
        CALL LCMGPD(KPLIB,'XI_V',IOFSET)
        CALL C_F_POINTER(IOFSET,XI_V(IG)%VECTOR,(/MI/))
        CALL LCMGPD(KPLIB,'GAMMA_V',IOFSET)
        CALL C_F_POINTER(IOFSET,GAMMA_V(IG)%VECTOR,(/MI/))
        DO JRES=1,NIRES
          IF(JRES.NE.IRES) THEN
            WRITE(TEXT12,'(3A4)') (IPPT2(JRES,I),I=2,4)
            CALL LCMGPD(KPLIB,TEXT12,IOFSET)
            CALL C_F_POINTER(IOFSET,SIGT_M(IG,JRES)%MATRIX,(/MI,MI/))
          ENDIF
          DO JG=IG-NJJ(IG)+1,IG
            MJ=MRANK(JG)
            CALL LCMGPL(JPLIB2(JRES),IPOS+IG-JG,IOFSET)
            CALL C_F_POINTER(IOFSET,SCAT_M(IG,JG,JRES)%MATRIX,(/MI,MJ/))
          ENDDO
        ENDDO
        IPOS=IPOS+NJJ(IG)
      ENDDO
*----
*  INITIALIZE THE SUBGROUP FLUX WITH FUNKNO$USS INFORMATION
*----
      IASM=0
      DO IG=1,NGRP
        IF(.NOT.MASKG(IG).OR.(GOLD(IRES,IG).NE.-1001.)) CYCLE
        MI=MRANK(IG)
        ALLOCATE(PSI_M(IG)%MATRIX(NUN,MI))
        DO IM=1,MI
          CALL LCMLEL(IPSYS,IASM+IM,ILONG,ITYLCM)
          IF(ILONG.EQ.-1) THEN
            KPSYS=LCMGIL(IPSYS,IASM+IM)
            CALL LCMGET(KPSYS,'FUNKNO$USS',PSI_M(IG)%MATRIX(:NUN,IM))
          ELSE
            PSI_M(IG)%MATRIX(:NUN,IM)=REAL(GAMMA_V(IG)%VECTOR(IM))
          ENDIF
        ENDDO
*----
*  COMPUTE NTOT0 AND SIGW00 CROSS SECTION OF ADMIXED RESONANT ISOTOPES.
*----
        IF(NIRES.GT.1) THEN
          ALLOCATE(XFLUX(NBNRS,MI))
          CALL LCMLEL(JPLI0,IG,ILONG,ITYLCM)
          IF(ILONG.EQ.0) CALL XABORT('USSIT3: MISSING NWT0-PT RECORD.')
          CALL LCMGDL(JPLI0,IG,XFLUX)
          T1=0.0D0
          DO JRES=1,NIRES
            IND=IPPT2(IRES,1)
            IF(JRES.NE.IRES) THEN
              DO IM=1,MI
                DXI=XI_V(IG)%VECTOR(IM)
                T1=T1+DXI*XFLUX(IND,IM)
                DO JM=1,MI
                  IF(JM.EQ.IM) CYCLE
                  SIGTAD(JRES,IG)=SIGTAD(JRES,IG)+DXI*
     1            SIGT_M(IG,JRES)%MATRIX(IM,JM)*XFLUX(IND,JM)
                ENDDO
              ENDDO
              SIGTAD(JRES,IG)=SIGTAD(JRES,IG)/T1
            ENDIF
          ENDDO
          DEALLOCATE(XFLUX)
        ENDIF
*----
*  COMPUTE GROUPWISE MACROSCOPIC CROSS SECTIONS.
*----
        ALLOCATE(SIGTXS(0:NBMIX),SIGS0X(0:NBMIX))
        DO IM=1,MI
          SIGTXS(0:NBMIX)=0.0
          SIGS0X(0:NBMIX)=0.0
          DO IBM=1,NBMIX
            IND=IREX(IBM)
            DO 10 JRES=0,NIRES
            IF(JRES.EQ.0) THEN
*             ADMIXED NON-RESONANT ISOTOPES.
              SIGTXS(IBM)=SIGTXS(IBM)+(SIGGAR(IBM,0,IG,1)-
     1        SIGGAR(IBM,0,IG,2))
              SIGS0X(IBM)=SIGS0X(IBM)-SIGGAR(IBM,0,IG,2)
            ELSE IF((JRES.NE.IRES).AND.(IND.GT.0).AND.(ICORR.EQ.1)) THEN
*             ECCO CORRELATION MODEL.
              IF((IPPT2(IRES,2).EQ.IPPT2(JRES,2)).AND.
     1           (IPPT2(IRES,3).EQ.IPPT2(JRES,3))) THEN
                DENSIT=CONR(IND,JRES)
                SIGTXS(IBM)=SIGTXS(IBM)+DENSIT*
     1                 REAL(SIGT_V(IG)%VECTOR(IM))
                SIGS0X(IBM)=SIGS0X(IBM)+DENSIT*
     1                 REAL(SCAT_M(IG,IG,JRES)%MATRIX(IM,IM))
              ENDIF
            ELSE IF((JRES.NE.IRES).AND.(IND.GT.0).AND.(ICORR.EQ.0)) THEN
*             MUTUAL SHIELDING MODEL OF CORRELATED RESONANT ISOTOPES.
              DENSIT=CONR(IND,JRES)
              SIGTXS(IBM)=SIGTXS(IBM)+DENSIT*REAL(SIGTAD(JRES,IG))
              SIGTXS(IBM)=SIGTXS(IBM)+DENSIT*
     1               REAL(SIGT_M(IG,JRES)%MATRIX(IM,IM))
              SIGS0X(IBM)=SIGS0X(IBM)+DENSIT*
     1               REAL(SCAT_M(IG,IG,JRES)%MATRIX(IM,IM))
            ENDIF
   10       CONTINUE
            IF(IND.GT.0) THEN
              DENSIT=CONR(IND,IRES)
              SIGTXS(IBM)=SIGTXS(IBM)+DENSIT*REAL(SIGT_V(IG)%VECTOR(IM))
              SIGS0X(IBM)=SIGS0X(IBM)+DENSIT*
     1               REAL(SCAT_M(IG,IG,IRES)%MATRIX(IM,IM))
            ENDIF
          ENDDO
          NPSYS(IASM+IM)=IASM+IM
          KPSYS=LCMDIL(IPSYS,IASM+IM)
          CALL LCMPUT(KPSYS,'DRAGON-TXSC',NBMIX+1,2,SIGTXS(0))
          CALL LCMPUT(KPSYS,'DRAGON-S0XSC',NBMIX+1,2,SIGS0X(0))
        ENDDO
        IASM=IASM+MI
        DEALLOCATE(SIGS0X,SIGTXS)
      ENDDO
*----
*  ASSEMBLY MATRIX OR REDUCED COLLISION PROBABILITIES CALCULATION.
*----
      NANI=1
      KNORM=1
      NALBP=0
      IMPY=MAX(0,IMPX-3)
      IF(IPHASE.EQ.1) THEN
*        USE A NATIVE DOOR.
         ISTRM=1
         NW=0
         CALL DOORAV(CDOOR,IPSYS,NPSYS,IPTRK,IFTRAK,IMPY,NASM,NREG,
     1   NBMIX,NANI,NW,MAT,VOL,KNORM,LEAKSW,TITR,NALBP,ISTRM)
      ELSE IF(IPHASE.EQ.2) THEN
*        USE A COLLISION PROBABILITY DOOR.
         IPIJK=1
         CALL DOORPV(CDOOR,IPSYS,NPSYS,IPTRK,IFTRAK,IMPY,NASM,NREG,
     1   NBMIX,NANI,MAT,VOL,KNORM,IPIJK,LEAKSW,.FALSE.,TITR,NALBP)
      ENDIF
*----
*  LOOP OVER ENERGY GROUPS FOR THE FLUX CALCULATION.
*----
      ALLOCATE(XFLUX2(NBNRS,MAXNOR,NGRP))
      XFLUX2(:NBNRS,:MAXNOR,:NGRP)=0.0
      IASM=0
      DO IG=1,NGRP
        MI=MRANK(IG)
        IF(.NOT.MASKG(IG).OR.(GOLD(IRES,IG).NE.-1001.)) CYCLE
        ITER=0
   20   ITER=ITER+1
        IF(ITER.GT.MAXST) GO TO 30
        ERR1=0.0
        ERR2=0.0
*----
*  COMPUTE THE AVERAGED SOURCE TAKING INTO ACCOUNT CORRELATION EFFECTS.
*----
        ALLOCATE(FUN(NUN,MI),SUN(NUN,MI))
        SUN(:NUN,:MI)=0.0
        DO IM=1,MI
          FUN(:NUN,IM)=PSI_M(IG)%MATRIX(:NUN,IM)
          NPSYS(IM)=IASM+IM
          DO I=1,NREG
            IBM=MAT(I)
            IF(IBM.EQ.0) CYCLE
            IUN=KEYFLX(I)
            IND=IREX(IBM)
            T1=0.0D0
            DO JRES=0,NIRES
              IF(JRES.EQ.0) THEN
                T1=T1+SIGGAR(IBM,0,IG,3)*GAMMA_V(IG)%VECTOR(IM)
              ELSE IF(IND.GT.0) THEN
                DENSIT=CONR(IND,JRES)
                DO JG=IG-NJJ(IG)+1,IG
                  IF(GOLD(IRES,JG).NE.-1001.) CYCLE
                  DO JM=1,MRANK(JG)
                    IF((JG.EQ.IG).AND.(JM.EQ.IM)) CYCLE
                    T1=T1+DENSIT*SCAT_M(IG,JG,JRES)%MATRIX(IM,JM)*
     1              PSI_M(JG)%MATRIX(IUN,JM)
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
            SUN(IUN,IM)=REAL(T1,4)
          ENDDO
        ENDDO
*----
*  SOLVE FOR THE MULTIBAND FLUX.
*----
        IDIR=0
        LEXAC=.FALSE.
        IPMACR=C_NULL_PTR
        IPSOU=C_NULL_PTR
        REBFLG=.FALSE.
        CALL DOORFV(CDOOR,IPSYS,NPSYS,IPTRK,IFTRAK,IMPX,MI,NBMIX,IDIR,
     1  NREG,NUN,IPHASE,LEXAC,MAT,VOL,KEYFLX,TITR,SUN,FUN,IPMACR,IPSOU,
     2  REBFLG)
*----
*  CONVERGENCE CONTROL.
*----
        DO IM=1,MI
          KPSYS=LCMGIL(IPSYS,IASM+IM)
          CALL LCMPUT(KPSYS,'FUNKNO$USS',NUN,2,FUN(1,IM))
          DO I=1,NREG
            IUN=KEYFLX(I)
            DELTA=FUN(IUN,IM)-PSI_M(IG)%MATRIX(IUN,IM)
            ERR1=MAX(ERR1,ABS(DELTA))
            ERR2=MAX(ERR2,ABS(FUN(IUN,IM)))
          ENDDO
          PSI_M(IG)%MATRIX(:NUN,IM)=FUN(:NUN,IM)
        ENDDO
        DEALLOCATE(SUN,FUN)
        IF(IMPX.GT.2) THEN
          WRITE(TEXT12,'(3A4)') (IPPT2(IRES,I),I=2,4)
          WRITE(6,'(15H USSIT3: GROUP=,I5,15H. RSE ITERATION,I4,
     1    11H. ISOTOPE='',A12,9H''. ERROR=,1P,E11.4,1H.)') IG,
     2    ITER,TEXT12,ERR1
        ENDIF
        IF(ERR1.GT.1.0E5) GO TO 30
        IF(ERR1.GT.1.0E-4*ERR2) GO TO 20
        IF(IMPX.GT.1) THEN
           WRITE(TEXT12,'(3A4)') (IPPT2(IRES,I),I=2,4)
           WRITE(6,'(15H USSIT3: GROUP=,I5,24H. RSE ITERATION CONVERGE,
     1     6HNCE IN,I4,22H ITERATIONS. ISOTOPE='',A12,2H''.)') IG,
     2     ITER,TEXT12
        ENDIF
*----
*  COMPUTE XFLUX2 FOR IRES IN GROUP IG.
*----
        XFLUX2(:NBNRS,:MI,IG)=0.0
        DO I=1,NREG
          IF(MAT(I).EQ.0) CYCLE
          IND=IREX(MAT(I))
          IF(IND.EQ.0) CYCLE
          IUN=KEYFLX(I)
          DO IM=1,MI
            XFLUX2(IND,IM,IG)=XFLUX2(IND,IM,IG)+VOL(I)*
     1      PSI_M(IG)%MATRIX(IUN,IM)
          ENDDO
        ENDDO
        DO IM=1,MI
          DO IND=1,NBNRS
            XFLUX2(IND,IM,IG)=XFLUX2(IND,IM,IG)/VOLMER(IND)
          ENDDO
        ENDDO
*----
* USE SNAPSHOT WEIGHTS AND HOMOGENIZE THE FLUX.
*----
        UNGAR(:NREG,IRES,IG)=0.0
        DO I=1,NREG
          IF(MAT(I).EQ.0) CYCLE
          IUN=KEYFLX(I)
          DO IM=1,MI
            UNGAR(I,IRES,IG)=UNGAR(I,IRES,IG)+REAL(XI_V(IG)%VECTOR(IM)*
     1      PSI_M(IG)%MATRIX(IUN,IM),4)/DELTAU(IG)
          ENDDO
        ENDDO
        GO TO 40
*----
*  ALTERNATIVE TREATMENT IN CASE OF FAILURE OF FIXED POINT ITERATIONS.
*  USE A NON-ITERATIVE RESPONSE MATRIX APPROACH.
*----
   30   IF(IMPX.GT.0) THEN
           WRITE(TEXT12,'(3A4)') (IPPT2(IRES,I),I=2,4)
           WRITE(6,'(15H USSIT3: GROUP=,I5,24H. SUBGROUP ITERATION FAI,
     1     16HLED FOR ISOTOPE ,A12,32H. USE AN ALTERNATIVE RESPONSE MA,
     2     14HTRIX APPROACH.)') IG,TEXT12
        ENDIF
        CALL USSEXD(MAXNOR,CDOOR,IPLI0,IPTRK,IFTRAK,IMPX,NGRP,IG,IASM,
     1  NBMIX,NREG,NUN,IPHASE,MAT,VOL,KEYFLX,IREX,SIGGAR,NJJ,TITR,
     2  NIRES,IRES,NBNRS,MRANK,CONR,GOLD,DELTAU,IPPT1,IPPT2,VOLMER,
     3  XFLUX2,UNGAR)
*----
* SAVE XFLUX2 FOR IRES IN GROUP IG.
*----
   40   CALL LCMPDL(JPLI0,IG,NBNRS*MI,2,XFLUX2(1,1,IG))
        IF(IMPX.GT.2) THEN
          DO IND=1,NBNRS
            T1=0.0D0
            DO IM=1,MI
              T1=T1+XI_V(IG)%VECTOR(IM)*XFLUX2(IND,IM,IG)/DELTAU(IG)
            ENDDO
            WRITE(6,'(31H USSIT3: AVERAGED FLUX IN GROUP,I4,9H AND RESO,
     1      11HNANT REGION,I4,21H FOR RESONANT ISOTOPE,I4,2H =,F9.5)')
     2      IG,IND,IRES,T1
          ENDDO
        ENDIF
        IASM=IASM+MI
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION.
*----
      DEALLOCATE(XFLUX2)
      DO IG=1,NGRP
        IF(.NOT.MASKG(IG).OR.(GOLD(IRES,IG).NE.-1001.)) CYCLE
        DEALLOCATE(PSI_M(IG)%MATRIX)
      ENDDO
   50 DEALLOCATE(SIGTAD,PSI_M,GAMMA_V,XI_V,SCAT_M,SIGT_M,SIGT_V)
      DEALLOCATE(MRANK,NPSYS,NJJ)
      DEALLOCATE(JPLIB2)
      RETURN
      END
