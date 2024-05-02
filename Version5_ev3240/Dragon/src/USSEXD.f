*DECK USSEXD
      SUBROUTINE USSEXD(MAXNOR,CDOOR,IPLI0,IPTRK,IFTRAK,IMPX,NGRP,IG,
     1 IASM,NBMIX,NREG,NUN,IPHASE,MAT,VOL,KEYFLX,IREX,SIGGAR,NJJ,TITR,
     2 NIRES,IRES,NBNRS,MRANK,CONR,GOLD,DELTAU,IPPT1,IPPT2,VOLMER,
     3 XFLUX,UNGAR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solution of the flux for the resonance spectrum expansion (RSE) method
* using the response matrix method. This is a non-iterative approach
* which is useful in exceptional cases where the fixed-point approach
* fails.
*
*Copyright:
* Copyright (C) 2024 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* MAXNOR  maximum order of the probability tables (RSE).
* CDOOR   name of the geometry/solution operator.
* IPLI0   pointer to the internal microscopic cross section library
*         builded by the self-shielding module.
* IPTRK   pointer to the tracking (L_TRACK signature).
* IFTRAK  file unit number used to store the tracks.
* IMPX    print flag (equal to zero for no print).
* NGRP    number of energy groups.
* IG      index of energy group being processed.
* IASM    offset of information computed by DOORAV or DOORPV.
* NBMIX   number of mixtures in the internal library.
* NREG    number of regions.
* NUN     number of unknowns in the flux or source vector in one
*         energy group and one band.
* IPHASE  type of flux solution (=1 use a native flux solution door;
*         =2 use collision probabilities).
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* KEYFLX  pointers of fluxes in unknown vector.
* IREX    fuel region index assigned to each mixture. Equal to zero
*         in non-resonant mixtures or in mixtures not used.
* SIGGAR  macroscopic x-s of the non-resonant isotopes in each mixture:
*         (*,*,*,1) total; (*,*,*,2) transport correction; 
*         (*,*,*,3) P0 scattering.
* NJJ     scattering bandwidth.
* TITR    title.
* NIRES   exact number of correlated resonant isotopes.
* IRES    index of the resonant isotope being processed.
* NBNRS   number of correlated fuel regions.
* MRANK   exact order of the probability table.
* CONR    number density of the resonant isotopes.
* GOLD    type of self-shielding model (=1.0 physical probability
*         tables; =-1001.0 resonance spectrum expansion method).
* IPPT1   pointer to LCM directory of each resonant isotope.
* IPPT2   information related to each resonant isotope:
*         IPPT2(:,1) index of a resonant region (used with infinite
*         dilution case);
*         IPPT2(:,2:4) alias name of resonant isotope.
* VOLMER  volumes of the resonant and non-resonant regions.
*
*Parameters: input/output
* XFLUX   subgroup flux.
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
      INTEGER MAXNOR,IFTRAK,IMPX,NGRP,IG,IASM,NBMIX,NREG,NUN,IPHASE,
     1 MAT(NREG),KEYFLX(NREG),IREX(NBMIX),NJJ(NGRP),NIRES,IRES,NBNRS,
     2 MRANK(NGRP),IPPT2(NIRES,4)
      REAL VOL(NREG),SIGGAR(NBMIX,0:NIRES,NGRP,3),CONR(NBNRS,NIRES),
     1 GOLD(NIRES,NGRP),DELTAU(NGRP),VOLMER(0:NBNRS),
     2 XFLUX(NBNRS,MAXNOR,NGRP),UNGAR(NREG,NIRES,NGRP)
      CHARACTER CDOOR*12,TITR*72
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) IPSYS,KPSYS,IPMACR,IPSOU,IPLIB,JPLIB1,KPLIB,IOFSET
      DOUBLE PRECISION QQQ,SSS,T1
      LOGICAL LEXAC,REBFLG
      CHARACTER CBDPNM*12,TEXT12*12
      TYPE VECTOR_ARRAY
        DOUBLE PRECISION, POINTER, DIMENSION(:) :: VECTOR
      END TYPE VECTOR_ARRAY
      TYPE(VECTOR_ARRAY) :: XI_V,GAMMA_V
*----
*  ALLOCATABLE ARRAYS
*----
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) ::  JPLIB2
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NPSYS
      REAL, ALLOCATABLE, DIMENSION(:) :: AWPHI,FUN,SUN
      REAL, ALLOCATABLE, DIMENSION(:,:) :: PAV
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: MATRIX
      TYPE MATRIX_ARRAY
        DOUBLE PRECISION, POINTER, DIMENSION(:,:) :: MATRIX
      END TYPE MATRIX_ARRAY
      TYPE(MATRIX_ARRAY), ALLOCATABLE, DIMENSION(:,:) :: SCAT_M
*----
*  STATEMENT FUNCTIONS
*----
      INM(IND,IM,NBNRS)=(IM-1)*NBNRS+IND
*----
*  SCRATCH STORAGE ALLOCATION
*----
      MI=MRANK(IG)
      ALLOCATE(PAV(0:NBNRS,0:NBNRS),AWPHI(0:NBNRS))
      ALLOCATE(MATRIX(NBNRS*MI,NBNRS*MI+1))
      ALLOCATE(JPLIB2(NIRES),SCAT_M(NGRP,NIRES))
*----
*  RECOVER RSE INFORMATION FROM MICROLIB
*----
      IPOS=1
      DO JG=1,IG-1
        IPOS=IPOS+NJJ(JG)
      ENDDO
      IPLIB=IPPT1(IRES)
      JPLIB1=LCMGID(IPLIB,'GROUP-RSE')
      DO JRES=1,NIRES
        WRITE(TEXT12,'(3A4)') (IPPT2(JRES,I),I=2,4)
        JPLIB2(JRES)=LCMGID(IPLIB,TEXT12) ! SCAT_M information
      ENDDO
      KPLIB=LCMGIL(JPLIB1,IG)
      CALL LCMGPD(KPLIB,'XI_V',IOFSET)
      CALL C_F_POINTER(IOFSET,XI_V%VECTOR,(/MI/))
      CALL LCMGPD(KPLIB,'GAMMA_V',IOFSET)
      CALL C_F_POINTER(IOFSET,GAMMA_V%VECTOR,(/MI/))
      DO JRES=1,NIRES
        DO JG=IG-NJJ(IG)+1,IG
          MJ=MRANK(JG)
          CALL LCMGPL(JPLIB2(JRES),IPOS+IG-JG,IOFSET)
          CALL C_F_POINTER(IOFSET,SCAT_M(JG,JRES)%MATRIX,(/MI,MJ/))
        ENDDO
      ENDDO
*----
*  RECOVER THE SPECIFIC DIRECTORY FOR IRES-TH RESONANT ISOTOPE.
*----
      WRITE(CBDPNM,'(3HCOR,I4.4,1H/,I4.4)') IRES,NIRES
      CALL LCMSIX(IPLI0,CBDPNM,1)
      IPSYS=LCMGID(IPLI0,'ASSEMB-RSE')
      CALL LCMSIX(IPLI0,' ',2)
*----
*  COMPUTE THE AVERAGED COLLISION PROBABILITY MATRIX.
*----
      ALLOCATE(NPSYS(MI*(NBNRS+1)))
      ALLOCATE(FUN(NUN*MI*(NBNRS+1)),SUN(NUN*MI*(NBNRS+1)))
      CALL XDRSET(FUN,NUN*MI*(NBNRS+1),0.0)
      CALL XDRSET(SUN,NUN*MI*(NBNRS+1),0.0)
      DO 50 IM=1,MI
      DO 40 JNBN=0,NBNRS
      NPSYS((IM-1)*(NBNRS+1)+JNBN+1)=IASM+IM
      T1=0.0D0
      DO 20 I=1,NREG
      IBM=MAT(I)
      IF(IBM.EQ.0) GO TO 20
      IOF=(IM-1)*NUN*(NBNRS+1)+JNBN*NUN+KEYFLX(I)-1
      IND=IREX(IBM)
      IF((JNBN.EQ.0).AND.(IND.EQ.0)) THEN
         SSS=SIGGAR(IBM,0,IG,3)*GAMMA_V%VECTOR(IM)
         T1=T1+SSS*VOL(I)
         SUN(IOF+1)=SUN(IOF+1)+REAL(SSS,4)
      ELSE IF(IND.EQ.JNBN) THEN
         T1=T1+VOL(I)
         SUN(IOF+1)=1.0
      ENDIF
   20 CONTINUE
      DO 30 I=1,NUN
      IOF=(IM-1)*NUN*(NBNRS+1)+JNBN*NUN+I-1
      IF(T1.NE.0.0) SUN(IOF+1)=SUN(IOF+1)/REAL(T1,4)
   30 CONTINUE
   40 CONTINUE
   50 CONTINUE
*----
*  SOLVE FOR THE MULTIBAND FLUX.
*----
      IDIR=0
      NABS=MI*(NBNRS+1)
      LEXAC=.FALSE.
      IPMACR=C_NULL_PTR
      IPSOU=C_NULL_PTR
      REBFLG=.FALSE.
      CALL DOORFV(CDOOR,IPSYS,NPSYS,IPTRK,IFTRAK,IMPX,NABS,NBMIX,
     1 IDIR,NREG,NUN,IPHASE,LEXAC,MAT,VOL,KEYFLX,TITR,SUN,FUN,IPMACR,
     2 IPSOU,REBFLG)
*----
*  HOMOGENIZE THE MULTIBAND FLUX.
*----
      DO 100 IM=1,MI
      CALL XDRSET(PAV,(NBNRS+1)**2,0.0)
      DO 70 JNBN=0,NBNRS
      DO 60 I=1,NREG
      IBM=MAT(I)
      IF(IBM.EQ.0) GO TO 60
      IOF=(IM-1)*NUN*(NBNRS+1)+JNBN*NUN+KEYFLX(I)-1
      PAV(IREX(IBM),JNBN)=PAV(IREX(IBM),JNBN)+FUN(IOF+1)*VOL(I)
   60 CONTINUE
   70 CONTINUE
      DO 90 I=0,NBNRS
      DO 80 J=0,NBNRS
      IF(VOLMER(I).NE.0.0) PAV(I,J)=PAV(I,J)*VOLMER(J)/VOLMER(I)
   80 CONTINUE
   90 CONTINUE
      KPSYS=LCMGIL(IPSYS,IASM+IM)
      CALL LCMPUT(KPSYS,'DRAGON-PAV',(NBNRS+1)**2,2,PAV(0,0))
  100 CONTINUE
      DEALLOCATE(SUN,FUN,NPSYS)
*----
*  RESPONSE MATRIX APPROACH. LOOP OVER THE SECONDARY SUBGROUPS.
*----
      MATRIX(:NBNRS*MI,:NBNRS*MI+1)=0.0D0
      DO 200 IM=1,MI
      KPSYS=LCMGIL(IPSYS,IASM+IM)
      CALL LCMGET(KPSYS,'DRAGON-PAV',PAV(0,0))
*----
*  LOOP OVER THE PRIMARY SUBGROUPS. MI+1 IS THE SOURCE.
*----
      DO 190 JM=1,MI+1
      IF(JM.LE.MI) THEN
         JNBMAX=NBNRS
      ELSE
         JNBMAX=1
      ENDIF
      DO 180 JNBN=1,JNBMAX
      CALL XDRSET(AWPHI(1),NBNRS,0.0)
      DO 160 I=1,NREG
      IBM=MAT(I)
      IF(IBM.EQ.0) GO TO 160
      JND=IREX(IBM)
      QQQ=0.0D0
      IF(JM.EQ.MI+1) THEN
         QQQ=QQQ+SIGGAR(IBM,0,IG,3)*GAMMA_V%VECTOR(IM)
         IF(JND.NE.0) THEN
           DO 130 JRES=1,NIRES
           DENSIT=CONR(JND,JRES)
           DO 120 JG=IG-NJJ(IG)+1,IG-1
           IF(GOLD(IRES,JG).NE.-1001.) CYCLE
           DO 110 KM=1,MRANK(JG)
           QQQ=QQQ+DENSIT*SCAT_M(JG,JRES)%MATRIX(IM,KM)*
     1     XFLUX(JND,KM,JG)
  110      CONTINUE
  120      CONTINUE
  130      CONTINUE
         ENDIF
      ELSE IF((JND.EQ.JNBN).AND.(JM.NE.IM)) THEN
         DO 140 JRES=1,NIRES
         DENSIT=CONR(JND,JRES)
         IF(GOLD(IRES,IG).NE.-1001.) CYCLE
         IF(JM.EQ.IM) CYCLE
         QQQ=QQQ-DENSIT*SCAT_M(IG,JRES)%MATRIX(IM,JM)
  140    CONTINUE
      ENDIF
      DO 150 IND=1,NBNRS
      AWPHI(IND)=AWPHI(IND)+PAV(IND,JND)*REAL(QQQ,4)*VOL(I)/VOLMER(JND)
  150 CONTINUE
  160 CONTINUE
      DO 170 IND=1,NBNRS
      MATRIX(INM(IND,IM,NBNRS),INM(JNBN,JM,NBNRS))=AWPHI(IND)
  170 CONTINUE
  180 CONTINUE
  190 CONTINUE
  200 CONTINUE
*
      DO 210 I=1,NBNRS*MI
      MATRIX(I,I)=MATRIX(I,I)+1.0D0
  210 CONTINUE
      CALL ALSBD(NBNRS*MI,1,MATRIX,IER,NBNRS*MI)
      IF(IER.NE.0) CALL XABORT('USSEXD: SINGULAR MATRIX.')
      CALL XDRSET(XFLUX(1,1,IG),NBNRS*MAXNOR,0.0)
      DO 230 IND=1,NBNRS
      DO 220 IM=1,MI
      I1=INM(IND,IM,NBNRS)
      XFLUX(IND,IM,IG)=REAL(MATRIX(I1,NBNRS*MI+1))
  220 CONTINUE
  230 CONTINUE
* END OF RESPONSE MATRIX APPROACH.
*----
*  COMPUTE THE AVERAGED SOURCE.
*----
      ALLOCATE(FUN(NUN*MI),SUN(NUN*MI))
      CALL XDRSET(SUN,NUN*MI,0.0)
      ALLOCATE(NPSYS(MI))
      DO 250 IM=1,MI
      NPSYS(IM)=IASM+IM
      KPSYS=LCMGIL(IPSYS,IASM+IM)
      CALL LCMLEN(KPSYS,'FUNKNO$USS',ILENG,ITYLCM)
      IF(ILENG.EQ.NUN) THEN
        CALL LCMGET(KPSYS,'FUNKNO$USS',FUN((IM-1)*NUN+1))
      ELSE
        CALL XDRSET(FUN((IM-1)*NUN+1),NUN,0.0)
      ENDIF
      IOF=(IM-1)*NUN
      DO 240 I=1,NREG
      IBM=MAT(I)
      IF(IBM.EQ.0) GO TO 240
      IND=IREX(IBM)
      QQQ=0.0D0
      DO JRES=0,NIRES
        IF(JRES.EQ.0) THEN
          QQQ=QQQ+SIGGAR(IBM,0,IG,3)*GAMMA_V%VECTOR(IM)
        ELSE IF(IND.GT.0) THEN
          DENSIT=CONR(IND,JRES)
          DO JG=IG-NJJ(IG)+1,IG
            IF(GOLD(IRES,JG).NE.-1001.) CYCLE
            DO JM=1,MRANK(JG)
              IF((JG.EQ.IG).AND.(JM.EQ.IM)) CYCLE
              QQQ=QQQ+DENSIT*SCAT_M(JG,JRES)%MATRIX(IM,JM)*
     1        XFLUX(IND,JM,JG)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
      SUN(IOF+KEYFLX(I))=REAL(QQQ,4)
  240 CONTINUE
  250 CONTINUE
*
      IF(IMPX.GT.0) THEN
         WRITE(TEXT12,'(3A4)') (IPPT2(IRES,I),I=2,4)
         WRITE(6,'(15H USSEXD: GROUP=,I5,24H. SUBGROUP CALCULATION B,
     1   37HASED ON RESPONSE MATRICES.  ISOTOPE='',A12,2H''.)') IG,
     2   TEXT12
      ENDIF
*----
*  SOLVE FOR THE MULTIBAND FLUX (VECTOR OF LENGTH NREG).
*----
      IPMACR=C_NULL_PTR
      IPSOU=C_NULL_PTR
      REBFLG=.FALSE.
      CALL DOORFV(CDOOR,IPSYS,NPSYS,IPTRK,IFTRAK,IMPX,MI,NBMIX,IDIR,
     1 NREG,NUN,IPHASE,LEXAC,MAT,VOL,KEYFLX,TITR,SUN,FUN,IPMACR,IPSOU,
     2 REBFLG)
      DEALLOCATE(NPSYS)
*----
*  INTEGRATE THE REGION-ORDERED FLUX OVER SUBGROUPS AND COMPUTE UNGAR,
*  THE REGION-ORDERED FLUX.
*----
      CALL XDRSET(UNGAR(1,IRES,IG),NREG,0.0)
      DO 270 IM=1,MI
      KPSYS=LCMGIL(IPSYS,IASM+IM)
      IOF=(IM-1)*NUN
      CALL LCMPUT(KPSYS,'FUNKNO$USS',NUN,2,FUN(IOF+1))
*
      DO 260 I=1,NREG
      IOF=(IM-1)*NUN+KEYFLX(I)
      UNGAR(I,IRES,IG)=UNGAR(I,IRES,IG)+REAL(XI_V%VECTOR(IM)*FUN(IOF)/
     1 DELTAU(IG),4)
  260 CONTINUE
  270 CONTINUE
      DEALLOCATE(SUN,FUN)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(SCAT_M,JPLIB2)
      DEALLOCATE(MATRIX)
      DEALLOCATE(AWPHI,PAV)
      RETURN
      END
