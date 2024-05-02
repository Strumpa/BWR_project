*DECK SNFP13
      SUBROUTINE SNFE3D(NUN,NGEFF,IMPX,INCONV,NGIND,LX,LY,LZ,
     1 IELEM,EELEM,NM,NME,NMX,NMY,NMZ,NMAT,NPQ,NSCT,MAT,VOL,TOTAL,
     2 ESTOPW,NCODE,ZCODE,DELTAE,QEXT,LFIXUP,DU,DE,DZ,W,MRMX,MRMY,MRMZ,
     3 DC,DB,DA,PL,FUNKNO,ISLG,FLUXC,ISBS,NBS,ISBSM,BS,MAXL,WX,WY,WZ,WE,
     4 CST,ISADPT,IBFP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform one inner iteration for solving SN equations in 3D Cartesian
* geometry for the HODD method. Energy-angle multithreading. Albedo
* boundary conditions. Boltzmann-Fokker-Planck (BFP) discretization.
*
*Copyright:
* Copyright (C) 2021 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert, A. A. Calloo and C. Bienvenue
*
*Parameters: input
* NUN     total number of unknowns in vector FUNKNO.
* NGEFF   number of energy groups processed in parallel.
* IMPX    print flag (equal to zero for no print).
* INCONV  energy group convergence flag (set to .FALSE. if converged).
* NGIND   energy group indices assign to the NGEFF set.
* LX      number of meshes along X axis.
* LY      number of meshes along Y axis.
* LZ      number of meshes along Z axis.
* IELEM   measure of order of the spatial approximation polynomial:
*         =1 constant - default for HODD;
*         =2 linear - default for DG;
*         >3 higher orders.
* EELEM   measure of order of the energy approximation polynomial:
*         =1 constant - default for HODD;
*         =2 linear - default for DG;
*         >3 higher orders.
* NM      number of moments in space and energy for flux components
* NME     number of moments for energy boundaries components
* NMX     number of moments for X axis boundaries components
* NMY     number of moments for Y axis boundaries components
* NMZ     number of moments for Z axis boundaries components
* NMAT    number of material mixtures.
* NPQ     number of SN directions in height octants.
* NSCT    maximum number of spherical harmonics moments of the flux.
* MAT     material mixture index in each region.
* VOL     volumes of each region.
* TOTAL   macroscopic total cross sections.
* ESTOPW  stopping power.
* NCODE   boundary condition indices.
* ZCODE   albedos.
* DELTAE  energy group width in MeV.
* QEXT    Legendre components of the fixed source.
* LFIXUP  flag to enable negative flux fixup.
* DU      first direction cosines ($\\mu$).
* DE      second direction cosines ($\\eta$).
* DZ      third direction cosines ($\\xi$).
* W       weights.
* MRMX    quadrature index.
* MRMY    quadrature index.
* MRMZ    quadrature index.
* DC      diamond-scheme parameter.
* DB      diamond-scheme parameter.
* DA      diamond-scheme parameter.
* PL      discrete values of the spherical harmonics corresponding
*         to the 3D SN quadrature.
* ISBS    flag to indicate the presence or not of boundary fixed
*         sources.
* NBS     number of boundary fixed sources.
* ISBSM   flag array to indicate the presence or not of boundary fixed
*         source in each unit surface.
* BS      boundary source array with their intensities.
* MAXL    maximum size of boundary source array.
* WX      spatial X axis closure relation weighting factors.
* WY      spatial Y axis closure relation weighting factors.
* WZ      spatial Z axis closure relation weighting factors.
* WE      energy closure relation weighting factors.
* CST     constants for the polynomial approximations.
* ISADPTX flag to enable/disable adaptive X axis flux calculations.
* ISADPTY flag to enable/disable adaptive Y axis flux calculations.
* ISADPTZ flag to enable/disable adaptive Z axis flux calculations.
* ISADPTE flag to enable/disable adaptive energy flux calculations.
* IBFP    type of energy proparation relation:
*         =1 Galerkin type;
*         =2 heuristic Przybylski and Ligou type.
*
*Parameters: input/output
* FUNKNO  Legendre components of the flux and boundary fluxes.
* FLUXC   flux at the cutoff energy.
*
*-----------------------------------------------------------------------
*
#if defined(_OPENMP)
      USE omp_lib
#endif
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NUN,NGEFF,IMPX,NGIND(NGEFF),LX,LY,LZ,IELEM,EELEM,NM,NME,
     1 NMX,NMY,NMZ,NMAT,NPQ,NSCT,MAT(LX,LY,LZ),NCODE(6),MRMX(NPQ),
     2 MRMY(NPQ),MRMZ(NPQ),ISLG(NGEFF),ISBS,NBS,
     3 ISBSM(6*ISBS,NPQ*ISBS,NGEFF*ISBS),MAXL
      LOGICAL INCONV(NGEFF)
      REAL VOL(LX,LY,LZ),TOTAL(0:NMAT,NGEFF),ESTOPW(0:NMAT,2,NGEFF),
     1 ZCODE(6),DELTAE(NGEFF),QEXT(NUN,NGEFF),DU(NPQ),DE(NPQ),DZ(NPQ),
     2 W(NPQ),DC(LX,LY,NPQ),DB(LX,LZ,NPQ),DA(LY,LZ,NPQ),PL(NSCT,NPQ),
     3 FUNKNO(NUN,NGEFF),FLUXC(LX,LY,LZ),BS(MAXL*ISBS,NBS*ISBS),
     4 WX(IELEM+1),WY(IELEM+1),WZ(IELEM+1),WE(EELEM+1),
     5 CST(MAX(IELEM,EELEM))
      LOGICAL LFIXUP,ISADPT(4)
*----
*  LOCAL VARIABLES
*----
      INTEGER NPQD(8),IIND(8)
      PARAMETER(IUNOUT=6,RLOG=1.0E-8,PI=3.141592654)
      REAL BM,BP,WX0(IELEM+1),WY0(IELEM+1),WZ0(IELEM+1),WE0(EELEM+1)
      DOUBLE PRECISION V,Q(NM),Q2(NM,NM+1),FEP(NME),XNK(NMZ)
      LOGICAL ISFIX(4)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: INDANG
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: FLUX,FLUX0
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: FLUX_G,
     1 FLUX0_G
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: XNI
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: XNJ
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(INDANG(NPQ,8))
      ALLOCATE(XNI(NMX,LY,LZ),XNJ(NMY,LZ))
      ALLOCATE(FLUX(NM,NSCT,LX,LY,LZ),
     1 FLUX0(NME,NPQ,LX,LY,LZ))
      ALLOCATE(FLUX_G(NM,NSCT,LX,LY,LZ,NGEFF),
     1 FLUX0_G(NME,NPQ,LX,LY,LZ,NGEFF))
*----
*  LENGTH OF FUNKNO COMPONENTS (IN ORDER)
*----
      LFLX=NM*LX*LY*LZ*NSCT
      LXNI=NMX*LY*LZ*NPQ
      LXNJ=NMY*LX*LZ*NPQ
      LXNK=NMZ*LX*LY*NPQ
      LFEP=NME*LX*LY*LZ*NPQ
*----
*  SET OCTANT SWAPPING ORDER.
*----
      NPQD(:8)=0
      INDANG(:NPQ,:8)=0
      DO 10 M=1,NPQ
        VU=DU(M)
        VE=DE(M)
        VZ=DZ(M)
        IF((VU.GE.0.0).AND.(VE.GE.0.0).AND.(VZ.GE.0.0)) THEN
          IND=1
          JND=8
        ELSE IF((VU.LE.0.0).AND.(VE.GE.0.0).AND.(VZ.GE.0.0)) THEN
          IND=2
          JND=7
        ELSE IF((VU.LE.0.0).AND.(VE.LE.0.0).AND.(VZ.GE.0.0)) THEN
          IND=3
          JND=5
        ELSE IF((VU.GE.0.0).AND.(VE.LE.0.0).AND.(VZ.GE.0.0)) THEN
          IND=4
          JND=6
        ELSE IF((VU.GE.0.0).AND.(VE.GE.0.0).AND.(VZ.LE.0.0)) THEN
          IND=5
          JND=4
        ELSE IF((VU.LE.0.0).AND.(VE.GE.0.0).AND.(VZ.LE.0.0)) THEN
          IND=6
          JND=3
        ELSE IF((VU.LE.0.0).AND.(VE.LE.0.0).AND.(VZ.LE.0.0)) THEN
          IND=7
          JND=1
        ELSE
          IND=8
          JND=2
        ENDIF
        IIND(JND)=IND
        NPQD(IND)=NPQD(IND)+1
        INDANG(NPQD(IND),IND)=M
   10 CONTINUE
*----
*  MAIN LOOP OVER OCTANTS.
*----

      FLUX_G(:NM,:NSCT,:LX,:LY,:LZ,:NGEFF)=0.0D0
      FLUX0_G(:NME,:NPQ,:LX,:LY,:LZ,:NGEFF)=0.0D0

      DO 420 JND=1,8
      IND=IIND(JND)
*----
*  PRELIMINARY LOOPS FOR SETTING BOUNDARY CONDITIONS.
*----

*$OMP  PARALLEL DO
*$OMP+ PRIVATE(M,IG,VU,VE,VZ,M1,IOF,JOF,IEL,I,J,K,IPQD,E1)
*$OMP+ SHARED(FUNKNO) COLLAPSE(2)

      DO 150 IG=1,NGEFF
      DO 140 IPQD=1,NPQD(IND)
      IF(.NOT.INCONV(IG)) GO TO 140
      M=INDANG(IPQD,IND)
      VU=DU(M)
      VE=DE(M)
      VZ=DZ(M)
      ! X-BOUNDARY
      IF(VU.GT.0.0)THEN
         M1=MRMX(M)
         IF(NCODE(1).NE.4)THEN
         DO IEL=1,NMX
            DO J=1,LY
            DO K=1,LZ
             IOF=(((M-1)*LZ+(K-1))*LY+(J-1))*NMX+IEL
             JOF=(((M1-1)*LZ+(K-1))*LY+(J-1))*NMX+IEL
             FUNKNO(LFLX+IOF,IG)=FUNKNO(LFLX+JOF,IG)
            ENDDO
            ENDDO
         ENDDO
         ENDIF
      ELSEIF(VU.LT.0.0)THEN
         M1=MRMX(M)
         IF(NCODE(2).NE.4)THEN
         DO IEL=1,NMX
            DO J=1,LY
            DO K=1,LZ
             IOF=(((M-1)*LZ+(K-1))*LY+(J-1))*NMX+IEL
             JOF=(((M1-1)*LZ+(K-1))*LY+(J-1))*NMX+IEL
             FUNKNO(LFLX+IOF,IG)=FUNKNO(LFLX+JOF,IG)
            ENDDO
            ENDDO
         ENDDO
         ENDIF
      ENDIF
      ! Y-BOUNDARY
      IF(VE.GT.0.0)THEN
         M1=MRMY(M)
         IF(NCODE(3).NE.4)THEN
         DO IEL=1,NMY
            DO I=1,LX
            DO K=1,LZ
             IOF=(((M-1)*LZ+(K-1))*LX+(I-1))*NMY+IEL
             JOF=(((M1-1)*LZ+(K-1))*LX+(I-1))*NMY+IEL
             FUNKNO(LFLX+LXNI+IOF,IG)=FUNKNO(LFLX+LXNI+JOF,IG)
            ENDDO
            ENDDO
         ENDDO
         ENDIF
      ELSEIF(VE.LT.0.0)THEN
         M1=MRMY(M)
         IF(NCODE(4).NE.4)THEN
         DO IEL=1,NMY
            DO I=1,LX
            DO K=1,LZ
             IOF=(((M-1)*LZ+(K-1))*LX+(I-1))*NMY+IEL
             JOF=(((M1-1)*LZ+(K-1))*LX+(I-1))*NMY+IEL
             FUNKNO(LFLX+LXNI+IOF,IG)=FUNKNO(LFLX+LXNI+JOF,IG)
            ENDDO
            ENDDO
         ENDDO
         ENDIF
      ENDIF
      ! Z-BOUNDARY
      IF(VZ.GT.0.0)THEN
         M1=MRMZ(M)
         IF(NCODE(5).NE.4)THEN
         DO IEL=1,NMZ
            DO I=1,LX
            DO J=1,LY
             IOF=(((M-1)*LY+(J-1))*LX+(I-1))*NMZ+IEL
             JOF=(((M1-1)*LY+(J-1))*LX+(I-1))*NMZ+IEL
             E1=FUNKNO(LFLX+LXNI+LXNJ+IOF,IG)
             FUNKNO(LFLX+LXNI+LXNJ+IOF,IG)=FUNKNO(LFLX+LXNI+LXNJ+JOF,IG)
             FUNKNO(LFLX+LXNI+LXNJ+JOF,IG)=E1
            ENDDO
            ENDDO
         ENDDO
         ENDIF
      ELSEIF(VZ.LT.0.0)THEN
         M1=MRMZ(M)
         IF(NCODE(6).NE.4)THEN
         DO IEL=1,NMZ
            DO I=1,LX
            DO J=1,LY
             IOF=(((M-1)*LY+(J-1))*LX+(I-1))*NMZ+IEL
             JOF=(((M1-1)*LY+(J-1))*LX+(I-1))*NMZ+IEL
             E1=FUNKNO(LFLX+LXNI+LXNJ+IOF,IG)
             FUNKNO(LFLX+LXNI+LXNJ+IOF,IG)=FUNKNO(LFLX+LXNI+LXNJ+JOF,IG)
             FUNKNO(LFLX+LXNI+LXNJ+JOF,IG)=E1
            ENDDO
            ENDDO
         ENDDO
         ENDIF
      ENDIF
  140 CONTINUE
  150 CONTINUE

*$OMP END PARALLEL DO

*----
*  MAIN SWAPPING LOOPS FOR SN FLUX CALCULATION
*----

*$OMP  PARALLEL DO
*$OMP+ PRIVATE(ITID,FLUX,M,IG,XNI,XNJ,XNK,Q,Q2,IOF,IER,II,JJ,IEL,I,J,K)
*$OMP+ PRIVATE(FEP,FLUX0,BM,BP,IIE,IIX,IIY,IIZ,IX,IY,IZ)
*$OMP+ PRIVATE(IE,ISFIX,JX,JY,JZ,JE,TB,V,SIGMA,IBM,I0,J0,K0,L,IPQD)
*$OMP+ FIRSTPRIVATE(WE,WX,WY,WZ,WE0,WX0,WY0,WZ0) SHARED(FUNKNO)
*$OMP+ REDUCTION(+:FLUX_G,FLUX0_G,FLUXC) COLLAPSE(2)
      
      ! LOOP FOR GROUPS TO EXECUTE IN PARALLEL
      DO 410 IG=1,NGEFF

      ! LOOP OVER ALL DIRECTIONS
      DO 400 IPQD=1,NPQD(IND)
      IF(.NOT.INCONV(IG)) GO TO 400
      M=INDANG(IPQD,IND)
      IF(W(M).EQ.0.0) GO TO 400

      ! GET AND PRINT THREAD NUMBER 
#if defined(_OPENMP)
      ITID=omp_get_thread_num()
#else
      ITID=0
#endif
      IF(IMPX.GT.5) WRITE(IUNOUT,500) ITID,NGIND(IG),IPQD

      ! INITIALIZE FLUX
      FLUX(:NM,:NSCT,:LX,:LY,:LZ)=0.0D0
      FLUX0(:NME,:NPQ,:LX,:LY,:LZ)=0.0D0
      WE0=WE
      WX0=WX
      WY0=WY
      WZ0=WZ

*----
*  LOOP OVER X-, Y- AND Z-DIRECTED AXES.
*----

      ! X-AXIS LOOP
      DO 350 I0=1,LX
      I=I0
      IF((IND.EQ.2).OR.(IND.EQ.3).OR.(IND.EQ.6).OR.(IND.EQ.7)) I=LX+1-I

      ! Y-AXIS LOOP
      DO 310 J0=1,LY
      J=J0
      IF((IND.EQ.3).OR.(IND.EQ.4).OR.(IND.EQ.7).OR.(IND.EQ.8)) J=LY+1-J

      ! Z-BOUNDARIES CONDITIONS
      DO IEL=1,NMZ
      IOF=(((M-1)*LY+(J-1))*LX+(I-1))*NMZ+IEL
      IF((IND.EQ.1).OR.(IND.EQ.2).OR.(IND.EQ.3).OR.(IND.EQ.4)) THEN
        XNK(IEL)=FUNKNO(LFLX+LXNI+LXNJ+IOF,IG)*ZCODE(5)
      ELSE
        XNK(IEL)=FUNKNO(LFLX+LXNI+LXNJ+IOF,IG)*ZCODE(6)
      ENDIF
      ENDDO

      ! Z-BOUNDARIES FIXED SOURCES
      IF(ISBS.EQ.1) THEN
        IF(((IND.EQ.5).OR.(IND.EQ.6).OR.(IND.EQ.7).OR.(IND.EQ.8))
     1  .AND.ISBSM(6,M,IG).NE.0) THEN
          XNK(1)=XNK(1)+BS((I-1)*LY+J,ISBSM(6,M,IG)) 
        ELSEIF(((IND.EQ.1).OR.(IND.EQ.2).OR.(IND.EQ.3).OR.(IND.EQ.4))
     1  .AND.ISBSM(5,M,IG).NE.0) THEN
          XNK(1)=XNK(1)+BS((I-1)*LY+J,ISBSM(5,M,IG)) 
        ENDIF
      ENDIF

      ! Z-AXIS LOOP
      DO 280 K0=1,LZ
      K=K0
      IF((IND.EQ.5).OR.(IND.EQ.6).OR.(IND.EQ.7).OR.(IND.EQ.8)) K=LZ+1-K

      ! Y-BOUNDARIES CONDITIONS
      IF(J0.EQ.1) THEN
       DO IEL=1,NMY
        IOF=(((M-1)*LZ+(K-1))*LX+(I-1))*NMY+IEL
        IF((IND.EQ.1).OR.(IND.EQ.2).OR.(IND.EQ.5).OR.(IND.EQ.6)) THEN
          XNJ(IEL,K)=FUNKNO(LFLX+LXNI+IOF,IG)*ZCODE(3)
        ELSE
          XNJ(IEL,K)=FUNKNO(LFLX+LXNI+IOF,IG)*ZCODE(4)
        ENDIF
        ENDDO

        !Y-BOUNDARIES FIXED SOURCES
        IF(ISBS.EQ.1) THEN
          IF(((IND.EQ.3).OR.(IND.EQ.4).OR.(IND.EQ.7).OR.(IND.EQ.8))
     1    .AND.ISBSM(4,M,IG).NE.0) THEN
            XNJ(1,K)=XNJ(1,K)+BS((I-1)*LZ+K,ISBSM(4,M,IG)) 
          ELSEIF(((IND.EQ.1).OR.(IND.EQ.2).OR.(IND.EQ.5).OR.(IND.EQ.6))
     1    .AND.ISBSM(3,M,IG).NE.0) THEN
            XNJ(1,K)=XNJ(1,K)+BS((I-1)*LZ+K,ISBSM(3,M,IG)) 
          ENDIF
        ENDIF
      ENDIF

      ! X-BOUNDARIES CONDITIONS
      IF(I0.EQ.1) THEN
        DO IEL=1,NMX
        IOF=(((M-1)*LZ+(K-1))*LY+(J-1))*NMX+IEL
        IF((IND.EQ.1).OR.(IND.EQ.4).OR.(IND.EQ.5).OR.(IND.EQ.8)) THEN
          XNI(IEL,J,K)=FUNKNO(LFLX+IOF,IG)*ZCODE(1)
        ELSE
          XNI(IEL,J,K)=FUNKNO(LFLX+IOF,IG)*ZCODE(2)
        ENDIF
        ENDDO

        ! X-BOUNDARIES FIXED SOURCES
        IF(ISBS.EQ.1) THEN
          IF(((IND.EQ.2).OR.(IND.EQ.3).OR.(IND.EQ.6).OR.(IND.EQ.7))
     1    .AND.ISBSM(2,M,IG).NE.0) THEN
            XNI(1,J,K)=XNI(1,J,K)+BS((J-1)*LZ+K,ISBSM(2,M,IG)) 
          ELSEIF(((IND.EQ.1).OR.(IND.EQ.4).OR.(IND.EQ.5).OR.(IND.EQ.8))
     1    .AND.ISBSM(1,M,IG).NE.0) THEN
            XNI(1,J,K)=XNI(1,J,K)+BS((J-1)*LZ+K,ISBSM(1,M,IG)) 
          ENDIF
        ENDIF
      ENDIF

      ! DATA
      IBM=MAT(I,J,K)
      IF(IBM.EQ.0) GO TO 280
      SIGMA=TOTAL(IBM,IG)
      BM=ESTOPW(IBM,1,IG)/DELTAE(IG)
      BP=ESTOPW(IBM,2,IG)/DELTAE(IG)
      V=VOL(I,J,K)

      ! TYPE OF ENERGY PROPAGATION FACTOR
      IF(IBFP.EQ.1) THEN ! GALERKIN TYPE
        TB=BM/BP
        WE(1)=WE(1)*TB
        WE(2:EELEM+1)=(WE(2:EELEM+1)-1)*TB+1
      ELSE ! PRZYBYLSKI AND LIGOU TYPE
        TB=1.0
      ENDIF
    
      ! SOURCE DENSITY TERM
      DO IEL=1,NM
      Q(IEL)=0.0D0
      DO L=1,NSCT
      IOF=((((K-1)*LY+(J-1))*LX+(I-1))*NSCT+(L-1))*NM+IEL
      Q(IEL)=Q(IEL)+QEXT(IOF,IG)*PL(L,M)/(4.0D0*PI)
      ENDDO
      ENDDO

      ! ENERGY GROUP UPPER BOUNDARY INCIDENT FLUX
      DO IEL=1,NME
      IOF=((((K-1)*LY+(J-1))*LX+(I-1))*NPQ+(M-1))*NME+IEL 
      FEP(IEL)=QEXT(LFLX+LXNI+LXNJ+LXNK+IOF,IG)
      ENDDO

      ISFIX=.FALSE.
      DO WHILE (.NOT.ALL(ISFIX)) ! LOOP FOR ADAPTIVE CALCULATION

      ! FLUX MOMENT COEFFICIENTS MATRIX
      CALL XDDSET(Q2,NM*(NM+1),0.0D0)

      DO IZ=1,IELEM
      DO JZ=1,IELEM
      DO IY=1,IELEM
      DO JY=1,IELEM
      DO IX=1,IELEM
      DO JX=1,IELEM
      DO IE=1,EELEM
      DO JE=1,EELEM 
        II=IELEM**2*EELEM*(IZ-1)+IELEM*EELEM*(IY-1)+EELEM*(IX-1)+IE
        JJ=IELEM**2*EELEM*(JZ-1)+IELEM*EELEM*(JY-1)+EELEM*(JX-1)+JE

        ! DIAGONAL TERMS
        IF(II.EQ.JJ) THEN
          Q2(II,JJ)=(SIGMA+CST(IE)**2*WE(JE+1)*BP+(IE-1)*(BM-BP))*V
     1              +CST(IX)**2*WX(JX+1)*ABS(DA(J,K,M))
     2              +CST(IY)**2*WY(JY+1)*ABS(DB(I,K,M))
     3              +CST(IZ)**2*WZ(JZ+1)*ABS(DC(I,J,M))

        ! UPPER DIAGONAL TERMS
        ELSEIF(II.LT.JJ) THEN
          IF(IZ.EQ.JZ) THEN
          IF(IY.EQ.JY) THEN
          IF(IX.EQ.JX) THEN
            ! ENERGY TERMS
            IF(MOD(IE+JE,2).EQ.1) THEN
              Q2(II,JJ)=-CST(IE)*CST(JE)*WE(JE+1)*BP*V
            ELSE
              Q2(II,JJ)=CST(IE)*CST(JE)*WE(JE+1)*BP*V
            ENDIF
          ELSEIF(IE.EQ.JE) THEN
            ! X-SPACE TERMS
            IF(MOD(IX+JX,2).EQ.1) THEN
              Q2(II,JJ)=CST(IX)*CST(JX)*WX(JX+1)*DA(J,K,M)
            ELSE
              Q2(II,JJ)=CST(IX)*CST(JX)*WX(JX+1)*ABS(DA(J,K,M))
            ENDIF
          ENDIF
          ELSEIF(IX.EQ.JX.AND.IE.EQ.JE) THEN
            ! Y-SPACE TERMS
            IF(MOD(IY+JY,2).EQ.1) THEN
              Q2(II,JJ)=CST(IY)*CST(JY)*WY(JY+1)*DB(I,K,M)
            ELSE
              Q2(II,JJ)=CST(IY)*CST(JY)*WY(JY+1)*ABS(DB(I,K,M))
            ENDIF
          ENDIF
          ELSEIF(IY.EQ.JY.AND.IX.EQ.JX.AND.IE.EQ.JE) THEN
            ! Z-SPACE TERMS
            IF(MOD(IZ+JZ,2).EQ.1) THEN
              Q2(II,JJ)=CST(IZ)*CST(JZ)*WZ(JZ+1)*DC(I,J,M)
            ELSE
              Q2(II,JJ)=CST(IZ)*CST(JZ)*WZ(JZ+1)*ABS(DC(I,J,M))
            ENDIF
          ENDIF
       
        ! UNDER DIAGONAL TERMS
        ELSE
          IF(IZ.EQ.JZ) THEN
          IF(IY.EQ.JY) THEN
          IF(IX.EQ.JX) THEN
            ! ENERGY TERMS
            IF(MOD(IE+JE,2).EQ.1) THEN
              Q2(II,JJ)=-CST(IE)*CST(JE)*(WE(JE+1)*BP-BM-BP)*V
            ELSE
              Q2(II,JJ)=CST(IE)*CST(JE)*(WE(JE+1)*BP+BM-BP)*V
            ENDIF
          ELSEIF(IE.EQ.JE) THEN
            ! X-SPACE TERMS
            IF(MOD(IX+JX,2).EQ.1) THEN
              Q2(II,JJ)=CST(IX)*CST(JX)*(WX(JX+1)-2)*DA(J,K,M)
            ELSE
              Q2(II,JJ)=CST(IX)*CST(JX)*WX(JX+1)*ABS(DA(J,K,M))
            ENDIF
          ENDIF
          ELSEIF(IX.EQ.JX.AND.IE.EQ.JE) THEN
            ! Y-SPACE TERMS
            IF(MOD(IY+JY,2).EQ.1) THEN
              Q2(II,JJ)=CST(IY)*CST(JY)*(WY(JY+1)-2)*DB(I,K,M)
            ELSE
              Q2(II,JJ)=CST(IY)*CST(JY)*WY(JY+1)*ABS(DB(I,K,M))
            ENDIF
          ENDIF
          ELSEIF(IY.EQ.JY.AND.IX.EQ.JX.AND.IE.EQ.JE) THEN
            ! Z-SPACE TERMS
            IF(MOD(IZ+JZ,2).EQ.1) THEN
              Q2(II,JJ)=CST(IZ)*CST(JZ)*(WZ(JZ+1)-2)*DC(I,J,M)
            ELSE
              Q2(II,JJ)=CST(IZ)*CST(JZ)*WZ(JZ+1)*ABS(DC(I,J,M))
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      ! FLUX SOURCE VECTOR
      DO IZ=1,IELEM
      DO IY=1,IELEM
      DO IX=1,IELEM
      DO IE=1,EELEM

        II=IELEM**2*EELEM*(IZ-1)+IELEM*EELEM*(IY-1)+EELEM*(IX-1)+IE
        IIE=IELEM**2*(IZ-1)+IELEM*(IY-1)+IX        
        IIX=IELEM*EELEM*(IZ-1)+EELEM*(IY-1)+IE
        IIY=IELEM*EELEM*(IZ-1)+EELEM*(IX-1)+IE
        IIZ=IELEM*EELEM*(IY-1)+EELEM*(IX-1)+IE
        Q2(II,NM+1)=Q(II)*V
        ! ENERGY TERMS
        IF(MOD(IE,2).EQ.1) THEN
          Q2(II,NM+1)=Q2(II,NM+1)+CST(IE)*(BM-WE(1)*BP)*FEP(IIE)*V
        ELSE
          Q2(II,NM+1)=Q2(II,NM+1)+CST(IE)*(BM+WE(1)*BP)*FEP(IIE)*V
        ENDIF
        ! X-SPACE TERMS
        IF(MOD(IX,2).EQ.1) THEN
          Q2(II,NM+1)=Q2(II,NM+1)+CST(IX)*(1-WX(1))
     1                *XNI(IIX,J,K)*ABS(DA(J,K,M))
        ELSE
          Q2(II,NM+1)=Q2(II,NM+1)-CST(IX)*(1+WX(1))
     1                *XNI(IIX,J,K)*DA(J,K,M)
        ENDIF
        ! Y-SPACE TERMS
        IF(MOD(IY,2).EQ.1) THEN
          Q2(II,NM+1)=Q2(II,NM+1)+CST(IY)*(1-WY(1))
     1                *XNJ(IIY,K)*ABS(DB(I,K,M))
        ELSE
          Q2(II,NM+1)=Q2(II,NM+1)-CST(IY)*(1+WY(1))
     1                *XNJ(IIY,K)*DB(I,K,M)
        ENDIF
        ! Z-SPACE TERMS
        IF(MOD(IZ,2).EQ.1) THEN
          Q2(II,NM+1)=Q2(II,NM+1)+CST(IZ)*(1-WZ(1))
     1                *XNK(IIZ)*ABS(DC(I,J,M))
        ELSE
          Q2(II,NM+1)=Q2(II,NM+1)-CST(IZ)*(1+WZ(1))
     1                *XNK(IIZ)*DC(I,J,M)
        ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      CALL ALSBD(NM,1,Q2,IER,NM)
      IF(IER.NE.0) CALL XABORT('SNFE2D: SINGULAR MATRIX.')

      ! ADAPTIVE CORRECTION OF WEIGHTING PARAMETERS
      IF(ANY(ISADPT)) THEN
        IF(ISADPT(1)) THEN
          CALL SNADPT(EELEM,NM,IELEM**3,Q2(1:EELEM:1,NM+1),
     1    FEP,TB,WE,ISFIX(1))
        ELSE
          ISFIX(1)=.TRUE.
        ENDIF
        IF(ISADPT(2)) THEN
          CALL SNADPT(IELEM,NM,EELEM*IELEM**2,
     1    Q2(1:IELEM*EELEM:IELEM,NM+1),XNI(:NMX,J,K),
     2    1.0,WX,ISFIX(2))
        ELSE
          ISFIX(2)=.TRUE.
        ENDIF
        IF(ISADPT(3)) THEN
          CALL SNADPT(IELEM,NM,EELEM*IELEM**2,
     1    Q2(1:EELEM*IELEM**2:IELEM*EELEM,NM+1),XNJ(:NMY,K),
     2    1.0,WY,ISFIX(3))
        ELSE
          ISFIX(3)=.TRUE.
        ENDIF
        IF(ISADPT(4)) THEN
          CALL SNADPT(IELEM,NM,EELEM*IELEM**2,
     1    Q2(1:NM:EELEM*IELEM**2,NM+1),XNK,1.0,WZ,ISFIX(4))
        ELSE
          ISFIX(4)=.TRUE.
        ENDIF
      ELSE
        ISFIX=.TRUE.
      ENDIF

      END DO ! END OF ADAPTIVE LOOP

      ! CLOSURE RELATIONS
      IF(IELEM.EQ.1.AND.LFIXUP.AND.(Q2(1,2).LE.RLOG)) Q2(1,2)=0.0
      XNI(:NMX,J,K)=WX(1)*XNI(:NMX,J,K)
      XNJ(:NMY,K)=WY(1)*XNJ(:NMY,K)
      XNK(:NMZ)=WZ(1)*XNK(:NMZ)
      FEP(:NME)=WE(1)*FEP(:NME)
      DO IZ=1,IELEM
      DO IY=1,IELEM
      DO IX=1,IELEM
      DO IE=1,EELEM
        
        II=IELEM**2*EELEM*(IZ-1)+IELEM*EELEM*(IY-1)+EELEM*(IX-1)+IE
        IIE=IELEM**2*(IZ-1)+IELEM*(IY-1)+IX        
        IIX=IELEM*EELEM*(IZ-1)+EELEM*(IY-1)+IE
        IIY=IELEM*EELEM*(IZ-1)+EELEM*(IX-1)+IE
        IIZ=IELEM*EELEM*(IY-1)+EELEM*(IX-1)+IE
        
        ! ENERGY  
        IF(MOD(IE,2).EQ.1) THEN
          FEP(IIE)=FEP(IIE)+CST(IE)*WE(IE+1)*Q2(II,NM+1)
        ELSE
          FEP(IIE)=FEP(IIE)-CST(IE)*WE(IE+1)*Q2(II,NM+1)
        ENDIF
        ! X-SPACE
        IF(MOD(IX,2).EQ.1) THEN
          XNI(IIX,J,K)=XNI(IIX,J,K)+CST(IX)*WX(IX+1)
     1                 *Q2(II,NM+1)
        ELSE
          XNI(IIX,J,K)=XNI(IIX,J,K)+CST(IX)*WX(IX+1)
     1                 *Q2(II,NM+1)*SIGN(1.0,DA(J,K,M))
        ENDIF
        ! Y-SPACE
        IF(MOD(IY,2).EQ.1) THEN
          XNJ(IIY,K)=XNJ(IIY,K)+CST(IY)*WY(IY+1)
     1                 *Q2(II,NM+1)
        ELSE
          XNJ(IIY,K)=XNJ(IIY,K)+CST(IY)*WY(IY+1)
     1                 *Q2(II,NM+1)*SIGN(1.0,DB(I,K,M))
        ENDIF
        ! Z-SPACE
        IF(MOD(IZ,2).EQ.1) THEN
          XNK(IIZ)=XNK(IIZ)+CST(IZ)*WZ(IZ+1)
     1                 *Q2(II,NM+1)
        ELSE
          XNK(IIZ)=XNK(IIZ)+CST(IZ)*WZ(IZ+1)
     1                 *Q2(II,NM+1)*SIGN(1.0,DC(I,J,M))
        ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      IF(IELEM.EQ.1.AND.LFIXUP.AND.(XNI(1,J,K).LE.RLOG))
     1 XNI(1,J,K)=0.0
      IF(IELEM.EQ.1.AND.LFIXUP.AND.(XNJ(1,K).LE.RLOG))
     1 XNJ(1,K)=0.0
      IF(IELEM.EQ.1.AND.LFIXUP.AND.(XNK(1).LE.RLOG)) XNK(1)=0.0
      WE=WE0
      WX=WX0
      WY=WY0
      WZ=WZ0

      ! SAVE ENERGY GROUP LOWER BOUNDARY OUTGOING FLUX
      FLUX0(:NME,M,I,J,K)=REAL(FEP(:NME))/DELTAE(IG)

      ! SAVE LAST GROUP LOWER BOUNDARY FLUX
      IF(ISLG(IG).EQ.1) THEN
      FLUXC(I,J,K)=FLUXC(I,J,K)+W(M)*REAL(FLUX0(1,M,I,J,K))*PL(1,M)
      ENDIF

      ! SAVE LEGENDRE MOMENT OF THE FLIX
      DO L=1,NSCT
      DO IEL=1,NM
      FLUX(IEL,L,I,J,K)=FLUX(IEL,L,I,J,K)+W(M)*Q2(IEL,NM+1)*
     1 PL(L,M)
      ENDDO
      ENDDO

  280 CONTINUE ! END OF Z-LOOP

      ! SAVE BOUNDARY CONDITIONS
      DO IEL=1,NMZ
      IOF=(((M-1)*LY+(J-1))*LX+(I-1))*NMZ+IEL
      FUNKNO(LFLX+LXNI+LXNJ+IOF,IG)=REAL(XNK(IEL))
      ENDDO

  310 CONTINUE ! END OF Y-LOOP

      ! SAVE BOUNDARY CONDITIONS
      DO K=1,LZ
      DO IEL=1,NMY
      IOF=(((M-1)*LZ+(K-1))*LX+(I-1))*NMY+IEL
      FUNKNO(LFLX+LXNI+IOF,IG)=REAL(XNJ(IEL,K))
      ENDDO
      ENDDO

  350 CONTINUE ! END OF X-LOOP

      ! SAVE BOUNDARY CONDITIONS
      DO K=1,LZ
      DO J=1,LY
      DO IEL=1,NMX
      IOF=(((M-1)*LZ+(K-1))*LY+(J-1))*NMX+IEL
      FUNKNO(LFLX+IOF,IG)=REAL(XNI(IEL,J,K))
      ENDDO
      ENDDO
      ENDDO

      ! SAVE FLUX INFORMATION
      FLUX_G(:,:,:,:,:,IG)=FLUX_G(:,:,:,:,:,IG)+FLUX(:,:,:,:,:)
      FLUX0_G(:,:,:,:,:,IG)=FLUX0_G(:,:,:,:,:,IG)+FLUX0(:,:,:,:,:)


  400 CONTINUE ! END OF DIRECTION LOOP
  410 CONTINUE ! END OF ENERGY LOOP
*$OMP END PARALLEL DO
  420 CONTINUE ! END OF OCTANT LOOP

      ! SAVE FLUX INFORMATION
      DO 430 IG=1,NGEFF
        IF(.NOT.INCONV(IG)) GO TO 430
        FUNKNO(:LFLX,IG)=
     1  RESHAPE(REAL(FLUX_G(:NM,:NSCT,:LX,:LY,:LZ,IG)), 
     2  (/LFLX/))
        FUNKNO(LFLX+LXNI+LXNJ+LXNK+1:LFLX+LXNI+LXNJ+LXNK+LFEP,IG)=
     1  RESHAPE(REAL(FLUX0_G(:NME,:NPQ,:LX,:LY,:LZ,IG)),
     2  (/ LFEP /) )
  430 CONTINUE

*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(FLUX0_G,FLUX_G,FLUX0,FLUX,XNJ,XNI,INDANG)
      RETURN
  500 FORMAT(16H SNFP13: thread=,I8,12H --->(group=,I4,7H angle=,I4,1H))
      END
