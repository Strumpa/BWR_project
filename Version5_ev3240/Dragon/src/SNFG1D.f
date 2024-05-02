*DECK SNFG1D
      SUBROUTINE SNFG1D(LX,NMAT,IELEM,NLF,NSCT,U,W,PL,MAT,VOL,TOTAL,
     1 NCODE,ZCODE,QEXT,LFIXUP,LSHOOT,ISBS,NBS,ISBSM,BS,WX,CST,
     2 ISADPTX,NUN,FUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform one inner iteration for solving SN equations in 1D slab
* geometry. Albedo boundary conditions.
*
*Copyright:
* Copyright (C) 2020 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert, A. A. Calloo and C. Bienvenue
*
*Parameters: input
* LX      number of regions.
* NMAT    number of material mixtures.
* IELEM   measure of order of the spatial approximation polynomial:
*         =1 constant - default for HODD;
*         =2 linear - default for DG;
*         >3 higher orders.
* NLF     number of $\\mu$ levels.
* NSCT    number of Legendre components in the flux:
*         =1: isotropic sources;
*         =2: linearly anisotropic sources.
* U       base points in $\\mu$ of the SN quadrature.
* W       weights of the SN quadrature.
* PL      discrete values of the Legendre polynomials corresponding
*         to the SN quadrature.
* MAT     material mixture index in each region.
* VOL     volumes of each region.
* TOTAL   macroscopic total cross sections.
* NCODE   boundary condition indices.
* ZCODE   albedos.
* QEXT    Legendre components of the fixed source.
* LFIXUP  flag to enable negative flux fixup.
* LSHOOT  flag to enable/disable shooting method.
* ISBS    flag to indicate the presence or not of boundary fixed
*         sources.
* NBS     number of boundary fixed sources.
* ISBSM   flag array to indicate the presence or not of boundary fixed
*         source in each unit surface.
* BS      boundary source array with their intensities.
* WX      spatial closure relation weighting factors.
* CST     constants for the polynomial approximations.
* ISADPTX flag to enable/disable spatial adaptive flux calculations.
* NUN     total number of unknowns in vector FUNKNO   
*
*Parameters: input/output
* FUNKNO  Legendre components of the flux and boundary fluxes.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER LX,NMAT,IELEM,NLF,NSCT,MAT(LX),NCODE(2),ISBS,NBS,
     1 ISBSM(2*ISBS,NLF*ISBS),NUN
      LOGICAL LFIXUP,LSHOOT,ISADPTX
      REAL U(NLF),W(NLF),PL(NSCT,NLF),VOL(LX),TOTAL(0:NMAT),ZCODE(2),
     1 QEXT(IELEM,NSCT,LX),FUNKNO(NUN),BS(NBS*ISBS),WX(IELEM+1),
     2 CST(IELEM)
*----
*  LOCAL VARIABLES
*----
      REAL WX0(IELEM+1)
      DOUBLE PRECISION XNI,XNI1,XNI2,XNIA,XNIB,XNIA1,XNIA2,XNIB1,XNIB2
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: Q
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: Q2
      PARAMETER(RLOG=1.0E-8)
      LOGICAL ISSHOOT,ISFIX
*----
*  ALLOCATABLE ARRAYS
*----
      ALLOCATE(Q(IELEM),Q2(IELEM,IELEM+1))
*----
*  LENGTH OF FUNKNO COMPONENTS (IN ORDER)
*----
      LFLX=IELEM*LX*NSCT
      LXNI=NLF
*----
*  INNER ITERATION
*----

      CALL XDRSET(FUNKNO(1:LFLX),LFLX,0.0)
      XNI=0.0D0
      XNI1=0.0D0
      XNI2=0.0D0
      XNIA=0.0D0
      XNIA2=0.0D0
      WX0=WX

      ! SHOOTING METHOD (ONLY IF THERE IS A NON-VACUUM RIGHT
      ! BOUNDARY CONDITION.
      ISSHOOT=(ZCODE(2).NE.0.0).AND.LSHOOT
      IF(ISSHOOT) THEN
        NS=6
      ELSE
        NS=2
      ENDIF

      ! LOOP OVER ALL DIRECTIONS
      DO 200 M0=1,NLF/2

      ! LOOP FOR SHOOTING METHOD
      DO 500 IS=1,NS

      ! CHOOSE DIRECTION
      IF(MOD(IS,2).EQ.0) THEN
        M=NLF-M0+1 ! FORWARD
      ELSE
        M=M0       ! BACKWARD
      ENDIF

      ! SHOOTING METHOD BOUNDARY CONDITIONS.
      IF(ISSHOOT) THEN
        ! 1ST BACKWARD SWEEP
        IF(IS.EQ.1) THEN
          XNI=0.0D0
          XNI1=0.0D0
          XNI2=0.0D0
        ! 1ST FORWARD SWEEP
        ELSEIF(IS.EQ.2) THEN
          XNIA1=0.0D0
          IF(NCODE(1).EQ.4) THEN
          XNIA1=REAL(XNI)
          XNI=0.0D0
          ELSE
          XNI=ZCODE(1)*REAL(XNI)
          ENDIF
        ! 2ND BACKWARD SWEEP
        ELSEIF(IS.EQ.3) THEN
          XNIA2=0.0D0
          XNIA=0.0D0
          IF(NCODE(1).EQ.4) THEN
          XNIA2=REAL(XNI)
          ELSE
          XNIA=REAL(XNI)
          ENDIF
          XNI=1.0D0
        ! 2ND FORWARD SWEEP
        ELSEIF(IS.EQ.4) THEN
          IF(NCODE(1).EQ.4) THEN
          XNIB1=REAL(XNI)
          XNI1=XNIA1/(1.0D0+XNIA1-XNIB1)
          XNI=1.0D0
          ELSE
          XNI=ZCODE(1)*REAL(XNI)
          ENDIF
        ! 3RD BACKWARD SWEEP
        ELSEIF(IS.EQ.5) THEN
          IF(NCODE(1).EQ.4) THEN
          XNIB2=REAL(XNI)
          XNI2=XNIA2/(1.0D0+XNIA2-XNIB2)
          XNI=XNI1
          ELSE
          XNIB=REAL(XNI)
          XNI=ZCODE(2)*XNIA/(1.0D0+ZCODE(2)*(XNIA-XNIB))
          ENDIF
        ! 3RD FORWARD SWEEP
        ELSEIF(IS.EQ.6) THEN
          XNI=ZCODE(1)*XNI
          IF(NCODE(1).EQ.4) XNI=XNI2
          ENDIF
      ! NO SHOOTING METHOD BOUNDARY CONDITIONS
      ELSE
        IF(.NOT.LSHOOT) THEN
        IF(U(M).GT.0.0) THEN
          IF(NCODE(1).NE.4) FUNKNO(LFLX+M)=FUNKNO(LFLX+NLF-M+1)
        ELSE
          IF(NCODE(2).NE.4) FUNKNO(LFLX+M)=FUNKNO(LFLX+NLF-M+1)
        ENDIF
        XNI=0.0D0
        ELSE
        IF(IS.EQ.1) THEN
          XNI=0.0D0
        ELSE
          XNI=ZCODE(1)*XNI
        ENDIF
        ENDIF
      ENDIF

      ! BOUNDARY FIXED SOURCES
      IF(U(M).GT.0.0) THEN
        IF(ISBS.EQ.1.AND.ISBSM(1,M).NE.0) XNI=XNI+BS(ISBSM(1,M))
      ELSE
        IF(ISBS.EQ.1.AND.ISBSM(2,M).NE.0) XNI=XNI+BS(ISBSM(2,M))
      ENDIF
      
      ! SWEEPING OVER ALL VOXELS
      DO 30 I0=1,LX
      I=I0
      IF(U(M).LT.0.0) I=LX+1-I0

      ! X-BOUNDARIES CONDITIONS (NO SHOOTING)
      IF(.NOT.LSHOOT.AND.I0.EQ.1) THEN
        IF(U(M).GT.0.0) THEN
          XNI=FUNKNO(LFLX+M)*ZCODE(1)
        ELSE
          XNI=FUNKNO(LFLX+M)*ZCODE(2)
        ENDIF
      ENDIF
      
      ! DATA
      IBM=MAT(I)
      SIGMA=TOTAL(IBM)

      ! SOURCE DENSITY TERM
      DO IEL=1,IELEM
      Q(IEL)=0.0
      DO L=1,NSCT
      Q(IEL)=Q(IEL)+QEXT(IEL,L,I)*PL(L,M)/2.0
      ENDDO
      ENDDO

      ISFIX=.FALSE.
      DO WHILE (.NOT.ISFIX) ! LOOP FOR ADAPTIVE CALCULATION

      ! FLUX MOMENTS CALCULATIONS
      CALL XDDSET(Q2,IELEM*(IELEM+1),0.0D0)
      DO II=1,IELEM
      DO JJ=1,IELEM

        ! MOMENT COEFFICIENTS
        IF(II.EQ.JJ) THEN
          Q2(II,JJ)=SIGMA*VOL(I)+CST(II)**2*WX(JJ+1)*ABS(U(M))
        ELSEIF(II.LT.JJ) THEN
          IF(MOD(II+JJ,2).EQ.1) THEN
            Q2(II,JJ)=CST(II)*CST(JJ)*WX(JJ+1)*U(M)
          ELSE
            Q2(II,JJ)=CST(II)*CST(JJ)*WX(JJ+1)*ABS(U(M))
          ENDIF
        ELSE
          IF(MOD(II+JJ,2).EQ.1) THEN
            Q2(II,JJ)=CST(II)*CST(JJ)*(WX(JJ+1)-2.0D0)*U(M)
          ELSE
            Q2(II,JJ)=CST(II)*CST(JJ)*WX(JJ+1)*ABS(U(M))
          ENDIF
        ENDIF 
      ENDDO
      ENDDO

      ! SOURCE TERMS
      DO II=1,IELEM
        IF(MOD(II,2).EQ.1) THEN
          Q2(II,IELEM+1)=Q(II)*VOL(I)+CST(II)*(1-WX(1))*ABS(U(M))*XNI
        ELSE
          Q2(II,IELEM+1)=Q(II)*VOL(I)-CST(II)*(1+WX(1))*U(M)*XNI
        ENDIF
      ENDDO
      !Q2(1,1)=SIGMA*VOL(I)+2.0D0*ABS(U(M))
      !Q2(1,2)=Q(1)*VOL(I)+2.0D0*ABS(U(M))*XNI
      
      CALL ALSBD(IELEM,1,Q2,IER,IELEM)
      IF(IER.NE.0) CALL XABORT('SNFE1D: SINGULAR MATRIX.')

      ! ADAPTIVE CORRECTION OF WEIGHTING PARAMETERS
      IF(ISADPTX) THEN
         CALL SNADPT(IELEM,IELEM,1,Q2(:IELEM,IELEM+1),XNI,
     1   1.0,WX,ISFIX)
      ELSE
        ISFIX=.TRUE. 
      ENDIF

      END DO ! END OF ADAPTIVE LOOP

      ! CLOSURE RELATIONS
      IF(IELEM.EQ.1.AND.LFIXUP.AND.(Q2(1,2).LE.RLOG)) Q2(1,2)=0.0D0
      XNI=WX(1)*XNI
      DO II=1,IELEM
        IF(MOD(II,2).EQ.1) THEN
          XNI=XNI+CST(II)*WX(II+1)*Q2(II,IELEM+1)
        ELSE
          XNI=XNI+CST(II)*WX(II+1)*Q2(II,IELEM+1)*SIGN(1.0,U(M))
        ENDIF
      ENDDO
      IF(IELEM.EQ.1.AND.LFIXUP.AND.(XNI.LE.RLOG)) XNI=0.0D0
      WX=WX0
       
      IF(ISSHOOT.AND.IS.LT.5) GO TO 30

      !IF(LFIXUP.AND.(Q2(1,2).LE.RLOG)) Q2(1,2)=0.0
      !XNI=2.0D0*Q2(1,2)-XNI

      ! SAVE LEGENDRE MOMENT OF THE FLUX
      DO L=1,NSCT
      DO IEL=1,IELEM
      IOF=(I-1)*NSCT*IELEM+(L-1)*IELEM+IEL
      FUNKNO(IOF)=FUNKNO(IOF)+W(M)*REAL(Q2(IEL,IELEM+1))*PL(L,M)
      ENDDO
      ENDDO

   30 CONTINUE ! END OF X-LOOP

      ! SAVE BOUNDARIES FLUX
      IF(.NOT.LSHOOT) FUNKNO(LFLX+M)=REAL(XNI)

  500 CONTINUE ! END OF SHOOTING METHOD LOOP 
  200 CONTINUE ! END OF DIRECTION LOOP

      DEALLOCATE(Q,Q2)
      RETURN
      END
