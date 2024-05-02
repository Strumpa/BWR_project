*DECK RESGEO
      SUBROUTINE RESGEO(IPMAP,IPMTX,LX,LY,LZ,NFUEL,IMPX,IGEO,NX,NY,NZ,
     1 NCH,NB,NTOT,LNAP,IPCPO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Create and check the fuel-map geometry.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* E. Varin, D. Sekki and V. Descotes
*
*Update(s):
* R. Chambon 2014
*
*Parameters: input
* IPMAP  pointer to fuel-map information.
* IPMTX  pointer to matex information.
* LX     number of elements along x-axis in geometry.
* LY     number of elements along y-axis in geometry.
* LZ     number of elements along z-axis in geometry.
* NFUEL  number of fuel types.
* IMPX   printing index (=0 for no print).
* IGEO   type of geometry (=7 or =9)
*
*Parameters: output
* NX     number of elements along x-axis in fuel map.
* NY     number of elements along y-axis in fuel map.
* NZ     number of elements along z-axis in fuel map.
* NCH    number of reactor channels.
* NB     number of fuel bundles per channel.
* NTOT   total number of fuel bundles.
* LNAP   Flag to call NAP: module to unfold geometry at assembly level
* IPCPO  pointer to multicompo information
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAP,IPMTX,IPCPO,IPGNW
      INTEGER LX,LY,LZ,NFUEL,IGEO,NX,NY,NZ,NCH,NB,NTOT
      LOGICAL LNAP
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40,IOUT=6,EPSI=1.0E-4)
      INTEGER ISTATE(NSTATE),JENT(1),IENT(1),JENT2(3),IENT2(3),NCODE(6),
     1 ICODE(6)
      TYPE(C_PTR) KENT(1),KENT2(3)
      REAL    GEOXX(LX+1),GEOYY(LY+1),GEOZZ(LZ+1),GEOSI,GMAPSI,ZCODE(6)
      CHARACTER HENT(1)*12,HENT2(3)*12,TEXT*12
      DOUBLE PRECISION DFLOT
*----
*  ALLOCATABLE STATEMENTS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ISPLX,ISPLY,ISPLZ,MAT
      REAL, ALLOCATABLE, DIMENSION(:) :: GMAPX,GMAPY,GMAPZ
*----
*  FUEL-MAP GEOMETRY
*----
      IF(IMPX.GT.1)WRITE(IOUT,*)'** CREATING FUEL-MAP GEOMETRY **'
      CALL LCMSIX(IPMAP,'GEOMAP',1)
      NENT=1
      JENT(1)=0
      HENT(1)='GEOMAP'
      IENT(1)=1
      KENT(1)=IPMAP
      CALL GEOD(NENT,HENT,IENT,JENT,KENT)
      IF(IMPX.GT.3)CALL LCMLIB(IPMAP)
*
      IF(LNAP) THEN
*----
*  FUEL-MAP GEOMETRY UNFOLDING WITH NAP:
*----
        IF(.NOT.C_ASSOCIATED(IPCPO)) THEN
          CALL XABORT('RESGEO: COMPO LCM OBJECT MISSING AT RHS.')
        ENDIF
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.3)CALL XABORT('@RESGEO: CHARACTER DATA EXPECTED.')
        IF(TEXT.NE.':::') CALL XABORT('@RESGEO: ::: keyword EXPECTED.')
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.3)CALL XABORT('@RESGEO: CHARACTER DATA EXPECTED.')
        IF(TEXT.NE.'NAP:') CALL XABORT('@RESGEO: NAP: keyword '
     1   //'EXPECTED.')
        CALL LCMOP(IPGNW,'GEONEW',0,1,0)
        CALL LCMSIX(IPMAP,' ',0)
        CALL LCMSIX(IPMAP,'GEOMAP',1)
        NENT2=3
        JENT2(1)=0
        JENT2(2)=2
        JENT2(3)=2
        HENT2(1)='GEONEW'
        HENT2(1)='GEOOLD'
        HENT2(1)='COMPO'
        IENT2(1)=1
        IENT2(2)=1
        IENT2(3)=1
        KENT2(1)=IPGNW
        KENT2(2)=IPMAP
        KENT2(3)=IPCPO
        CALL NAP(NENT2,HENT2,IENT2,JENT2,KENT2)
        CALL LCMSIX(IPMAP,' ',0)
        IF(IMPX.GT.3)CALL LCMLIB(IPMAP)
        CALL LCMDEL(IPMAP,'GEOMAP')
        IF(IMPX.GT.3)CALL LCMLIB(IPMAP)
        CALL LCMSIX(IPMAP,'GEOMAP',1)
        IF(IMPX.GT.3)CALL LCMLIB(IPMAP)
        CALL LCMEQU(IPGNW,IPMAP)
        IF(IMPX.GT.3)CALL LCMLIB(IPMAP)
        CALL LCMCL(IPGNW,1)
      ENDIF
****
      CALL LCMSIX(IPMAP,' ',0)
      IF(IMPX.GT.3)CALL LCMLIB(IPMAP)
      CALL LCMSIX(IPMAP,'GEOMAP',1)
      IF(IMPX.GT.3)CALL LCMLIB(IPMAP)
****
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPMAP,'STATE-VECTOR',ISTATE)
      IF(ISTATE(1).NE.IGEO) CALL XABORT('@RESGEO: THE GEOMETRY '
     1 // 'IN FUEL-MAP MUST HAVE THE SAME TYPE AS IN THE MATEX-OBJECT')
      IGEO=ISTATE(1)
      NX=ISTATE(3)
      NY=ISTATE(4)
      NZ=ISTATE(5)
*----
*  READ FUEL-MAP GEOMETRY AND PERFORM MESH-SPLITTING
*----
      IMPX0=MAX(0,IMPX-1)
      NX2=NX
      NY2=NY
      IF(IGEO.GE.8) NY2=1
      NZ2=NZ
      ALLOCATE(ISPLX(NX2),ISPLY(NY2),ISPLZ(NZ2))
      ISPLTL=0
      ISPLTH=0
      IHEX=0
      CALL LCMLEN(IPMAP,'SPLITL',ILEN,ITYLCM)
      IF(ILEN.GT.0) CALL LCMGET(IPMAP,'SPLITL',ISPLTL)
      CALL LCMLEN(IPMAP,'SPLITH',ILEN,ITYLCM)
      IF(ILEN.GT.0) CALL LCMGET(IPMAP,'SPLITH',ISPLTH)
      IF(IGEO.LT.8) THEN
        CALL LCMLEN(IPMAP,'SPLITX',ILONG,ITYLCM)
        IF(ILONG.GT.0) THEN
          CALL LCMGET(IPMAP,'SPLITX',ISPLX)
          NX2=0
          DO IOLD=1,NX
            NX2=NX2+ISPLX(IOLD)
          ENDDO
        ENDIF
        CALL LCMLEN(IPMAP,'SPLITY',ILONG,ITYLCM)
        IF(ILONG.GT.0) THEN
          CALL LCMGET(IPMAP,'SPLITY',ISPLY)
          NY2=0
          DO IOLD=1,NY
            NY2=NY2+ISPLY(IOLD)
          ENDDO
        ENDIF
      ELSEIF((ISPLTH.NE.0).AND.((IGEO.EQ.8).OR.(IGEO.EQ.9))) THEN
        NX2=NX*6*(ISPLTH**2)
        CALL LCMGET(IPMAP,'IHEX',IHEX)
      ELSEIF((ISPLTL.NE.0).AND.((IGEO.EQ.8).OR.(IGEO.EQ.9))) THEN
        NX2=NX*3*(ISPLTL**2)
        CALL LCMGET(IPMAP,'IHEX',IHEX)
      ENDIF
      CALL LCMLEN(IPMAP,'SPLITZ',ILONG,ITYLCM)
      IF(ILONG.GT.0) THEN
        CALL LCMGET(IPMAP,'SPLITZ',ISPLZ)
        NZ2=0
        DO IOLD=1,NZ
          NZ2=NZ2+ISPLZ(IOLD)
        ENDDO
      ENDIF
      MAXPTS=NX2*NY2*NZ2
      MAXX=NX2
      IF(IHEX.EQ.1) THEN
        MAXPTS=12*MAXPTS
        MAXX=12*MAXX
      ELSE IF((IHEX.EQ.2).OR.(IHEX.EQ.3)) THEN
        MAXPTS=6*MAXPTS
        MAXX=6*MAXX
      ELSE IF(IHEX.EQ.4) THEN
        MAXPTS=4*MAXPTS
        MAXX=4*MAXX
      ELSE IF(IHEX.EQ.5) THEN
        MAXPTS=3*MAXPTS
        MAXX=3*MAXX
      ELSE IF((IHEX.GE.6).AND.(IHEX.LE.8)) THEN
        MAXPTS=2*MAXPTS
        MAXX=2*MAXX
      ENDIF
      ALLOCATE(MAT(MAXPTS),GMAPX(MAXX+1),GMAPY(NY2+1),GMAPZ(NZ2+1))
      CALL READ3D(MAXX,NY2,NZ2,MAXPTS,IPMAP,IHEX,IR,ILK,SIDE,GMAPX,
     1 GMAPY,GMAPZ,IMPX0,NX2,NY2,NZ2,MAT,NEL,NCODE,ICODE,ZCODE,ISPLX,
     2 ISPLY,ISPLZ,ISPLH,ISPLL)
      IF((NEL.NE.NX2*NY2*NZ2).AND.(IHEX.EQ.0))CALL XABORT('@RESGEO: WR'
     1 // 'ONG GEOMETRY.')
      IF((NEL.NE.NX2*NZ2).AND.(IHEX.NE.0))CALL XABORT('@RESGEO: WRONG ' 
     1 // 'HEXAGONAL GEOMETRY, WRONG NUMBER OF ELEMENTS.')
      DEALLOCATE(MAT,ISPLZ,ISPLY,ISPLX)
      IF(IMPX.GT.2)WRITE(IOUT,*)'CHECKING FUEL-MAP GEOMETRY'
      IF((IGEO.NE.7).AND.(IGEO.NE.9))CALL XABORT('@RESGEO: ONLY '
     1 //'3D-CARTESIAN OR 3D-HEXAGONAL GEOMETRY ALLOWED.')
      IF(IHEX.EQ.0) THEN
        IF((LX.LT.NX).OR.(LY.LT.NY).OR.(LZ.LT.NZ)) THEN
        WRITE(IOUT,*) 'Geometry LX=',LX,', LY=',LY,' and LZ=',LZ, 
     1   ' must be greater or equal to map ',
     2   'NX=',NX,' NY=',NY,' and NZ=',NZ
          CALL XABORT('@RESGEO: WRONG GEOMETRY DEFINITION.')
        ENDIF  
      ELSE
        IF((LX.LT.NX).OR.(LZ.LT.NZ)) THEN
        WRITE(IOUT,*) 'Geometry LX=',LX,' and LZ=',LZ, 
     1   ' must be greater or equal to map ',
     2   'NX=',NX,' and NZ=',NZ
          CALL XABORT('@RESGEO: WRONG GEOMETRY DEFINITION.')
        ENDIF  
      ENDIF 
      IF(NZ.LT.NB)THEN
        WRITE(IOUT,*)'@RESGEO: FOUND NZ=',NZ,' LESS THAN NB=',NB
        CALL XABORT('@RESGEO: WRONG FUEL-MAP GEOMETRY DEFINITION.')
      ENDIF
*----
*  CHECK MESHX OR SIDE
*----
      IF(IGEO.EQ.7) THEN
        CALL XDRSET(GEOXX,LX+1,0.)
        CALL LCMGET(IPMTX,'MESHX',GEOXX)
        DO 10 IMP=1,NX+1
        DO IGM=1,LX+1
          IF(ABS(GMAPX(IMP)-GEOXX(IGM)).LT.EPSI)THEN
            GEOXX(IGM)=GMAPX(IMP)
            GOTO 10
          ENDIF
        ENDDO
        WRITE(IOUT,*)'@RESGEO: MESHX IN L_MAP ',GMAPX(IMP)
        CALL XABORT('@RESGEO: UNABLE TO FIND THIS MESHX IN L_GEOM.')
   10   CONTINUE
        CALL LCMPUT(IPMTX,'MESHX',LX+1,2,GEOXX)
      ELSE IF(IGEO.EQ.9) THEN
        ISPLTL=0
        NY=1
        CALL LCMGET(IPMAP,'SIDE',GMAPSI)
        CALL LCMLEN(IPMAP,'SPLITL',ILONG,ITYLCM)
        IF(ILONG.GT.0) CALL LCMGET(IPMAP,'SPLITL',ISPLTL)
        IF(ISPLTL.EQ.0) ISPLTL=1
        GMAPSI=GMAPSI/REAL(ISPLTL)
        CALL LCMGET(IPMTX,'SIDE',GEOSI)
        IF(ABS(GMAPSI-GEOSI).LT.EPSI)THEN
          GEOSI=GMAPSI
          GOTO 20
        ENDIF
        WRITE(IOUT,*)'@RESGEO: SIDE IN L_MAP ',GMAPSI, GEOSI
        CALL XABORT('@RESGEO: UNABLE TO FIND THIS SIDE IN L_GEOM.')
   20   CONTINUE
        CALL LCMPUT(IPMTX,'SIDE',1,2,GEOSI)
      ENDIF
*----
*  CHECK MESHY (ONLY IF 3D-CARTESIAN GEOMETRY)
*----
      IF(IGEO.EQ.7) THEN
        CALL XDRSET(GEOYY,LY+1,0.)
        CALL LCMGET(IPMTX,'MESHY',GEOYY)
        DO 30 IMP=1,NY+1
        DO IGM=1,LY+1
          IF(ABS(GMAPY(IMP)-GEOYY(IGM)).LT.EPSI)THEN
            GEOYY(IGM)=GMAPY(IMP)
            GOTO 30
          ENDIF
        ENDDO
        WRITE(IOUT,*)'@RESGEO: MESHY IN FUEL MAP ',GMAPY(IMP)
        CALL XABORT('@RESGEO: UNABLE TO FIND THIS MESHY IN L_GEOM.')
   30   CONTINUE
        CALL LCMPUT(IPMTX,'MESHY',LY+1,2,GEOYY)
      ENDIF
*----
*  CHECK MESHZ
*----
      CALL XDRSET(GEOZZ,LZ+1,0.)
      CALL LCMGET(IPMTX,'MESHZ',GEOZZ)
      DO 50 IMP=1,NZ+1
      DO IGM=1,LZ+1
        IF(ABS(GMAPZ(IMP)-GEOZZ(IGM)).LT.EPSI)THEN
          GEOZZ(IGM)=GMAPZ(IMP)
          GOTO 50
        ENDIF
      ENDDO
      WRITE(IOUT,*)'@RESGEO: MESHZ IN FUEL MAP ',GMAPZ(IMP)
      CALL XABORT('@RESGEO: UNABLE TO FIND THIS MESHZ IN L_GEOM.')
   50 CONTINUE
      CALL LCMPUT(IPMTX,'MESHZ',LZ+1,2,GEOZZ)
      DEALLOCATE(GMAPZ,GMAPY,GMAPX)
*----
*  CHECK FUEL MIXTURES
*----
      CALL RESPFM(IPMAP,IPMTX,NX,NY,NZ,LX,LY,LZ,NFUEL,IMPX,IGEO,NCH,NB,
     1 NTOT)
      RETURN
      END
