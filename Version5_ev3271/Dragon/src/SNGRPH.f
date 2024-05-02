*DECK SNGRPH
      SUBROUTINE SNGRPH (IMPX,NMZ,NHEX,LZ,NBC,MAXCELLCNT,MAXWAVECNT,
     1 IZGLOB,CONNEC,TASKLIST,TASKSINWAVE,NWAVES,MAXLOZ,CON_WHOAMI,
     2 CON_TOWHOM)
*
*-----------------------------------------------------------------------
*
*Purpose:
* TBBuild graphs for parallel sweep on hexagonal geometries.
*
*Copyright:
* Copyright (C) 2021 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. A. Calloo
*
*Parameters: input
* MAXCELLCNT is APPCELLCNT for 3D !!!
*
*
*Parameters: output
* 
*
*Comments:
* The lozenge under consideration is given by the position within the
*         the matrix. See user manual and/or data manual and/or thesis
*                                 _____
*                                /   / \
*                               / B /   \
*                         ,----(----  A  )----.
*                        /      \ C \   /      \
*                       /        \___\_/        \
*                       \   4    /     \   2    /
*                        \      /       \      /
*                         )----(   1     )----(
*                        /      \       /      \
*                       /        \_____/        \
*                       \   5    /     \   7    /
*                        \      /       \      /
*                         `----(   6     )----'
*                               \       /
*                                \_____/
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IMPX,NMZ,NHEX,LZ,NBC,MAXCELLCNT,MAXWAVECNT,IZGLOB(NHEX,6)
      INTEGER, DIMENSION(3,(NHEX*3)*2,6) :: CONNEC
      INTEGER, DIMENSION(2,MAXCELLCNT,MAXWAVECNT,6) :: TASKLIST
      INTEGER, DIMENSION(MAXWAVECNT,6) :: TASKSINWAVE
*----
*  LOCAL VARIABLES
*----
      INTEGER IND,REM,IFACE(3),LOZ(3),SIDE(3),PRIHEX,STRTHEX,HEXCW,
     1 HEXAW,CELLCNT,WAVECNT,CNTR1,CNTR2,NWAVES,MAXLOZ
      INTEGER CSIDx,CLOZx,CHEXx,CSIDy,CLOZy,CHEXy
      INTEGER TOTTASK
      LOGICAL, DIMENSION(3,3,NHEX,NMZ) :: DOMSWP

      INTEGER :: CON_WHOAMI(2,(NHEX*3)*2,6), CON_TOWHOM(2,(NHEX*3)*2,6)
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: HEXSWP2
      INTEGER, ALLOCATABLE, DIMENSION(:)     :: HEXINCOL
      INTEGER BFPOS,CBFPOS1,CBFPOS2,SIZEBFLUX


*
*                               4
*                           ________
*                          /        \
*                     3   /          \   2
*                        /            \
*                       (              )
*                        \            /
*                     5   \          /   1
*                          \________/
*                                         
*                              6 
      
      IF((NMZ.GT.LZ)) CALL XABORT('SNGRPH: MORE MACROCELLS' 
     1   //' ALONG Z AXIS THAN ACTUAL CELLS. (1)')
      NCELZ=CEILING(REAL(LZ)/NMZ)
      IF((NCELZ*(NMZ-1)).EQ.LZ) CALL XABORT('SNGRPH: THIS WILL ONLY'
     1   //' CORRESPOND TO (NUMBER MACROCELLS ASKED MINUS ONE). EITHER '
     2   //' INCREASE SPLIT OR DECREASE MACROCELL NUMBER BY 1.')
      
      NCOL=2*NBC -1
      ALLOCATE(HEXINCOL(NCOL))
      DO ICOL=1,NCOL
         HEXINCOL(ICOL)=NBC-1+ICOL
         IF(ICOL.EQ.CEILING(REAL(NCOL)/2)) HEXINCOL(ICOL)=NCOL
         IF(ICOL.GT.CEILING(REAL(NCOL)/2)) HEXINCOL(ICOL)=
     1      NBC-1+(NCOL+1-ICOL)
      ENDDO

      ALLOCATE(HEXSWP2(NCOL,NCOL,6))
      HEXSWP2(:NCOL,:NCOL,:6)=0

      DO IND=1,6
         ICOUNT=0
         DO ICOL=1,NCOL
            DO IH=1,HEXINCOL(ICOL)
               ICOUNT=ICOUNT+1
               HEXSWP2(IH,ICOL,IND)=IZGLOB(ICOUNT,IND)
            ENDDO
         ENDDO
      ENDDO
      CON_WHOAMI(:2,:((NHEX*3)*2),:6)=0
      CON_TOWHOM(:2,:((NHEX*3)*2),:6)=0
*----
* MAIN LOOP OVER DIRECTIONS
*----
      DO IND=1,6

         IFACE(:3) = 0
         DOMSWP(:3,:3,:NHEX,:NMZ)=.FALSE.
*----
* Identify incoming sides for fluxes
*----
         IF(IND.EQ.1) THEN
            IFACE(1) = 3
            IFACE(2) = 5
            IFACE(3) = 6
         ELSE IF(IND.EQ.2) THEN
            IFACE(1) = 5
            IFACE(2) = 6
            IFACE(3) = 1
         ELSE IF(IND.EQ.3) THEN
            IFACE(1) = 6
            IFACE(2) = 1
            IFACE(3) = 2
         ELSE IF(IND.EQ.4) THEN
            IFACE(1) = 1
            IFACE(2) = 2
            IFACE(3) = 4
         ELSE IF(IND.EQ.5) THEN
            IFACE(1) = 2
            IFACE(2) = 4
            IFACE(3) = 3
         ELSE IF(IND.EQ.6) THEN
            IFACE(1) = 4
            IFACE(2) = 3
            IFACE(3) = 5
         ENDIF

* PRIHEX  - hexagon where sweep starts, with three sides of incoming fluxes
* STRTHEX - hexagon with lowest numbering in outermost shell

         PRIHEX  = IZGLOB(1,IND)
         STRTHEX = NHEX -(NBC -1)*6 +1
*        * PRIMARY INCOMING HEXAGON SIDES
         DO IFC=1,3
            REM      = MOD(IFACE(IFC),2)
            SIDE(IFC) = 1
            LOZ(IFC)  = IFACE(IFC)/2
            IF(REM.EQ.0) SIDE(IFC) = 2
            IF(REM.EQ.1) LOZ(IFC)  = LOZ(IFC) +1
            DOMSWP(SIDE(IFC),LOZ(IFC),PRIHEX,:) = .TRUE.
         ENDDO 

*        * SECONDARY INCOMING HEXAGONS SIDES
         HEXCW  = PRIHEX
         HEXAW  = PRIHEX
         DO IH=1,NBC-1
            HEXCW  = HEXCW -1
            IF(HEXCW.EQ.STRTHEX -1) HEXCW = NHEX
            DO IFC=1,2
               DOMSWP(SIDE(IFC),LOZ(IFC),HEXCW,:) = .TRUE.
            ENDDO
*           *
            HEXAW  = HEXAW +1
            IF(HEXAW.EQ.NHEX +1) HEXAW = STRTHEX
            DO IFC=2,3
               DOMSWP(SIDE(IFC),LOZ(IFC),HEXAW,:) = .TRUE.
            ENDDO
         ENDDO

*        * TERTIARY INCOMING HEXAGONS SIDES
         DO IH=1,NBC-1
            HEXCW  = HEXCW -1
            IF(HEXCW.EQ.STRTHEX-1) HEXCW = NHEX
            DOMSWP(SIDE(1),LOZ(1),HEXCW,:) = .TRUE.
*           *
            HEXAW  = HEXAW +1
            IF(HEXAW.EQ.NHEX+1) HEXAW = STRTHEX
            DOMSWP(SIDE(3),LOZ(3),HEXAW,:) = .TRUE.
         ENDDO
*        * TOP (OR BOTTOM) INCOMING HEXAGONS SIDES
         DOMSWP(3,:,:,1) = .TRUE.
*----
* Identify which (and how many) lozenges can be processed in each 'wave'
*----
         TOTTASK = 0
         WAVECNT = 0
         DO WHILE (TOTTASK.LT.(NHEX*NMZ))
            CELLCNT = 0
            WAVECNT = WAVECNT +1
            IF(WAVECNT.GT.MAXWAVECNT) CALL XABORT ('SNGRPH: MORE'
     1         //'WAVES OF MACROCELLS THAN EXPECTED. THIS SHOULD '
     2         //'NOT HAPPEN. IS PROBABLY INDICATIVE OF A BUG.')

            ! Iterate over domain to find next cells in current
            ! wavefront.
            DO IZ=1,NMZ
            DO IH=1,NHEX

               NUM_INC_SIDES = 0
               DO ILOZ=1,3
               DO ISIDE=1,3
                  IF((DOMSWP(ISIDE,ILOZ,IH,IZ).EQV..TRUE.))THEN
                     NUM_INC_SIDES=NUM_INC_SIDES+1
                  ENDIF
               ENDDO
               ENDDO

               IF(NUM_INC_SIDES.EQ.6)THEN
                  TOTTASK = TOTTASK +1
                  CELLCNT = CELLCNT +1
                  TASKLIST(1,CELLCNT,WAVECNT,IND)=IH
                  TASKLIST(2,CELLCNT,WAVECNT,IND)=IZ
                  DO ILOZ=1,3
                  DO ISIDE=1,3
                     DOMSWP(ISIDE,ILOZ,IH,IZ)=.FALSE.
                     DOMSWP(ISIDE,ILOZ,IH,IZ)=.FALSE.
                     DOMSWP(ISIDE,ILOZ,IH,IZ)=.FALSE.
                  ENDDO
                  ENDDO
               ENDIF
               
               IF(CELLCNT.GT.MAXCELLCNT) CALL XABORT ('SNGRPH: MORE '
     1            //'MACROCELLS THAN EXPECTED IN CURRENT WAVEFRONT. '
     2            //'THIS SHOULD NOT HAPPEN. IS PROBABLY INDICATIVE '
     3            //'OF A BUG.')

            ENDDO
            ENDDO 

            TASKSINWAVE(WAVECNT,IND)=CELLCNT

            DO ITASK=1,CELLCNT
            DO ILOZ=1,3
               IH   = TASKLIST(1,ITASK,WAVECNT,IND)
               IZ   = TASKLIST(2,ITASK,WAVECNT,IND)

               CHEXx = CONNEC(1,(IH-1)*3*2 +(ILOZ-1)*2 +1,IND)
               CLOZx = CONNEC(2,(IH-1)*3*2 +(ILOZ-1)*2 +1,IND)
               CSIDx = CONNEC(3,(IH-1)*3*2 +(ILOZ-1)*2 +1,IND)

               CHEXy = CONNEC(1,(IH-1)*3*2 +(ILOZ-1)*2 +2,IND)
               CLOZy = CONNEC(2,(IH-1)*3*2 +(ILOZ-1)*2 +2,IND)
               CSIDy = CONNEC(3,(IH-1)*3*2 +(ILOZ-1)*2 +2,IND)

               IF(CHEXx.LE.NHEX) DOMSWP(CSIDx,CLOZx,CHEXx,IZ)=.TRUE.
               IF(CHEXy.LE.NHEX) DOMSWP(CSIDy,CLOZy,CHEXy,IZ)=.TRUE.
               IF(IZ.LT.NMZ) DOMSWP(3,ILOZ,IH,IZ+1)=.TRUE.
            ENDDO
            ENDDO

         ENDDO
  

******************************************************************
******************************************************************
******************************************************************
! This was designed for 2D KBA-style sweep. However, because the
! connectivity in each plane of hexagons obviously remains the 
! same in 3D, this works in 3D too. It just grinds a bit more 
! uselessly on the loops but that's negligible compared to overall
! computational time. 

         SIZEBFLUX = (MAXCELLCNT*2) -1

         DO IWAVE=1,MAXWAVECNT
            DO ITASK=1,TASKSINWAVE(IWAVE,IND)
               BFPOS=0
               DO ILOZ=1,3

                  ISIDEx=0
                  ISIDEy=0
                  IF((IND.EQ.1).OR.(IND.EQ.4))THEN
                     IF(ILOZ.EQ.1)THEN
                        ISIDEx=1
                        ISIDEy=2
                     ELSEIF(ILOZ.EQ.2)THEN
                        ISIDEx=1
                        ISIDEy=3
                     ELSEIF(ILOZ.EQ.3)THEN
                        ISIDEx=2
                        ISIDEy=3
                     ENDIF
                  ELSEIF((IND.EQ.2).OR.(IND.EQ.5))THEN
                     IF(ILOZ.EQ.1)THEN
                        ISIDEx=3
                        ISIDEy=1
                     ELSEIF(ILOZ.EQ.2)THEN
                        ISIDEx=3
                        ISIDEy=2
                     ELSEIF(ILOZ.EQ.3)THEN
                        ISIDEx=1
                        ISIDEy=2
                     ENDIF
                  ELSEIF((IND.EQ.3).OR.(IND.EQ.6))THEN
                     IF(ILOZ.EQ.1)THEN
                        ISIDEx=2
                        ISIDEy=3
                     ELSEIF(ILOZ.EQ.2)THEN
                        ISIDEx=2
                        ISIDEy=1
                     ELSEIF(ILOZ.EQ.3)THEN
                        ISIDEx=3
                        ISIDEy=1
                     ENDIF
                  ENDIF

                  IH   = TASKLIST(1,ITASK,IWAVE,IND)
                  IZ   = TASKLIST(2,ITASK,IWAVE,IND)

                  BFPOS=0
                  DO ICOL=1,NCOL
                     DO IHIC=1,HEXINCOL(ICOL)
                        
                        IND2=IND+1
                        IF(IND2.GT.6) IND2=1

                        IF(IH.EQ.HEXSWP2(IHIC,ICOL,IND2))THEN
                           BFPOS=NCOL-ICOL+1
                           EXIT
                        ENDIF
                     ENDDO
                     IF(BFPOS.NE.0) EXIT
                  ENDDO

                  CON_WHOAMI(1,(IH-1)*3*2 +(ILOZ-1)*2 +1,IND)=BFPOS
                  CON_WHOAMI(2,(IH-1)*3*2 +(ILOZ-1)*2 +1,IND)=ISIDEx
                  CON_WHOAMI(1,(IH-1)*3*2 +(ILOZ-1)*2 +2,IND)=BFPOS
                  CON_WHOAMI(2,(IH-1)*3*2 +(ILOZ-1)*2 +2,IND)=ISIDEy


                  CHEXx = CONNEC(1,(IH-1)*3*2 +(ILOZ-1)*2 +1,IND)
                  CHEXy = CONNEC(1,(IH-1)*3*2 +(ILOZ-1)*2 +2,IND)

                  CBFPOS1=0
                  CBFPOS2=0
                  DO ICOL=1,NCOL
                     DO IHIC=1,HEXINCOL(ICOL)
                        
                        IND2=IND+1
                        IF(IND2.GT.6) IND2=1

                        IF(CHEXx.EQ.HEXSWP2(IHIC,ICOL,IND2))THEN
                           CBFPOS1=NCOL-ICOL+1
                        ENDIF
                        IF(CHEXy.EQ.HEXSWP2(IHIC,ICOL,IND2))THEN
                           CBFPOS2=NCOL-ICOL+1
                           EXIT
                        ENDIF
                     ENDDO
                     IF((CBFPOS1.NE.0).AND.(CBFPOS2.NE.0)) EXIT
                  ENDDO

                  CON_TOWHOM(1,(IH-1)*3*2 +(ILOZ-1)*2 +1,IND)=CBFPOS1
                  CON_TOWHOM(2,(IH-1)*3*2 +(ILOZ-1)*2 +1,IND)=ISIDEx
                  CON_TOWHOM(1,(IH-1)*3*2 +(ILOZ-1)*2 +2,IND)=CBFPOS2
                  CON_TOWHOM(2,(IH-1)*3*2 +(ILOZ-1)*2 +2,IND)=ISIDEy

               ENDDO
            ENDDO
         ENDDO

******************************************************************
******************************************************************
******************************************************************

      ENDDO

*----
* Identify actual number of 'waves' and largest number of tasks/lozenges
* to be dealt with in a wave
*----
      MAXLOZ =0
      NWAVES =0
      DO IND=1,6
         CNTR1=0
         CNTR2=MAXVAL(TASKSINWAVE(:,IND))
         DO I=1,MAXWAVECNT
            IF(TASKSINWAVE(I,IND).EQ.0) EXIT
            CNTR1=CNTR1 +1
         ENDDO
         IF(IND.EQ.1)THEN
            NWAVES = CNTR1
            MAXLOZ = CNTR2
         ELSE
*           * INCONSISTENT VALUE ERROR CHECKS
            IF(NWAVES.NE.CNTR1) THEN
               WRITE(*,*) NWAVES, CNTR1
               CALL XABORT('SNGRPH: NUMBER OF WAVES '
     1        //'DIFFERENT PER DIRECTION. UNEXPECTED. THIS SHOULD BE '
     2        //'INVESTIGATED.')
            ENDIF
            IF(MAXLOZ.NE.CNTR2) CALL XABORT('SNGRPH: MAXIMUM NUMBER '
     1        //'OF TASKS DIFFERENT PER DIRECTION. UNEXPECTED. THIS '
     2        //'SHOULD BE INVESTIGATED.')
         ENDIF
*        * ZERO VALUE ERROR CHECKS
         IF(NWAVES.EQ.0) CALL XABORT('SNGRPH: GRAPH NOT WORKING '
     1    //'PROPERLY. THIS SHOULD BE INVESTIGATED. (1)')
         IF(MAXLOZ.EQ.0) CALL XABORT('SNGRPH: GRAPH NOT WORKING '
     1    //'PROPERLY. THIS SHOULD BE INVESTIGATED. (2)')
      ENDDO

      DO IND=1,6
         DO I=1,MAXWAVECNT
            WRITE(*,*) TASKLIST(:,:,I,IND)
         ENDDO
         WRITE(*,*) ''
      ENDDO

      DO I=1,MAXWAVECNT
         WRITE(*,*) TASKSINWAVE(I,1)
      ENDDO
      WRITE(*,*) ''
      WRITE(*,*) ''
      WRITE(*,*) ''


*----
*  PRINT A FEW GEOMETRY CHARACTERISTICS
*----
      IF(IMPX.GT.3) THEN
         write(6,*) ' '
         write(6,*) 'NBC   =',NBC
      ENDIF
*
      RETURN
      END
