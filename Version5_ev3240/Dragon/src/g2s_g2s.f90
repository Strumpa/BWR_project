!
!-----------------------------------------------------------------------
!
!Purpose:
! Generate a surfacic 2D geometry following the TDT specification.
!
!Copyright:
! Copyright (C) 2001 Ecole Polytechnique de Montreal.
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
!Author(s):
! G. Civario (CS-SI)
!
!Parameters: input/output
!  NENTRY : NUMBER OF LINKED LISTS AND FILES USED BY THE MODULE.
!  HENTRY : CHARACTER*12 NAME OF EACH LINKED LIST OR FILE.
!  IENTRY : =1 LINKED LIST; =2 XSM FILE; =3 SEQUENTIAL BINARY FILE;
!           =4 SEQUENTIAL ASCII FILE; =5 DIRECT ACCESS FILE.
!  JENTRY : =0 THE LINKED LIST OR FILE IS CREATED;
!           =1 THE LINKED LIST OR FILE IS OPEN FOR MODIFICATIONS;
!           =2 THE LINKED LIST OR FILE IS OPEN IN READ-ONLY MODE.
!  KENTRY : FILE UNIT NUMBER OR LINKED LIST ADDRESS.
!
!-----------------------------------------------------------------------
!
subroutine G2S(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
  use GANLIB
  use SALGET_FUNS_MOD
  use celluleBase
  use cellulePlaced
  use boundCond
  use ptNodes
  use pretraitement
  use derivedPSPLOT
  use track
  use segArc
  use generTabSegArc
  use generSAL

  implicit none

  integer NENTRY
  integer IENTRY,JENTRY
  type(c_ptr) KENTRY
  character*12 HENTRY
  dimension IENTRY(*),JENTRY(*),KENTRY(*),HENTRY(*)


  integer,parameter :: dimTabCelluleBase = 20000
  integer,parameter :: dimTabSegArc = 100000

  type(c_ptr)  :: ipGeo,ipGeo_1
  integer      :: sizeB,sizeP,sizeSA,nbNode,nbCLP,ilong,nbFlux,ipSal,ipPs,ipAl, &
                  ipZa,indic,nitma,impx
  character(len=12) :: name_geom,nammy,text12
  logical      :: empty,lcm,drawNod,drawMix
  real,dimension(2) :: zoomx,zoomy
  integer,allocatable,dimension(:) :: gig,merg
  integer,dimension(10) :: datain
  real :: flott
  double precision :: dflott
  integer      :: lgMaxGig=0
  !
  ipAl=-1 ! no Alamos file read
  ipZa=-1 ! no PropertyMap file read
  ipGeo_1=c_null_ptr ! no geometry read
  if ((NENTRY == 2).and.(IENTRY(2) == 4)) then
     !generating ps file from Salomon file
     ipPs   = FILUNIT(KENTRY(1)) ! output psfile
     g_psp_isEpsFile = (index(HENTRY(1),'.eps')/=0) !is it an eps file ?
     ipSal  = FILUNIT(KENTRY(2)) ! input Salomon/Alamos file
     ! check argument types and permissions
     if ((IENTRY(1) /= 4) .or. (JENTRY(1) /= 0)) &
          call XABORT('G2S: a new file was expected for the postscript file')
     if ((IENTRY(2) /= 4) .or. (JENTRY(2) /= 2)) &
          call XABORT('G2S: expecting Salomon file in read-only mode')
  else if ((NENTRY == 3).and.(IENTRY(3) == 4)) then
     !generating Salomon and ps files from Alamos file
     ipSal  = FILUNIT(KENTRY(1)) ! output Salomon file
     ipPs   = FILUNIT(KENTRY(2)) ! output psfile
     ipAl   = FILUNIT(KENTRY(3)) ! input Alamos file
     g_psp_isEpsFile = (index(HENTRY(2),'.eps')/=0) !is it an eps file ?
     ! check argument types and permissions
     if ((IENTRY(1) /= 4) .or. (JENTRY(1) /= 0)) &
          call XABORT('G2S: a new file was expected for the Salomon file')
     if ((IENTRY(2) /= 4) .or. (JENTRY(2) /= 0)) &
          call XABORT('G2S: a new file was expected for the postscript file')
     if ((IENTRY(3) /= 4) .or. (JENTRY(3) /= 2)) &
          call XABORT('G2S: expecting Alamos file in read-only mode')
  else if ((NENTRY == 4).and.(IENTRY(3) == 4).and.(IENTRY(4) == 4)) then
     !generating Salomon and ps files from Alamos and PropertyMap file
     ipSal  = FILUNIT(KENTRY(1)) ! output Salomon file
     ipPs   = FILUNIT(KENTRY(2)) ! output psfile
     ipAl   = FILUNIT(KENTRY(3)) ! input Alamos file
     ipZa   = FILUNIT(KENTRY(4)) ! input PropertyMap file
     g_psp_isEpsFile = (index(HENTRY(2),'.eps')/=0) !is it an eps file ?
     ! check argument types and permissions
     if ((IENTRY(1) /= 4) .or. (JENTRY(1) /= 0)) &
          call XABORT('G2S: a new file was expected for the Salomon file')
     if ((IENTRY(2) /= 4) .or. (JENTRY(2) /= 0)) &
          call XABORT('G2S: a new file was expected for the postscript file')
     if ((IENTRY(3) /= 4) .or. (JENTRY(3) /= 2)) &
          call XABORT('G2S: expecting Alamos file in read-only mode')
     if ((IENTRY(4) /= 4) .or. (JENTRY(4) /= 2)) &
          call XABORT('G2S: expecting PropertyMap file in read-only mode')
  else if ((NENTRY == 2).and.(IENTRY(2) <= 2)) then
     !generating Salomon file from LCM geometry
     ipPs   = -1        ! no postscript file
     ipSal  = FILUNIT(KENTRY(1)) ! output Salomon file
     ipGeo_1= KENTRY(2) ! geometry read
     ! check argument types and permissions
     if ((IENTRY(1) /= 4) .or. (JENTRY(1) /= 0)) &
          call XABORT('G2S: a new ASCII file was expected for writing geometry')
     if ((IENTRY(2) > 2) .or. (JENTRY(2) /= 2)) &
          call XABORT('G2S: expecting LCM geometry in read-only mode(1)')
  else if ((NENTRY == 3).and.(IENTRY(3) <= 2)) then
     !generating Salomon file and ps file from LCM geometry
     ipSal  = FILUNIT(KENTRY(1)) ! output Salomon file
     ipPs   = FILUNIT(KENTRY(2)) ! output psfile
     g_psp_isEpsFile = (index(HENTRY(2),'.eps')/=0) !is it an eps file ?
     ipGeo_1= KENTRY(3) ! geometry read
     ! check argument types and permissions
     if ((IENTRY(1) /= 4) .or. (JENTRY(1) /= 0)) &
          call XABORT('G2S: a new file was expected for writing geometry')
     if ((IENTRY(2) /= 4) .or. (JENTRY(2) /= 0)) &
          call XABORT('G2S: a new file was expected for the postscript file')
     if ((IENTRY(3) > 2) .or. (JENTRY(3) /= 2)) &
          call XABORT('G2S: expecting LCM geometry in read-only mode(2)')
  else
     call XABORT('G2S: you must provide 2, 3 or 4 arguments')
  end if

  impx=1
  drawNod = .false.
  drawMix = .false.
  zoomx = (/ 0.0, 1.0 /)
  zoomy = (/ 0.0, 1.0 /)
  typgeo=0
  10 call REDGET(indic,nitma,flott,text12,dflott)
  if (indic == 10) go to 20
  if (indic /= 3) call XABORT('G2S: character data expected.')
  if (text12 == 'EDIT') then
    ! read the print index.
    call REDGET(indic,impx,flott,text12,dflott)
    if (indic /= 1) call XABORT('G2S: integer data expected(1).')
  else if (text12 == 'DRAWNOD') then
     drawNod=.true.
     drawmix=.true.
  else if (text12 == 'DRAWMIX') then
     drawNod=.true.
     drawmix=.false.
  else if (text12 == 'ZOOMX') then
    call REDGET(indic,nitma,zoomx(1),text12,dflott)
    if (indic /= 2) call XABORT('G2S: real data expected(1).')
    call REDGET(indic,nitma,zoomx(2),text12,dflott)
    if (indic /= 2) call XABORT('G2S: real data expected(2).')
    if ((zoomx(1).lt.0.0).or.(zoomx(2).le.zoomx(1)).or.(zoomx(2).gt.1.0)) then
      call XABORT('G2S: invalid zoom factors in x.')
    endif
  else if (text12 == 'ZOOMY') then
    call REDGET(indic,nitma,zoomy(1),text12,dflott)
    if (indic /= 2) call XABORT('G2S: real data expected(3).')
    call REDGET(indic,nitma,zoomy(2),text12,dflott)
    if (indic /= 2) call XABORT('G2S: real data expected(4).')
    if ((zoomy(1).lt.0.0).or.(zoomy(2).le.zoomy(1)).or.(zoomy(2).gt.1.0)) then
      call XABORT('G2S: invalid zoom factors in y.')
    endif
  else if (text12 == 'ALAMOS') then
    if (NENTRY == 2) call XABORT('G2S: three entries required.')
    if (ipAl == -1) call XABORT('G2S: no RHS Salomon file.')
    call REDGET(indic,typgeo,flott,text12,dflott)
    if (indic /= 1) call XABORT('G2S: integer data expected(2).')
  else if (text12 == ';') then
     go to 20
  else
     call XABORT('G2S: '//text12//' is an invalid keyword.')
  end if
  go to 10

  !conversion of Alamos file into a Salomon file
  20 if (ipAl /= -1) call g2s_convert(impx,ipAl,ipZa,ipSal)

  sizeB = 0   !cellules de base
  sizeP = 0   !cellules placees
  sizeSA = 0  !elements geometriques
  if (c_associated(ipGeo_1)) then
     ! lecture of the geometry name
     call LCMINF(ipGeo_1,name_geom,nammy,empty,ilong,lcm)
  
     ! copy the input geometric object
     call lcmop(ipGeo,'geom_copy',0,1,0)
     call lcmequ(ipGeo_1,ipGeo)

     !initialisation des differents tableaux
     call initializeData(dimTabCelluleBase,dimTabSegArc)

     !unfold the geometry
     call g2s_unfold(ipGeo,impx)

     !pretraitement des donnees lues (remplace la partie python)
     !+completion des cellules de base et remplissage du tableau
     !des cellules placees
     call prepareData(ipGeo,sizeB,sizeP,lgMaxGig)
     if (impx > 0) write(*,*) 'fin   : prepareData lgMaxGig=',lgMaxGig

     !en sortie, toutes les cellules de base ont tous leurs
     !champs remplis, et le tableau des cellules placees est pret

     !eclatement des cellules
     call splitCells(sizeP,sizeSA)
 
     !creation de nouveaux segements aux interfaces des cellules
     !et elimination des doublons
     call addSegsAndClean(sizeSA)

     !prise en compte des conditions aux limites
     call appliBoundariConditions(ipGeo,sizeSA,nbCLP)

     !calcul des nodes delimites par les elements
     allocate(merg(dimTabCelluleBase),stat=alloc_ok)
     if (alloc_ok /= 0) call XABORT("G2S: g2s_g2s(1) => allocation pb(1)")
     call createNodes(sizeSA,dimTabCelluleBase,nbNode,merg)
     if (sizeSA > dimTabSegArc) call XABORT('g2s_g2s: sizeSA overflow')
     
     !calcul des arrays gig et merg
     allocate(gig(nbNode*lgMaxGig),stat=alloc_ok)
     if (alloc_ok /= 0) call XABORT("G2S: g2s_g2s(1) => allocation pb(2)")
     call generateTrack(sizeP,sizeSA,nbNode,lgMaxGig,gig,merg)
     nbFlux=maxval(merg(:nbNode))
  else
     if (JENTRY(NENTRY) == 0) call XABORT('G2S: a RHS ascii file is expected')
     !initialisation de TabSegArc
     call SALGET(datain,4,ipSal,0,'dimensions for geometry')
     nbNode=datain(3)
     sizeSA=datain(4)
     rewind(ipSal)
     allocate(TabSegArc(sizeSA),medium(nbNode),stat=alloc_ok)
     if (alloc_ok /= 0) call XABORT("G2S: generateTabSegArc => allocation pb")
     call initializebCData()  
     allocate(merg(nbNode),stat=alloc_ok)
     if (alloc_ok /= 0) call XABORT("G2S: g2s_g2s(2) => allocation pb")
     call generateTabSegArc(ipSal,sizeSA,nbNode,nbCLP,nbFlux,merg,name_geom,impx)
  endif

  !impression des segArc charges
  if (ipPs /= -1) call drawSegArc(ipPs,sizeSA,drawMix,drawNod,zoomx,zoomy)

  if (c_associated(ipGeo_1)) then
     !creation du fichier de commande SAL
     call generateSALFile(ipSal,sizeSA,nbNode,nbCLP,nbFlux,merg,name_geom)
     deallocate(gig)
     call LCMCL(ipGeo,2)
  endif
  deallocate(merg)

  write(6,*) "  At end of G2S:"
  write(6,*) "    ",sizeSA,"segs or arcs"
  write(6,*) "    ",nbNode,"nodes"
  write(6,*) "    ",nbCLP,"boundary conditions other than default"

  !liberation de la memoire allouee
  call destroyData(sizeB,sizeP)
end subroutine G2S

subroutine initializeData(dimTabCelluleBase,dimTabSegArc)
  use celluleBase
  use cellulePlaced
  use boundCond
  use segArc
  implicit none
  integer,intent(in) :: dimTabCelluleBase,dimTabSegArc

  call initializeTabCelluleBase(dimTabCelluleBase)
  call initializeTabCellulePlaced()
  allocate(tabSegArc(dimTabSegArc),medium(dimTabCelluleBase))
  call initializebCData()  
end subroutine initializeData

subroutine destroyData(szB,szP)
  use celluleBase
  use cellulePlaced
  use boundCond
  use segArc
  implicit none
  integer,intent(in) :: szB,szP

!  if (szB /= 0) call destroyTabCelluleBase(szB)
  deallocate(medium,tabSegArc)
  if (szB /= 0) deallocate(TabCelluleBase)
  if (szP /= 0) call destroyTabCellulePlaced(szP)
    call destroybCData()
end subroutine destroyData
