!
!-----------------------------------------------------------------------
!
!Purpose:
! Generate a dataset for use in a Monte Carlo code.
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
subroutine G2MC(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
  use SALGET_FUNS_MOD
  use celluleBase
  use cellulePlaced
  use boundCond
  use ptNodes
  use pretraitement
  use derivedPSPLOT
  use monteCarlo
  use track
  use segArc
  use GANLIB
  use generTabSegArc

  implicit none

  integer NENTRY
  integer IENTRY,JENTRY
  type(c_ptr) KENTRY
  character*12 HENTRY
  dimension IENTRY(*),JENTRY(*),KENTRY(*),HENTRY(*)

  integer,parameter :: dimTabCelluleBase = 20000
  integer,parameter :: dimTabSegArc = 100000

  type(c_ptr) :: ipGeo,ipGeo_1
  integer     :: ipMC,ipSal,ipPs,sizeB,sizeP,sizeSA,nbNode,nbCLP,nbFlux,indic, &
                 & nitma,impx
  real :: flott
  double precision :: dflott
  integer      :: lgMaxGig=0
  integer,dimension(10) :: datain
  integer,allocatable,dimension(:) :: merg
  character(len=12) :: name_geom,text12
  logical      :: drawNod,drawMix
  real,dimension(2) :: zoomx,zoomy

  ipGeo_1=c_null_ptr ! no geometry read
  ipSal=-1 ! no Salomon file read
  ipMC   = FILUNIT(KENTRY(1)) ! Monte-Carlo file generated
  if ((NENTRY == 2).and.(IENTRY(2) == 4)) then
     !generating Monte-Carlo file from Salomon file
     ipSal  = FILUNIT(KENTRY(2)) ! input Salomon file (surfacic elements)
     ipPs   = -1         ! no postscript file
     ! check that second argumnet is file to write
     ! then the tracking object
     if ((IENTRY(1) /= 4) .or. (JENTRY(1) /= 0)) &
        call XABORT('G2MC: a new ascii file expected at LHS for containing MC info')
     if ((IENTRY(2) /= 4) .or. (JENTRY(2) /= 2)) &
        call XABORT('G2MC: read-only ascii file expected at RHS with surfacic elements')
  else if ((NENTRY == 3).and.(IENTRY(3) == 4)) then
     !generating Monte-Carlo and ps files from Salomon file
     ipPs   = FILUNIT(KENTRY(2)) ! output psfile
     ipSal  = FILUNIT(KENTRY(3)) ! input Salomon file (surfacic elements)
     g_psp_isEpsFile = (index(HENTRY(2),'.eps')/=0) !is it an eps file ?
     ! check argument types and permissions
     if ((IENTRY(1) /= 4) .or. (JENTRY(1) /= 0)) &
          call XABORT('G2MC: a new file was expected for the Monte-Carlo file')
     if ((IENTRY(2) /= 4) .or. (JENTRY(2) /= 0)) &
          call XABORT('G2MC: a new file was expected for the postscript file')
     if ((IENTRY(3) /= 4) .or. (JENTRY(3) /= 2)) &
          call XABORT('G2MC: expecting Salomon file in read-only mode')
  else if ((NENTRY == 2).and.(IENTRY(2) <= 2)) then
     !generating Monte-Carlo file from LCM geometry
     ipPs   = -1        ! no postscript file
     ipGeo_1= KENTRY(2) ! input geometry
     ! check argument types and permissions
     if ((IENTRY(1) /= 4) .or. (JENTRY(1) /= 0)) &
        call XABORT('G2MC: a new ascii file expected at LHS for containing MC info')
  else if ((NENTRY == 3).and.(IENTRY(3) <= 2)) then
     !generating Monte-Carlo file and ps file from LCM geometry
     ipPs   = FILUNIT(KENTRY(2)) ! output psfile
     g_psp_isEpsFile = (index(HENTRY(2),'.eps')/=0) !is it an eps file ?
     ipGeo_1= KENTRY(3) ! input geometry
     ! check argument types and permissions
     if ((IENTRY(1) /= 4) .or. (JENTRY(1) /= 0)) &
        call XABORT('G2MC: a new ascii file expected at LHS for containing MC info')
     if ((IENTRY(2) /= 4) .or. (JENTRY(2) /= 0)) &
        call XABORT('G2MC: a new file was expected for the postscript file')
  else
     call XABORT('G2MC: you must provide 2 or 3 arguments')
  end if
  !
  impx=1
  drawNod = .false.
  drawMix = .false.
  zoomx = (/ 0.0, 1.0 /)
  zoomy = (/ 0.0, 1.0 /)
  typgeo=0
  10 call REDGET(indic,nitma,flott,text12,dflott)
  if (indic == 10) go to 20
  if (indic /= 3) call XABORT('G2MC: character data expected.')
  if (text12 == 'EDIT') then
    ! read the print index.
    call REDGET(indic,impx,flott,text12,dflott)
    if (indic /= 1) call XABORT('G2MC: integer data expected.')
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
  else if (text12 == ';') then
     go to 20
  else
     call XABORT('G2MC: '//text12//' is an invalid keyword.')
  end if
  go to 10

  20 sizeB = 0   !cellules de base
  sizeP = 0   !cellules placees
  sizeSA = 0  !elements geometriques
  if (c_associated(ipGeo_1)) then
     ! copy the input geometric object
     call lcmop(ipGeo,'geom_copy',0,1,0)
     call lcmequ(ipGeo_1,ipGeo)

     !initialisation des differents tableaux
     call initializeData(dimTabCelluleBase,dimTabSegArc)

     !unfold the geometry
     call g2s_unfold(ipGeo,0)

     !pretraitement des donnees lues (remplace la partie python)
     !+completion des cellules de base et remplissage du tableau
     !des cellules placees
     call prepareData(ipGeo,sizeB,sizeP,lgMaxGig)

     !en sortie, toutes les cellules de base ont tous leurs
     !champs remplis, et le tableau des cellules placees est pret

     !eclatement des cellules
     call splitCells(sizeP,sizeSA)

     !creation de nouveaux segments aux interfaces des cellules
     !et elimination des doublons
     call addSegsAndClean(sizeSA)

     !prise en compte des conditions aux limites
     call appliBoundariConditions(ipGeo,sizeSA,nbCLP)

     !calcul des nodes delimites par les elements
     allocate(merg(dimTabCelluleBase),stat=alloc_ok)
     if (alloc_ok /= 0) call XABORT("G2MC: g2s_g2mc(1) => allocation pb(1)")
     call createNodes(sizeSA,dimTabCelluleBase,nbNode,merg)
     if (sizeSA > dimTabSegArc) call XABORT('g2s_g2mc: sizeSA overflow')
  else
     if (JENTRY(nentry) == 0) call XABORT('G2M: an existing Salomon file is expected')
     !initialisation de TabSegArc
     call SALGET(datain,4,ipSal,0,'dimensions for geometry')
     nbNode=datain(3)
     sizeSA=datain(4)
     rewind(ipSal)
     allocate(tabSegArc(sizeSA),medium(nbNode))
     call initializebCData()  
     allocate(merg(nbNode),stat=alloc_ok)
     if (alloc_ok /= 0) call XABORT("G2MC: g2s_g2mc => allocation pb")
     call generateTabSegArc(ipSal,sizeSA,nbNode,nbCLP,nbFlux,merg,name_geom,impx)
  endif
  deallocate(merg)

  !impression des segArc charges
  if (ipPs /= -1) call drawSegArc(ipPs,sizeSA,drawMix,drawNod,zoomx,zoomy)

  !creation du fichier de commande Monte-Carlo
  if (index(HENTRY(1),'.tp')/=0) then
     ! generate a Tripoli4 datafile
     call generateTripoliFile(ipMC,sizeSA,nbNode)
  else if (index(HENTRY(1),'.sp')/=0) then
     ! generate a Serpent datafile
     call generateSerpentFile(ipMC,sizeSA,nbNode)
  else
     ! generate a MCNP datafile
     call generateMCNPFile(ipMC,sizeSA,nbNode)
  end if

  !liberation de la memoire allouee
  call destroyData(sizeB,sizeP)
end subroutine G2MC
