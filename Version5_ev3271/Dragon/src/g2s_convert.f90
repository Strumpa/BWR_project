!
!-----------------------------------------------------------------------
!
!Purpose:
! Convert an Alamos surfacic file towards Salomon format.
!
!Copyright:
! Copyright (C) 2022 Ecole Polytechnique de Montreal.
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
!Author(s):
! A. Hebert
!
!Parameters: input
! ipAl   Alamos ascii file index
!
!Parameters: output
! ipSal  Salomon ascii file index
!
!-----------------------------------------------------------------------
!
subroutine g2s_convert(impx,ipAl,ipZa,ipSal)
    use constUtiles
    use segArc
    ! typgeo type of geometry (0: TISO tracking; 5: rectangle TRAN;
    !         6: rectangle REFL; 7: eight;8: SA60; 9: hexagon; 10: RA60;
    !         11: R120)
    integer,intent(in) :: impx,ipAl,ipZa,ipSal

    ! local variables
    integer :: i, j, k, nbfold, defautCl, nangles, saltype, idummy
    double precision :: aStart, bStart, rStart, anglStart, aEnd, bEnd, rEnd, anglEnd, &
                       & radius, angle, angle0, dx, dy, dx0, dy0, r_xmin, r_xmax, &
                       & r_ymin, r_ymax
    character(len=4) :: strDCL
    character(len=12) :: text12, name_geom
    character(len=40) :: text40
    character(len=131) :: hsmg
    integer, allocatable, dimension(:) :: nber,merg,milTab,milTab2,iper
    integer, allocatable, dimension(:,:) :: bcElems,iboundary
    double precision, allocatable, dimension(:) :: ccx,ccy,ttx,tty
    double precision, allocatable, dimension(:,:) :: apnodes
    character(len=24), allocatable, dimension(:) :: propertyNames,mixNames,mixNames2
    type(t_segArc) :: sa
    double precision, allocatable, dimension(:) :: list_angles
    !
    ! set angles
    nangles=-99
    if(typgeo == 0) then
      nangles=0
    else if((abs(typgeo) == 5).or.(abs(typgeo) == 6).or.(abs(typgeo) == 11)) then
      nangles=4
    else if((abs(typgeo) == 7).or.(abs(typgeo) == 8).or.(abs(typgeo) == 10)) then
      nangles=3
    else if(abs(typgeo) == 9) then
      nangles=6
    else
      write(hsmg,'(35hg2s_convert: invalid symmetry type=,i3,1h.)') typgeo
      call XABORT(hsmg)
    endif
    allocate(list_angles(nangles))
    if((abs(typgeo) == 5).or.(abs(typgeo) == 6)) then
      list_angles=(/ 0.0, 0.0, 90.0, 90.0 /)
    else if(abs(typgeo) == 7) then
      list_angles=(/ 0.0, 45.0, 90.0 /)
    else if((abs(typgeo) == 8).or.(abs(typgeo) == 10)) then
      list_angles=(/ 0.0, 60.0, 120.0 /)
    else if(abs(typgeo) == 9) then
      list_angles=(/ 0.0, 120.0, 60.0, 0.0, 120.0, 60.0 /)
    else if(abs(typgeo) == 11) then
      list_angles=(/ 0.0, 0.0, 60.0, 60.0 /)
    endif

    ! ----------------------------------
    ! Read the Alamos surfacic file
    ! ----------------------------------
    read(ipAl,'(a12,a24)',end=100) text12,name_geom
    10 read(ipAl,'(a12,a24)',end=100) text12,name_geom
    if (text12 /= ' Geometrie:') go to 10
    20 read(ipAl,'(a12)',end=100) text12
    if (text12 /= ' Milieux:') go to 20
    read(ipAl,"(i7)",advance='no') nbmil
    allocate(mixNames(nbmil))
    read(ipAl,*) mixNames(:)
    30 read(ipAl,'(a12)',end=100) text12
    if (text12 /= ' Mailles:') go to 30
    read(ipAl,"(i7)") nbNode
    allocate(milTab(nbNode),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT('g2s_convert: allocation problem.')
    read(ipAl,'(10i7)',end=100) milTab(:)
    nbFlux=nbNode
    if(ipZa /= -1) then
      ! ----------------------------------
      ! Read the PropertyMap file
      ! ----------------------------------
      allocate(mixNames2(nbNode),milTab2(nbNode),propertyNames(nbNode),iper(nbNode),stat=alloc_ok)
      if (alloc_ok /= 0) call XABORT('g2s_convert: allocation problem.')
      40 read(ipZa,'(a12,a24)',end=110) text12
      if (text12 /= 'PropertyMap:') go to 40
      read(ipZA,*) propertyNames(:)
      nbmil=0
      loop1: do i=1,nbNode
        do j=1,nbmil
          if(propertyNames(i) == mixNames2(j)) then
            milTab2(i)=j; cycle loop1
          endif
        enddo
        nbmil=nbmil+1
        mixNames2(nbmil)=propertyNames(i)
        milTab2(i)=nbmil
        iper(nbmil)=milTab(i)
      enddo loop1
      write(6,'(1x,a)') name_geom
      do i=1,nbmil
        write(6,'(1x,i5,2x,a,4h--> ,a)') i,mixNames2(i),mixNames(iper(i))
      enddo
      deallocate(mixNames); allocate(mixNames(nbmil))
      milTab(:nbNode)=milTab2(:nbNode); mixNames(:nbmil)=mixNames2(:nbmil)
      deallocate(iper,propertyNames,milTab2,mixNames2)
    else
      write(6,'(1x,a)') name_geom
      do i=1,nbmil
        write(6,'(1x,i5,2x,a)') i,mixNames(i)
      enddo
    endif
    if(impx > 1) then
      write(6,'(5h--cut,75(1h-))')
      i=1
      do j=1,1+(nbmil-1)/4
        write(text40,'(15h(10h  INTEGER ,,i1,12h(a,1x),2h:=,,i1,8hi5,2h ;))') &
        & min(i+3,nbmil)-i+1,min(i+3,nbmil)-i+1
        write(6,text40) (trim(mixNames(k)),k=i,min(i+3,nbmil)),(k,k=i,min(i+3,nbmil))
        i=i+4
      enddo
      write(6,'(5h--cut,75(1h-))')
    endif
    deallocate(mixNames)

    ! read surfacic elements
    rewind(ipAl)
    50 read(ipAl,'(a12)',end=100) text12
    if (text12 /= ' Noeuds:') go to 50
    read(ipAl,'(i7)') nbPoints
    allocate(apnodes(2,nbPoints))
    r_xmin=999.0
    r_ymin=999.0
    do i=1,nbPoints
      if((abs(typgeo) == 8).or.(abs(typgeo) == 10).or.(abs(typgeo) == 11)) then
        read(ipAl,*,end=100) apnodes(2,i),apnodes(1,i)
        apnodes(2,i)=-apnodes(2,i) ! rotate Alamos geometry by 90 degres
      else
        read(ipAl,*,end=100) apnodes(1,i),apnodes(2,i)
      endif
      r_xmin=min(r_xmin,apnodes(1,i))
      r_ymin=min(r_ymin,apnodes(2,i))
    enddo
    do i=1,nbPoints
      apnodes(1,i)=apnodes(1,i)-r_xmin
      apnodes(2,i)=apnodes(2,i)-r_ymin
    enddo
    r_xmin=999.0
    r_xmax=-999.0
    r_ymin=999.0
    r_ymax=-999.0
    do i=1,nbPoints
      r_xmin=min(r_xmin,apnodes(1,i))
      r_ymin=min(r_ymin,apnodes(2,i))
      r_xmax=max(r_xmax,apnodes(1,i))
      r_ymax=max(r_ymax,apnodes(2,i))
    enddo
    basex=real(r_xmax-r_xmin)
    basey=real(r_ymax-r_ymin)
    
    ! compute boundary elements characteristics
    if(abs(typgeo) > 0) then
      allocate(ttx(nangles),tty(nangles),ccx(nangles),ccy(nangles))
      if(abs(typgeo) == 9) basex=basex/2.0
      if(abs(typgeo) == 11) basex=2.0*basex/3.0
      do i = 1,nangles
        ccx(i)=0.0 ; ccy(i)=0.0 ; ttx(i)=0.0 ; tty(i)=0.0
        if((abs(typgeo) == 5).or.(abs(typgeo) == 6)) then
          select case(i)
            case(2)
              ccy(i)=basey
            case(4)
              ccx(i)=basex
          end select
        else if((abs(typgeo) == 7).or.(abs(typgeo) == 8).or.(abs(typgeo) == 10)) then
          select case(i)
            case(3)
              ccx(i)=basex
          end select
        else if(abs(typgeo) == 9) then
          select case(i)
            case(1:2)
              ccx(i)=basex/2.0
            case(3)
              ccy(i)=basex*sqrt(0.75)
            case(4)
              ccx(i)=basex/2.0
              ccy(i)=2.0*basex*sqrt(0.75)
            case(5)
              ccx(i)=2.0*basex
              ccy(i)=basex*sqrt(0.75)
            case(6)
              ccx(i)=1.5*basex
          end select
        else if(abs(typgeo) == 11) then
          select case(i)
            case(2)
              ccx(i)=basex/2.0
              ccy(i)=basex*sqrt(0.75)
            case(4)
              ccx(i)=basex
          end select
        endif
        if(typgeo == 5) then
          ! pure translation in Cartesian geometry
          select case(i)
            case(1)
              tty(i)=basey
            case(2)
              tty(i)=-basey
            case(3)
              ttx(i)=basex
            case(4)
              ttx(i)=-basex
          end select
        else if(typgeo == 9) then
          ! pure translation in hexagonal geometry
          select case(i)
            case(1)
              tty(i)=2.0*basex*sqrt(0.75)
            case(2)
              ttx(i)=2.0*basex
              tty(i)=basex*sqrt(0.75)
            case(3)
              ttx(i)=2.0*basex
              tty(i)=-basex*sqrt(0.75)
            case(4)
              tty(i)=-2.0*basex*sqrt(0.75)
            case(5)
              ttx(i)=-2.0*basex
              tty(i)=-basex*sqrt(0.75)
            case(6)
              ttx(i)=-2.0*basex
              tty(i)=basex*sqrt(0.75)
          end select
        else
          ttx(i)=ccx(i) ; tty(i)=ccy(i)
        endif
      end do
    endif

    ! read surfacic elements
    read(ipAl,'(a12)') text12
    if (text12 /= ' Aretes:') call XABORT('g2s_convert: keyword Aretes: expected.')
    read(ipAl,'(i7)',end=100) iszSA
    allocate(tabSegArc(iszSA),iboundary(iszSA,2))
    ibd=0
    do i=1,iszSA
      read(ipAl,*,end=100) itype
      backspace(ipAl)
      if(itype == 0) then
        ! line segment
        tabSegArc(i)%typ=1
        read(ipAl,*) itype,ipt1,ipt2,nmoins,nplus
        tabSegArc(i)%x=real(apnodes(1,ipt1))
        tabSegArc(i)%y=real(apnodes(2,ipt1))
        tabSegArc(i)%dx=real(apnodes(1,ipt2))
        tabSegArc(i)%dy=real(apnodes(2,ipt2))
      else if(itype == 1) then
        ! full circle
        tabSegArc(i)%typ=2
        read(ipAl,*) itype,ipt3,radius,nplus,nmoins
        tabSegArc(i)%x=real(apnodes(1,ipt3))
        tabSegArc(i)%y=real(apnodes(2,ipt3))
        tabSegArc(i)%r=real(radius)
        tabSegArc(i)%b=0.0 ; tabSegArc(i)%a=0.0
      else if(itype == 2) then
        ! circular arc
        tabSegArc(i)%typ=3
        read(ipAl,*) itype,ipt1,ipt2,ipt3,nplus,nmoins
        aStart=real(apnodes(1,ipt1)-apnodes(1,ipt3))
        bStart=real(apnodes(2,ipt1)-apnodes(2,ipt3))
        rStart=sqrt(aStart**2+bStart**2)
        if(bStart >= 0.0) then
          anglStart=acos(aStart/rStart)
        else
          anglStart=dpi_c-acos(aStart/rStart)
        endif
        aEnd=real(apnodes(1,ipt2)-apnodes(1,ipt3))
        bEnd=real(apnodes(2,ipt2)-apnodes(2,ipt3))
        rEnd=sqrt(aEnd**2+bEnd**2)
        if(bEnd >= 0.0) then
          anglEnd=acos(aEnd/rEnd)
        else
          anglEnd=dpi_c-acos(aEnd/rEnd)
        endif
        tabSegArc(i)%x=real(apnodes(1,ipt3))
        tabSegArc(i)%y=real(apnodes(2,ipt3))
        tabSegArc(i)%r=real(0.5d0*(rStart+rEnd))
        tabSegArc(i)%b=real(anglEnd)
        tabSegArc(i)%a=real(anglStart)
      else
        write(hsmg,'(35hg2s_convert: invalid element type=,i3,1h.)') itype
        call XABORT(hsmg)
      endif
      tabSegArc(i)%nodeg=nplus
      tabSegArc(i)%noded=nmoins
      if(abs(typgeo) > 0) then
        if(((nplus==0).or.(nmoins==0)).and.(itype == 0)) then
          ibd=ibd+1
          if(ibd > iszSA) call XABORT('g2s_convert: boundary overflow')
          iboundary(ibd,2)=0
          dx=tabSegArc(i)%x-tabSegArc(i)%dx
          dy=tabSegArc(i)%y-tabSegArc(i)%dy
          angle=atan(dy/dx)*rad2deg
          if(angle < 0.0) angle=angle+180.
          do j=1,nangles
            if(abs(angle-list_angles(j)) > 0.5) cycle
            if((abs(tabSegArc(i)%x-ccx(j)) < epsilon).and.(abs(tabSegArc(i)%y-ccy(j)) < epsilon)) then
              iboundary(ibd,2)=j
            else if((abs(tabSegArc(i)%dx-ccx(j)) < epsilon).and.(abs(tabSegArc(i)%dy-ccy(j)) < epsilon)) then
              iboundary(ibd,2)=j
            else
              dx0=tabSegArc(i)%x-ccx(j)
              dy0=tabSegArc(i)%y-ccy(j)
              angle0=atan(dy0/dx0)*rad2deg
              if(angle0 < 0.0) angle0=angle0+180.
              if(abs(angle-angle0) <= 1.0) iboundary(ibd,2)=j ! look for 1 degree agreement
            endif
          enddo
          iboundary(ibd,1)=i
        endif
      endif
    enddo
    if(abs(typgeo) > 0) deallocate(ccy,ccx)
    
    ! list boundary elements
    allocate(bcElems(iszSA,nangles),nber(nangles))
    if(abs(typgeo) > 0) then
      nber(:)=0
      bcElems(:,:)=-1
      do ibound=1,ibd
        i1=iboundary(ibound,1)
        i2=iboundary(ibound,2)
        if(i2 == 0) then
          write(hsmg,'(26hg2s_convert: boundary side,i7,18h is not allocated., &
          & 33h Try a different value of typgeo.)') i1
          call XABORT(hsmg)
        endif
        if(i2 > nangles) call XABORT('g2s_convert: boundary angle overflow.')
        nber(i2)=nber(i2)+1
        if(nber(i2) > iszSA) call XABORT('g2s_convert: boundary element overflow.')
        bcElems(nber(i2),i2)=i1
      enddo
    endif
    deallocate(iboundary)
    
    ! read Bords data in some specific Alamos files
    read(ipAl,'(a12)') text12
    if (text12 == ' Bords:') then
      read(ipAl,*) nbBords
      do ib=1,nbBords
        read(ipAl,*) itype,ib1,ib2
      enddo
      read(ipAl,'(a12)') text12
    endif
    
    ! skip milTab data
    if (text12 /= ' Mailles:') call XABORT('g2s_convert: keyword Mailles: expected.')
    read(ipAl,'(i7)') nbNode
    read(ipAl,'(10i7)',end=100) (idummy, i=1,nbNode)
    read(ipAl,'(a12)') text12
    if (text12 == ' ') read(ipAl,'(a12)') text12
      print *,'read=',text12
    if (text12 /= ' Fin:') call XABORT('g2s_convert: keyword Fin: expected.')
    !
    ! set nbfold
    nbfold=0
    if(typgeo == -5) then
      typgeo=1 ; nbfold=4
    else if(typgeo == -6) then
      typgeo=1 ; nbfold=4
    else if(typgeo == -7) then
      typgeo=1 ; nbfold=8
    else if(typgeo == -8) then
      typgeo=1 ; nbfold=6
    else if(typgeo == -9) then
      typgeo=1 ; nbfold=1
    else if(typgeo == -10) then
      typgeo=2 ; nbfold=6
    else if(typgeo == -11) then
      typgeo=2 ; nbfold=3
    endif
    ! ----------------------------------
    ! Generate the Salomon surfacic file
    ! ----------------------------------
    write(ipSal,'(a5)') 'BEGIN'

    write(ipSal,'(/a24)') 'DEFINE DOMAINE'
    write(ipSal,'(a24/)') '=============='

    write(ipSal,'(a28/)') '1.main dimensions:'
    write(ipSal,'(a32)')  '*typgeo  nbfold  nbNode  nbelem '
    write(ipSal,'(10'//formati//')') typgeo,nbfold,nbNode,iszSA,1,nbFlux

    write(ipSal,'(/a37)') '2.impression and precision:'
    write(ipSal,'(/a21)') '*index   kndex   prec'
    write(ipSal,'(10'//formati//')') 0,0,1

    write(ipSal,'(/a40)') '3.precision of geometry data:'
    write(ipSal,'(/a4)')  '*eps'
    write(ipSal,'('//formatr//')') gSALeps

    write(ipSal,'(/a58)') '4.flux region number per geometry region (mesh):'
    write(ipSal,'(/a6)')  '*merge'
    allocate(merg(nbNode))
    do i=1,nbNode
      merg(i)=i
    enddo
    write(ipSal,'(10'//formati//')') merg(:nbNode)
    deallocate(merg)

    write(ipSal,'(/a29)') '5.name of geometry:'
    write(ipSal,'(/a11)')  '*macro_name'
    write(ipSal,'(4(3x,a12))') name_geom

    write(ipSal,'(/a46)') '6.macro order number per flux region:'
    write(ipSal,'(/a14)')  '*macro_indices'
    write(ipSal,'(10'//formati//')') (1,i=1,nbFlux)
    
    write(ipSal,'(/a57)') '7.read integer and real data for each elements:'
    do i=1,iszSA
      sa = tabSegArc(i)
      cx=real(sa%x)
      cy=real(sa%y)
      if(sa%typ == 1) then
        ! line segment
        tx= real(sa%dx-sa%x)
        ty= real(sa%dy-sa%y)
        delta=0.
      else if(sa%typ == 2) then
        ! full circle
        tx=real(sa%r)
        ty=0.
        delta=0.
      else if(sa%typ == 3) then
        ! circular arc
        tx=real(sa%r)
        ty=real(sa%a*rad2deg)
        if (sa%b>sa%a) then
          delta=real((sa%b-sa%a)*rad2deg)
        else
          delta=real((sa%b-sa%a)*rad2deg+360.d0)
        endif
      endif
      write(ipSal,*)
      write(ipSal,*) 'elem =',i 
      write(ipSal,'(a22)')    '*type    node-   node+' 
      write(ipSal,'(10'//formati//')') sa%typ,sa%nodeg,sa%noded
      write(ipSal,'(a63)')    '*cx            cy            ex or R       &
            &ey or theta1  theta2' 
      write(ipSal,'(5'//formatr//')')  cx,cy,tx,ty,delta
    enddo
    deallocate(tabSegArc)

    ! write boundary elements lists
    write(ipSal,'(/a63)') '8.read integer and real data for boundary conditions:'
    defautCl = 0
    albedo = 0.0
    strDCL = 'VOID'
    write(ipSal,'(/a40)') '*defaul  nbbcda  allsur  divsur  ndivsur'
    write(ipSal,'(10'//formati//')') defautCl,nangles,0,0,0
    write(ipSal,'(/a24)') 'DEFAULT = ' // strDCL
    write(ipSal,'(a24)')  '=============='
    write(ipSal,'(/a17)') '*albedo  deltasur'
    write(ipSal,'(5'//formatr//')') albedo,0.0
    if(typgeo /= 0) then
      do i = 1,nangles
        write(ipSal,*)
        write(ipSal,*) 'particular boundary condition number',i
        write(ipSal,'(/a13)') '*type    nber'
        if(typgeo <= 2) then
          saltype = 1 ! isotropic reflexion
        else if((typgeo == 5).or.(typgeo == 9)) then
          saltype = 2 ! translation
        else if(typgeo < 9) then
          saltype = 4 ! axial symmetry (specular reflexion)
        else if(typgeo == 10) then ! RA60
          if(i == 1) then
            saltype = 2 ! translation
          else
            saltype = 3 ! rotation
          endif
        else if(typgeo == 11) then ! R120
          if((i == 1).or.(i == 4)) then
            saltype = 2 ! translation
          else
            saltype = 3 ! rotation
          endif
        else
          call XABORT('g2s_convert: unknown type of geometry.')
        endif
        write(ipSal,'(10'//formati//')') saltype, nber(i)
        write(ipSal,'(/a14)') '*elems(1,nber)'
        write(ipSal,'(10'//formati//')') bcElems(:nber(i),i)
        write(ipSal,'(/a22)') '*cx      cy      angle'
        write(ipSal,'(5'//formatr//')') ttx(i),tty(i), list_angles(i)
      end do
      deallocate(tty,ttx)
    endif
    deallocate(nber,bcElems,list_angles)

    write(ipSal,'(/a28)') '9.medium per node:'
    write(ipSal,'(/a11)') '*mil(nbreg)'
    write(ipSal,'(10'//formati//')') milTab(:nbNode)
    deallocate(milTab)

    write(ipSal,'(/a3)') 'END'
    rewind(ipSal)
    return
    !
    100 call XABORT('g2s_convert: end of Alamos surfacic file encountered.')
    110 call XABORT('g2s_convert: end of PropertyMap file encountered.')
 end subroutine g2s_convert
