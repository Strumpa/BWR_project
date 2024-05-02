subroutine NSSANM2(nunkn,nx,ny,ll4f,ll4x,ng,bndtl,npass,nmix,idl,kn,iqfr, &
& qfr,mat,xxx,yyy,keff,diff,sigr,chi,sigf,scat,fd,savg)
!
!-----------------------------------------------------------------------
!
!Purpose:
! Compute the ANM boundary fluxes and currents using a solution of
! one- and two-node relations in Cartesian 2D geometry.
!
!Copyright:
! Copyright (C) 2022 Ecole Polytechnique de Montreal
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version
!
!Author(s): A. Hebert
!
!Parameters: input
! nunkn   number of unknowns per energy group.
! nx      number of x-nodes in the nodal calculation.
! ny      number of y-nodes in the nodal calculation.
! ll4f    number of averaged flux unknowns.
! ll4x    number of X-directed net currents.
! ng      number of energy groups.
! bndtl   set to 'flat' or 'quadratic'.
! npass   number of transverse current iterations.
! nmix    number of material mixtures in the nodal calculation.
! idl     position of averaged fluxes in unknown vector.
! kn      element-ordered interface net current unknown list.
! iqfr    node-ordered physical albedo indices.
! qfr     albedo function information.
! mat     material mixture index in eacn node.
! xxx     Cartesian coordinates along the X axis.
! yyy     Cartesian coordinates along the Y axis.
! keff    effective multiplication facctor.
! diff    diffusion coefficients
! sigr    removal cross sections.
! chi     fission spectra.
! sigf    nu times fission cross section.
! scat    scattering cross section.
! fd      discontinuity factors.
! savg    nodal fluxes and net currents.
!
!Parameters: output
! savg    nodal fluxes, boundary fluxes and net currents.
!
!-----------------------------------------------------------------------
  !
  !----
  !  subroutine arguments
  !----
  integer,intent(in) :: nunkn,nx,ny,ll4f,ll4x,ng,npass,nmix,idl(nx,ny),kn(6,nx,ny), &
  & iqfr(6,nx,ny),mat(nx,ny)
  real,intent(in) :: qfr(6,nx,ny,ng),xxx(nx+1),yyy(ny+1),keff,diff(nmix,ng), &
  & sigr(nmix,ng),chi(nmix,ng),sigf(nmix,ng),scat(nmix,ng,ng),fd(nmix,4,ng,ng)
  real, dimension(nunkn,ng),intent(inout) :: savg
  character(len=12), intent(in) :: bndtl
  !----
  ! local and allocatable arrays
  !----
  real :: xyz(4)
  real, allocatable, dimension(:) :: work1,work2,work4,work5
  real, allocatable, dimension(:,:) :: A,Lambda,work3
  real(kind=8), allocatable, dimension(:,:) :: LLR
  real(kind=8), allocatable, dimension(:,:,:,:) :: Lx,Rx,Ly,Ry
  !----
  ! scratch storage allocation
  !----
  allocate(A(ng,ng+1),Lambda(ng,ng),LLR(ng,8*ng))
  allocate(work1(ng),work2(ng),work3(ng,ng),work4(ng),work5(ng))
  allocate(Lx(ng,8*ng,nx,ny),Rx(ng,8*ng,nx,ny))
  allocate(Ly(ng,8*ng,nx,ny),Ry(ng,8*ng,nx,ny))
  !----
  !  compute 2D ANM coupling matrices for each single node
  !----
  do j=1,ny
    dely=yyy(j+1)-yyy(j)
    do i=1,nx
      delx=xxx(i+1)-xxx(i)
      ibm=mat(i,j)
      if(ibm == 0) cycle
      work1(:ng)=diff(ibm,:ng)
      work2(:ng)=sigr(ibm,:ng)
      work3(:ng,:ng)=scat(ibm,:ng,:ng)
      work4(:ng)=chi(ibm,:ng)
      work5(:ng)=sigf(ibm,:ng)
      !
      kk1=iqfr(1,i,j)
      kk2=iqfr(2,i,j)
      xyz(2:3)=xxx(i:i+1)-xxx(i)
      if(kk1 == -2) then
        ! reflection boundary condition
        xyz(1)=2.0*xyz(2)-xyz(3)
      else if(kk1 < 0) then
        ! zero/void/albedo boundary condition
        xyz(1)=-99999.
      else
        ! left neighbour
        xyz(1)=xxx(i-1)-xxx(i)
      endif
      if(kk2 == -2) then
        ! reflection boundary condition
        xyz(4)=2.0*xyz(3)-xyz(2)
      else if(kk2 < 0) then
        ! zero/void/albedo boundary condition
        xyz(4)=-99999.
      else
        ! right neighbour
        xyz(4)=xxx(i+2)-xxx(i)
      endif
      call NSSLR2(keff,ng,bndtl,xyz,dely,work1,work2,work3,work4,work5,Lx(1,1,i,j),Rx(1,1,i,j))
      !
      kk3=iqfr(3,i,j)
      kk4=iqfr(4,i,j)
      xyz(2:3)=yyy(j:j+1)-yyy(j)
      if(kk3 == -2) then
        ! reflection boundary condition
        xyz(1)=2.0*xyz(2)-xyz(3)
      else if(kk3 < 0) then
        ! zero/void/albedo boundary condition
        xyz(1)=-99999.
      else
        ! left neighbour
        xyz(1)=yyy(j-1)-yyy(j)
      endif
      if(kk4 == -2) then
        ! reflection boundary condition
        xyz(4)=2.0*xyz(3)-xyz(2)
      else if(kk4 < 0) then
        ! zero/void/albedo boundary condition
        xyz(4)=-99999.0
      else
        ! right neighbour
        xyz(4)=yyy(j+2)-yyy(j)
      endif
      call NSSLR2(keff,ng,bndtl,xyz,delx,work1,work2,work3,work4,work5, &
      & Ly(1,1,i,j),Ry(1,1,i,j))
    enddo
  enddo
  !----
  !  perform transverse current iterations
  !----
  do ipass=1,npass
    !----
    !  one- and two-node relations along X axis
    !----
    do j=1,ny
      nxmin=1
      do i=1,nx
        if(mat(i,j) > 0) exit
        nxmin=i+1
      enddo
      if(nxmin > nx) cycle
      nxmax=nx
      do i=nx,1,-1
        if(mat(i,j) > 0) exit
        nxmax=i-1
      enddo
      ! one-node relation at left
      ind1=idl(nxmin,j)
      if(ind1 == 0) call XABORT('NSSANM2: invalid idl index.(1)')
      iqf1=iqfr(1,nxmin,j)
      jxm=kn(1,nxmin,j) ; jxp=kn(2,nxmin,j) ; jym=kn(3,nxmin,j) ; jyp=kn(4,nxmin,j)
      jym_p=0 ; jyp_p=0
      if(nxmin < nx) then
        jym_p=kn(3,nxmin+1,j) ; jyp_p=kn(4,nxmin+1,j)
      endif
      A(:ng,:ng+1)=0.0
      LLR(:ng,:8*ng)=0.0
      if((iqf1 > 0).or.(iqf1 == -1)) then
        ! physical albedo
        Lambda(:ng,:ng)=0.0
        do ig=1,ng
          Lambda(ig,ig)=qfr(1,nxmin,j,ig)
        enddo
        LLR(:ng,:8*ng)=matmul(Lambda(:ng,:ng),Lx(:ng,:8*ng,nxmin,j))
        A(:ng,:ng)=-real(LLR(:ng,ng+1:2*ng),4)
        do ig=1,ng
          A(ig,ig)=-1.0+A(ig,ig)
        enddo
      else if(iqf1 == -2) then
        ! zero net current
        do ig=1,ng
          A(ig,ig)=1.0
        enddo
      else if(iqf1 == -3) then
        ! zero flux
        LLR(:ng,:8*ng)=Lx(:ng,:8*ng,nxmin,j)
        A(:ng,:ng)=real(-LLR(:ng,ng+1:2*ng),4)
      else if(iqf1 == -4) then
        call XABORT('NSSANM2: SYME boundary condition is not supported.(1)')
      else
        call XABORT('NSSANM2: illegal left X-boundary condition.')
      endif
      if(iqf1 /= -2) then
        A(:ng,ng+1)=real(matmul(LLR(:ng,:ng),savg(ind1,:ng)),4)
        do ig=1,ng
          do jg=1,ng
            if(jym /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,3*ng+jg)*savg(5*ll4f+ll4x+jym,jg),4)
            if(jym_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,4*ng+jg)*savg(5*ll4f+ll4x+jym_p,jg),4)
            if(jyp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,6*ng+jg)*savg(5*ll4f+ll4x+jyp,jg),4)
            if(jyp_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,7*ng+jg)*savg(5*ll4f+ll4x+jyp_p,jg),4)
          enddo
        enddo
      endif
      call ALSB(ng,1,A,ier,ng)
      if(ier /= 0) call XABORT('NSSANM2: singular matrix.(1)')
      savg(5*ll4f+jxm,:ng)=A(:ng,ng+1)
      !
      ! two-node relations
      do i=nxmin,nxmax-1
        ind1=idl(i,j)
        if(ind1 == 0) call XABORT('NSSANM2: invalid idl index.(2)')
        ind2=idl(i+1,j)
        if(ind2 == 0) call XABORT('NSSANM2: invalid idl index.(3)')
        if(kn(1,i+1,j) /= kn(2,i,j)) call XABORT('NSSANM2: invalid kn index.(1)')
        if(iqfr(2,i,j) /= 0) call XABORT('NSSANM2: invalid iqfr index.(1)')
        if(iqfr(1,i+1,j) /= 0) call XABORT('NSSANM2: invalid iqfr index.(2)')
        jxm=kn(1,i,j) ; jxp=kn(2,i,j) ; jym=kn(3,i,j) ; jyp=kn(4,i,j)
        jym_m=0 ; jyp_m=0 ; jym_pp=0 ; jyp_pp=0
        if((i == 1).and.(iqfr(1,1,j) == -2)) then
          jym_m=kn(3,1,j) ; jyp_m=kn(4,1,j)
        else if(i > 1) then
          jym_m=kn(3,i-1,j) ; jyp_m=kn(4,i-1,j)
        endif
        jym_p=kn(3,i+1,j) ; jyp_p=kn(4,i+1,j)
        if((i == nx-1).and.(iqfr(2,nx,j) == -2)) then
          jym_pp=kn(3,nx,j) ; jyp_pp=kn(4,nx,j)
        else if(i < nx-1) then
          jym_pp=kn(3,i+2,j) ; jyp_pp=kn(4,i+2,j)
        endif
        !
        A(:ng,:ng+1)=0.0
        ! node i
        LLR(:ng,:8*ng)=matmul(fd(mat(i,j),2,:ng,:ng),Rx(:ng,:8*ng,i,j))
        do ig=1,ng
          A(:ng,ig)=A(:ng,ig)+real(LLR(:ng,ng+ig),4)
          do jg=1,ng
            A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,jg)*savg(ind1,jg),4)
            if(jym_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,2*ng+jg)*savg(5*ll4f+ll4x+jym_m,jg),4)
            if(jym /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,3*ng+jg)*savg(5*ll4f+ll4x+jym,jg),4)
            if(jym_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,4*ng+jg)*savg(5*ll4f+ll4x+jym_p,jg),4)
            if(jyp_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,5*ng+jg)*savg(5*ll4f+ll4x+jyp_m,jg),4)
            if(jyp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,6*ng+jg)*savg(5*ll4f+ll4x+jyp,jg),4)
            if(jyp_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,7*ng+jg)*savg(5*ll4f+ll4x+jyp_p,jg),4)
          enddo
        enddo
        ! node i+1
        LLR(:ng,:8*ng)=matmul(fd(mat(i+1,j),1,:ng,:ng),Lx(:ng,:8*ng,i+1,j))
        do ig=1,ng
          A(:ng,ig)=A(:ng,ig)+real(-LLR(:ng,ng+ig),4)
          do jg=1,ng
            A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,jg)*savg(ind2,jg),4)
            if(jym /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,2*ng+jg)*savg(5*ll4f+ll4x+jym,jg),4)
            if(jym_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,3*ng+jg)*savg(5*ll4f+ll4x+jym_p,jg),4)
            if(jym_pp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,4*ng+jg)*savg(5*ll4f+ll4x+jym_pp,jg),4)
            if(jyp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,5*ng+jg)*savg(5*ll4f+ll4x+jyp,jg),4)
            if(jyp_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,6*ng+jg)*savg(5*ll4f+ll4x+jyp_p,jg),4)
            if(jyp_pp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,7*ng+jg)*savg(5*ll4f+ll4x+jyp_pp,jg),4)
          enddo
        enddo
        call ALSB(ng,1,A,ier,ng)
        if(ier /= 0) call XABORT('NSSANM2: singular matrix.(2)')
        if(jxp /= 0) savg(5*ll4f+jxp,:ng)=A(:ng,ng+1)
      enddo
      !
      ! one-node relation at right
      ind1=idl(nxmax,j)
      if(ind1 == 0) call XABORT('NSSANM2: invalid idl index.(4)')
      iqf2=iqfr(2,nxmax,j)
      jxm=kn(1,nxmax,j) ; jxp=kn(2,nxmax,j) ; jym=kn(3,nxmax,j) ; jyp=kn(4,nxmax,j)
      jym_m=0 ; jyp_m=0
      if(nxmax > 1) then
        jym_m=kn(3,nxmax-1,j) ; jyp_m=kn(4,nxmax-1,j)
      endif
      A(:ng,:ng+1)=0.0
      LLR(:ng,:8*ng)=0.0
      if((iqf2 > 0).or.(iqf2 == -1)) then
        ! physical albedo
        Lambda(:ng,:ng)=0.0
        do ig=1,ng
          Lambda(ig,ig)=qfr(2,nxmax,j,ig)
        enddo
        LLR(:ng,:8*ng)=matmul(Lambda(:ng,:ng),Rx(:ng,:8*ng,nxmax,j))
        A(:ng,:ng)=real(LLR(:ng,ng+1:2*ng),4)
        do ig=1,ng
          A(ig,ig)=-1.0+A(ig,ig)
        enddo
      else if(iqf2 == -2) then
        ! zero net current
        do ig=1,ng
          A(ig,ig)=1.0
        enddo
      else if(iqf2 == -3) then
        ! zero flux
        LLR(:ng,:8*ng)=Rx(:ng,:8*ng,nxmax,j)
        A(:ng,:ng)=real(LLR(:ng,ng+1:2*ng),4)
      else if(iqf2 == -4) then
        call XABORT('NSSANM2: SYME boundary condition is not supported.(2)')
      else
        call XABORT('NSSANM2: illegal right X-boundary condition.')
      endif
      if(iqf2 /= -2) then
        A(:ng,ng+1)=real(matmul(-LLR(:ng,:ng),savg(ind1,:ng)),4)
        do ig=1,ng
          do jg=1,ng
            if(jym_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,2*ng+jg)*savg(5*ll4f+ll4x+jym_m,jg),4)
            if(jym /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,3*ng+jg)*savg(5*ll4f+ll4x+jym,jg),4)
            if(jyp_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,5*ng+jg)*savg(5*ll4f+ll4x+jyp_m,jg),4)
            if(jyp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,6*ng+jg)*savg(5*ll4f+ll4x+jyp,jg),4)
          enddo
        enddo
      endif
      call ALSB(ng,1,A,ier,ng)
      if(ier /= 0) call XABORT('NSSANM2: singular matrix.(3)')
      if(jxp /= 0) savg(5*ll4f+jxp,:ng)=A(:ng,ng+1)
    enddo
    !----
    !  one- and two-node relations along Y axis
    !----
    do i=1,nx
      nymin=1
      do j=1,ny
        if(mat(i,j) > 0) exit
        nymin=j+1
      enddo
      if(nymin > ny) cycle
      nymax=ny
      do j=ny,1,-1
        if(mat(i,j) > 0) exit
        nymax=j-1
      enddo
      ! one-node relation at left
      ind1=idl(i,nymin)
      if(ind1 == 0) call XABORT('NSSANM2: invalid idl index.(5)')
      iqf3=iqfr(3,i,nymin)
      jxm=kn(1,i,nymin) ; jxp=kn(2,i,nymin) ; jym=kn(3,i,nymin) ; jyp=kn(4,i,nymin)
      jxm_p=0 ; jxp_p=0
      if(nymin < ny) then
        jxm_p=kn(1,i,nymin+1) ; jxp_p=kn(2,i,nymin+1)
      endif
      A(:ng,:ng+1)=0.0
      LLR(:ng,:8*ng)=0.0
      if((iqf3 > 0).or.(iqf3 == -1)) then
        ! physical albedo
        Lambda(:ng,:ng)=0.0
        do ig=1,ng
          Lambda(ig,ig)=qfr(3,i,nymin,ig)
        enddo
        LLR(:ng,:8*ng)=matmul(Lambda(:ng,:ng),Ly(:ng,:8*ng,i,nymin))
        A(:ng,:ng)=-real(LLR(:ng,ng+1:2*ng),4)
        do ig=1,ng
          A(ig,ig)=-1.0+A(ig,ig)
        enddo
      else if(iqf3 == -2) then
        ! zero net current
        do ig=1,ng
          A(ig,ig)=1.0
        enddo
      else if(iqf3 == -3) then
      ! zero flux
        LLR(:ng,:8*ng)=Ly(:ng,:8*ng,i,nymin)
        A(:ng,:ng)=real(-LLR(:ng,ng+1:2*ng),4)
      else if(iqf3 == -4) then
        call XABORT('NSSANM2: SYME boundary condition is not supported.(3)')
      else
        call XABORT('NSSANM2: illegal left Y-boundary condition.')
      endif
      if(iqf3 /= -2) then
        A(:ng,ng+1)=real(matmul(LLR(:ng,:ng),savg(ind1,:ng)),4)
        do ig=1,ng
          do jg=1,ng
            if(jxm /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,3*ng+jg)*savg(5*ll4f+jxm,jg),4)
            if(jxm_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,4*ng+jg)*savg(5*ll4f+jxm_p,jg),4)
            if(jxp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,6*ng+jg)*savg(5*ll4f+jxp,jg),4)
            if(jxp_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,7*ng+jg)*savg(5*ll4f+jxp_p,jg),4)
          enddo
        enddo
      endif
      call ALSB(ng,1,A,ier,ng)
      if(ier /= 0) call XABORT('NSSANM2: singular matrix.(4)')
      if(jym /= 0) savg(5*ll4f+ll4x+jym,:ng)=A(:ng,ng+1)
      !
      ! two-node relations
      do j=nymin,nymax-1
        ind1=idl(i,j)
        if(ind1 == 0) call XABORT('NSSANM2: invalid idl index.(6)')
        ind2=idl(i,j+1)
        if(ind2 == 0) call XABORT('NSSANM2: invalid idl index.(7)')
        if(kn(3,i,j+1) /= kn(4,i,j)) call XABORT('NSSANM2: invalid kn index.(2)')
        if(iqfr(4,i,j) /= 0) call XABORT('NSSANM2: invalid iqfr index.(3)')
        if(iqfr(3,i,j+1) /= 0) call XABORT('NSSANM2: invalid iqfr index.(4)')
        jxm=kn(1,i,j) ; jxp=kn(2,i,j) ; jym=kn(3,i,j) ; jyp=kn(4,i,j)
        jxm_m=0 ; jxp_m=0 ; jxm_pp=0 ; jxp_pp=0
        if((j == 1).and.(iqfr(3,i,1) == -2)) then
          jxm_m=kn(1,i,1) ; jxp_m=kn(2,i,1)
        else if(j > 1) then
          jxm_m=kn(1,i,j-1) ; jxp_m=kn(2,i,j-1)
        endif
        jxm_p=kn(1,i,j+1) ; jxp_p=kn(2,i,j+1)
        if((j == ny-1).and.(iqfr(4,i,ny) == -2)) then
          jxm_pp=kn(1,i,ny) ; jxp_pp=kn(2,i,ny)
        else if(j < ny-1) then
          jxm_pp=kn(1,i,j+2) ; jxp_pp=kn(2,i,j+2)
        endif
        !
        A(:ng,:ng+1)=0.0
        ! node j
        LLR(:ng,:8*ng)=matmul(fd(mat(i,j),4,:ng,:ng),Ry(:ng,:8*ng,i,j))
        do ig=1,ng
          A(:ng,ig)=A(:ng,ig)+real(LLR(:ng,ng+ig),4)
          do jg=1,ng
            A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,jg)*savg(ind1,jg),4)
            if(jxm_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,2*ng+jg)*savg(5*ll4f+jxm_m,jg),4)
            if(jxm /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,3*ng+jg)*savg(5*ll4f+jxm,jg),4)
            if(jxm_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,4*ng+jg)*savg(5*ll4f+jxm_p,jg),4)
            if(jxp_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,5*ng+jg)*savg(5*ll4f+jxp_m,jg),4)
            if(jxp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,6*ng+jg)*savg(5*ll4f+jxp,jg),4)
            if(jxp_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,7*ng+jg)*savg(5*ll4f+jxp_p,jg),4)
          enddo
        enddo
        ! node j+1
        LLR(:ng,:8*ng)=matmul(fd(mat(i,j+1),3,:ng,:ng),Ly(:ng,:8*ng,i,j+1))
        do ig=1,ng
          A(:ng,ig)=A(:ng,ig)+real(-LLR(:ng,ng+ig),4)
          do jg=1,ng
            A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,jg)*savg(ind2,jg),4)
            if(jxm /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,2*ng+jg)*savg(5*ll4f+jxm,jg),4)
            if(jxm_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,3*ng+jg)*savg(5*ll4f+jxm_p,jg),4)
            if(jxm_pp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,4*ng+jg)*savg(5*ll4f+jxm_pp,jg),4)
            if(jxp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,5*ng+jg)*savg(5*ll4f+jxp,jg),4)
            if(jxp_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,6*ng+jg)*savg(5*ll4f+jxp_p,jg),4)
            if(jxp_pp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,7*ng+jg)*savg(5*ll4f+jxp_pp,jg),4)
          enddo
        enddo
        call ALSB(ng,1,A,ier,ng)
        if(ier /= 0) call XABORT('NSSANM2: singular matrix.(5)')
        if(jyp /= 0) savg(5*ll4f+ll4x+jyp,:ng)=A(:ng,ng+1)
      enddo
      !
      ! one-node relation at right
      ind1=idl(i,nymax)
      if(ind1 == 0) call XABORT('NSSANM2: invalid idl index.(8)')
      iqf4=iqfr(4,i,nymax)
      jxm=kn(1,i,nymax) ; jxp=kn(2,i,nymax) ; jym=kn(3,i,nymax) ; jyp=kn(4,i,nymax)
      jxm_m=0 ; jxp_m=0
      if(nymax > 1) then
        jxm_m=kn(1,i,nymax-1) ; jxp_m=kn(2,i,nymax-1)
      endif
      A(:ng,:ng+1)=0.0
      LLR(:ng,:8*ng)=0.0
      if((iqf4 > 0).or.(iqf4 == -1)) then
        ! physical albedo
        Lambda(:ng,:ng)=0.0
        do ig=1,ng
          Lambda(ig,ig)=qfr(4,i,nymax,ig)
        enddo
        LLR(:ng,:8*ng)=matmul(Lambda(:ng,:ng),Ry(:ng,:8*ng,i,nymax))
        A(:ng,:ng)=real(LLR(:ng,ng+1:2*ng),4)
        do ig=1,ng
          A(ig,ig)=-1.0+A(ig,ig)
        enddo
      else if(iqf4 == -2) then
        ! zero net current
        do ig=1,ng
          A(ig,ig)=1.0
        enddo
      else if(iqf4 == -3) then
        ! zero flux
        LLR(:ng,:8*ng)=Ry(:ng,:8*ng,i,nymax)
        A(:ng,:ng)=real(LLR(:ng,ng+1:2*ng),4)
      else if(iqf4 == -4) then
        call XABORT('NSSANM2: SYME boundary condition is not supported.(4)')
      else
        call XABORT('NSSANM2: illegal right Y-boundary condition.')
      endif
      if(iqf4 /= -2) then
        A(:ng,ng+1)=real(matmul(-LLR(:ng,:ng),savg(ind1,:ng)),4)
        do ig=1,ng
          do jg=1,ng
            if(jxm_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,2*ng+jg)*savg(5*ll4f+jxm_m,jg),4)
            if(jxm /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,3*ng+jg)*savg(5*ll4f+jxm,jg),4)
            if(jxp_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,5*ng+jg)*savg(5*ll4f+jxp_m,jg),4)
            if(jxp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,6*ng+jg)*savg(5*ll4f+jxp,jg),4)
          enddo
        enddo
      endif
      call ALSB(ng,1,A,ier,ng)
      if(ier /= 0) call XABORT('NSSANM2: singular matrix.(6)')
      if(jyp /= 0) savg(5*ll4f+ll4x+jyp,:ng)=A(:ng,ng+1)
    enddo
    !----
    !  end of transverse current iterations
    !----
  enddo
  !----
  ! compute boundary fluxes
  !----
  do j=1,ny
    do i=1,nx
      ind1=idl(i,j)
      if(ind1 == 0) cycle
      jxm=kn(1,i,j) ; jxp=kn(2,i,j) ; jym=kn(3,i,j) ; jyp=kn(4,i,j)
      jym_m=0 ; jyp_m=0 ; jym_p=0 ; jyp_p=0
      if((i == 1).and.(iqfr(1,1,j) == -2)) then
        jym_m=kn(3,1,j) ; jyp_m=kn(4,1,j)
      else if(i > 1) then
        jym_m=kn(3,i-1,j) ; jyp_m=kn(4,i-1,j)
      endif
      if((i == nx).and.(iqfr(2,nx,j) == -2)) then
        jym_p=kn(3,nx,j) ; jyp_p=kn(4,nx,j)
      else if(i < nx) then
        jym_p=kn(3,i+1,j) ; jyp_p=kn(4,i+1,j)
      endif
      ! x- relations
      savg(ll4f+ind1,:ng)=real(matmul(Lx(:ng,:ng,i,j),savg(ind1,:ng)),4)
      if(jxm /= 0) savg(ll4f+ind1,:ng)=savg(ll4f+ind1,:ng)+ &
        & real(matmul(Lx(:ng,ng+1:2*ng,i,j),savg(5*ll4f+jxm,:ng)),4)
      if(jym_m /= 0) savg(ll4f+ind1,:ng)=savg(ll4f+ind1,:ng)+ &
        & real(matmul(Lx(:ng,2*ng+1:3*ng,i,j),savg(5*ll4f+ll4x+jym_m,:ng)),4)
      if(jym /= 0) savg(ll4f+ind1,:ng)=savg(ll4f+ind1,:ng)+ &
        & real(matmul(Lx(:ng,3*ng+1:4*ng,i,j),savg(5*ll4f+ll4x+jym,:ng)),4)
      if(jym_p /= 0) savg(ll4f+ind1,:ng)=savg(ll4f+ind1,:ng)+ &
        & real(matmul(Lx(:ng,4*ng+1:5*ng,i,j),savg(5*ll4f+ll4x+jym_p,:ng)),4)
      if(jyp_m /= 0) savg(ll4f+ind1,:ng)=savg(ll4f+ind1,:ng)+ &
        & real(matmul(Lx(:ng,5*ng+1:6*ng,i,j),savg(5*ll4f+ll4x+jyp_m,:ng)),4)
      if(jyp /= 0) savg(ll4f+ind1,:ng)=savg(ll4f+ind1,:ng)+ &
        & real(matmul(Lx(:ng,6*ng+1:7*ng,i,j),savg(5*ll4f+ll4x+jyp,:ng)),4)
      if(jyp_p /= 0) savg(ll4f+ind1,:ng)=savg(ll4f+ind1,:ng)+ &
        & real(matmul(Lx(:ng,7*ng+1:8*ng,i,j),savg(5*ll4f+ll4x+jyp_p,:ng)),4)
      !
      ! x+ relations
      savg(2*ll4f+ind1,:ng)=real(matmul(Rx(:ng,:ng,i,j),savg(ind1,:ng)),4)
      if(jxp /= 0) savg(2*ll4f+ind1,:ng)=savg(2*ll4f+ind1,:ng)+ &
        & real(matmul(Rx(:ng,ng+1:2*ng,i,j),savg(5*ll4f+jxp,:ng)),4)
      if(jym_m /= 0) savg(2*ll4f+ind1,:ng)=savg(2*ll4f+ind1,:ng)+ &
        & real(matmul(Rx(:ng,2*ng+1:3*ng,i,j),savg(5*ll4f+ll4x+jym_m,:ng)),4)
      if(jym /= 0) savg(2*ll4f+ind1,:ng)=savg(2*ll4f+ind1,:ng)+ &
        & real(matmul(Rx(:ng,3*ng+1:4*ng,i,j),savg(5*ll4f+ll4x+jym,:ng)),4)
      if(jym_p /= 0) savg(2*ll4f+ind1,:ng)=savg(2*ll4f+ind1,:ng)+ &
        & real(matmul(Rx(:ng,4*ng+1:5*ng,i,j),savg(5*ll4f+ll4x+jym_p,:ng)),4)
      if(jyp_m /= 0) savg(2*ll4f+ind1,:ng)=savg(2*ll4f+ind1,:ng)+ &
        & real(matmul(Rx(:ng,5*ng+1:6*ng,i,j),savg(5*ll4f+ll4x+jyp_m,:ng)),4)
      if(jyp /= 0) savg(2*ll4f+ind1,:ng)=savg(2*ll4f+ind1,:ng)+ &
        & real(matmul(Rx(:ng,6*ng+1:7*ng,i,j),savg(5*ll4f+ll4x+jyp,:ng)),4)
      if(jyp_p /= 0) savg(2*ll4f+ind1,:ng)=savg(2*ll4f+ind1,:ng)+ &
        & real(matmul(Rx(:ng,7*ng+1:8*ng,i,j),savg(5*ll4f+ll4x+jyp_p,:ng)),4)
      !
      jxm_m=0 ; jxp_m=0 ; jxm_p=0 ; jxp_p=0
      jxm=kn(1,i,j) ; jxp=kn(2,i,j) ; jym=kn(3,i,j) ; jyp=kn(4,i,j)
      if((j == 1).and.(iqfr(3,i,1) == -2)) then
        jxm_m=kn(1,i,1) ; jxp_m=kn(2,i,1)
      else if(j > 1) then
        jxm_m=kn(1,i,j-1) ; jxp_m=kn(2,i,j-1)
      endif
      if((j == ny).and.(iqfr(4,i,ny) == -2)) then
        jxm_p=kn(1,i,ny) ; jxp_p=kn(2,i,ny)
      else if(j < ny) then
        jxm_p=kn(1,i,j+1) ; jxp_p=kn(2,i,j+1)
      endif
      ! y- relations
      savg(3*ll4f+ind1,:ng)=real(matmul(Ly(:ng,:ng,i,j),savg(ind1,:ng)),4)
      if(jym /= 0) savg(3*ll4f+ind1,:ng)=savg(3*ll4f+ind1,:ng)+ &
        & real(matmul(Ly(:ng,ng+1:2*ng,i,j),savg(5*ll4f+ll4x+jym,:ng)),4)
      if(jxm_m /= 0) savg(3*ll4f+ind1,:ng)=savg(3*ll4f+ind1,:ng)+ &
        & real(matmul(Ly(:ng,2*ng+1:3*ng,i,j),savg(5*ll4f+jxm_m,:ng)),4)
      if(jxm /= 0) savg(3*ll4f+ind1,:ng)=savg(3*ll4f+ind1,:ng)+ &
        & real(matmul(Ly(:ng,3*ng+1:4*ng,i,j),savg(5*ll4f+jxm,:ng)),4)
      if(jxm_p /= 0) savg(3*ll4f+ind1,:ng)=savg(3*ll4f+ind1,:ng)+ &
        & real(matmul(Ly(:ng,4*ng+1:5*ng,i,j),savg(5*ll4f+jxm_p,:ng)),4)
      if(jxp_m /= 0) savg(3*ll4f+ind1,:ng)=savg(3*ll4f+ind1,:ng)+ &
        & real(matmul(Ly(:ng,5*ng+1:6*ng,i,j),savg(5*ll4f+jxp_m,:ng)),4)
      if(jxp /= 0) savg(3*ll4f+ind1,:ng)=savg(3*ll4f+ind1,:ng)+ &
        & real(matmul(Ly(:ng,6*ng+1:7*ng,i,j),savg(5*ll4f+jxp,:ng)),4)
      if(jxp_p /= 0) savg(3*ll4f+ind1,:ng)=savg(3*ll4f+ind1,:ng)+ &
        & real(matmul(Ly(:ng,7*ng+1:8*ng,i,j),savg(5*ll4f+jxp_p,:ng)),4)
      !
      ! y+ relations
      savg(4*ll4f+ind1,:ng)=real(matmul(Ry(:ng,:ng,i,j),savg(ind1,:ng)),4)
      if(jyp /= 0) savg(4*ll4f+ind1,:ng)=savg(4*ll4f+ind1,:ng)+ &
        & real(matmul(Ry(:ng,ng+1:2*ng,i,j),savg(5*ll4f+ll4x+jyp,:ng)),4)
      if(jxm_m /= 0) savg(4*ll4f+ind1,:ng)=savg(4*ll4f+ind1,:ng)+ &
        & real(matmul(Ry(:ng,2*ng+1:3*ng,i,j),savg(5*ll4f+jxm_m,:ng)),4)
      if(jxm /= 0) savg(4*ll4f+ind1,:ng)=savg(4*ll4f+ind1,:ng)+ &
        & real(matmul(Ry(:ng,3*ng+1:4*ng,i,j),savg(5*ll4f+jxm,:ng)),4)
      if(jxm_p /= 0) savg(4*ll4f+ind1,:ng)=savg(4*ll4f+ind1,:ng)+ &
        & real(matmul(Ry(:ng,4*ng+1:5*ng,i,j),savg(5*ll4f+jxm_p,:ng)),4)
      if(jxp_m /= 0) savg(4*ll4f+ind1,:ng)=savg(4*ll4f+ind1,:ng)+ &
        & real(matmul(Ry(:ng,5*ng+1:6*ng,i,j),savg(5*ll4f+jxp_m,:ng)),4)
      if(jxp /= 0) savg(4*ll4f+ind1,:ng)=savg(4*ll4f+ind1,:ng)+ &
        & real(matmul(Ry(:ng,6*ng+1:7*ng,i,j),savg(5*ll4f+jxp,:ng)),4)
      if(jxp_p /= 0) savg(4*ll4f+ind1,:ng)=savg(4*ll4f+ind1,:ng)+ &
        & real(matmul(Ry(:ng,7*ng+1:8*ng,i,j),savg(5*ll4f+jxp_p,:ng)),4)
    enddo
  enddo
  !----
  ! scratch storage deallocation
  !----
  deallocate(Ry,Ly,Rx,Lx)
  deallocate(work5,work4,work3,work2,work1)
  deallocate(LLR,Lambda,A)
end subroutine NSSANM2
