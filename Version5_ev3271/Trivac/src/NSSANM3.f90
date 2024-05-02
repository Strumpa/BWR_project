subroutine NSSANM3(nunkn,nx,ny,nz,ll4f,ll4x,ll4y,ng,bndtl,npass,nmix,idl, &
& kn,iqfr,qfr,mat,xxx,yyy,zzz,keff,diff,sigr,chi,sigf,scat,fd,savg)
!
!-----------------------------------------------------------------------
!
!Purpose:
! Compute the ANM boundary fluxes and currents using a solution of
! one- and two-node relations in Cartesian 3D geometry.
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
! nz      number of z-nodes in the nodal calculation.
! ll4f    number of averaged flux unknowns.
! ll4x    number of X-directed net currents.
! ll4y    number of Y-directed net currents.
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
! zzz     Cartesian coordinates along the Z axis.
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
  integer,intent(in) :: nunkn,nx,ny,nz,ll4f,ll4x,ll4y,ng,npass,nmix,idl(nx,ny,nz), &
  & kn(6,nx,ny,nz),iqfr(6,nx,ny,nz),mat(nx,ny,nz)
  real,intent(in) :: qfr(6,nx,ny,nz,ng),xxx(nx+1),yyy(ny+1),zzz(nz+1),keff,diff(nmix,ng), &
  & sigr(nmix,ng),chi(nmix,ng),sigf(nmix,ng),scat(nmix,ng,ng),fd(nmix,6,ng,ng)
  real, dimension(nunkn,ng),intent(inout) :: savg
  character(len=12), intent(in) :: bndtl
  !----
  ! local and allocatable arrays
  !----
  real :: xyz(4)
  real, allocatable, dimension(:) :: work1,work2,work4,work5
  real, allocatable, dimension(:,:) :: A,Lambda,work3
  real(kind=8), allocatable, dimension(:,:) :: LLR
  real(kind=8), allocatable, dimension(:,:,:,:,:) :: Lx,Rx,Ly,Ry,Lz,Rz
  !----
  ! scratch storage allocation
  !----
  allocate(A(ng,ng+1),Lambda(ng,ng),LLR(ng,14*ng))
  allocate(work1(ng),work2(ng),work3(ng,ng),work4(ng),work5(ng))
  allocate(Lx(ng,14*ng,nx,ny,nz),Rx(ng,14*ng,nx,ny,nz))
  allocate(Ly(ng,14*ng,nx,ny,nz),Ry(ng,14*ng,nx,ny,nz))
  allocate(Lz(ng,14*ng,nx,ny,nz),Rz(ng,14*ng,nx,ny,nz))
  !----
  !  compute 3D ANM coupling matrices for each single node
  !----
  do k=1,nz
    delz=zzz(k+1)-zzz(k)
    do j=1,ny
      dely=yyy(j+1)-yyy(j)
      do i=1,nx
        delx=xxx(i+1)-xxx(i)
        ibm=mat(i,j,k)
        if(ibm == 0) cycle
        work1(:ng)=diff(ibm,:ng)
        work2(:ng)=sigr(ibm,:ng)
        work3(:ng,:ng)=scat(ibm,:ng,:ng)
        work4(:ng)=chi(ibm,:ng)
        work5(:ng)=sigf(ibm,:ng)
        !
        kk1=iqfr(1,i,j,k)
        kk2=iqfr(2,i,j,k)
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
        call NSSLR3(keff,ng,bndtl,xyz,dely,delz,work1,work2,work3, &
        & work4,work5,Lx(1,1,i,j,k),Rx(1,1,i,j,k))
        !
        kk3=iqfr(3,i,j,k)
        kk4=iqfr(4,i,j,k)
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
        call NSSLR3(keff,ng,bndtl,xyz,delz,delx,work1,work2,work3, &
        & work4,work5,Ly(1,1,i,j,k),Ry(1,1,i,j,k))
        !
        kk5=iqfr(5,i,j,k)
        kk6=iqfr(6,i,j,k)
        xyz(2:3)=zzz(k:k+1)-zzz(k)
        if(kk5 == -2) then
          ! reflection boundary condition
          xyz(1)=2.0*xyz(2)-xyz(3)
        else if(kk5 < 0) then
          ! zero/void/albedo boundary condition
          xyz(1)=-99999.
        else
          ! left neighbour
          xyz(1)=zzz(k-1)-zzz(k)
        endif
        if(kk6 == -2) then
          ! reflection boundary condition
          xyz(4)=2.0*xyz(3)-xyz(2)
        else if(kk6 < 0) then
          ! zero/void/albedo boundary condition
          xyz(4)=-99999.0
        else
          ! right neighbour
          xyz(4)=zzz(k+2)-zzz(k)
        endif
        call NSSLR3(keff,ng,bndtl,xyz,delx,dely,work1,work2,work3, &
        & work4,work5,Lz(1,1,i,j,k),Rz(1,1,i,j,k))
      enddo
    enddo
  enddo
  !----
  !  perform transverse current iterations
  !----
  do ipass=1,npass
    !----
    !  one- and two-node relations along X axis
    !----
    do k=1,nz
      do j=1,ny
        nxmin=1
        do i=1,nx
          if(mat(i,j,k) > 0) exit
          nxmin=i+1
        enddo
        if(nxmin > nx) cycle
        nxmax=nx
        do i=nx,1,-1
          if(mat(i,j,k) > 0) exit
          nxmax=i-1
        enddo
        ! one-node relation at left
        ind1=idl(nxmin,j,k)
        if(ind1 == 0) call XABORT('NSSANM3: invalid idl index.(1)')
        iqf1=iqfr(1,nxmin,j,k)
        jxm=kn(1,nxmin,j,k) ; jxp=kn(2,nxmin,j,k) ; jym=kn(3,nxmin,j,k) ; jyp=kn(4,nxmin,j,k)
        jzm=kn(5,nxmin,j,k) ; jzp=kn(6,nxmin,j,k)
        jym_p=0 ; jyp_p=0 ; jzm_p=0 ; jzp_p=0
        if(nxmin < nx) then
          jym_p=kn(3,nxmin+1,j,k) ; jyp_p=kn(4,nxmin+1,j,k)
          jzm_p=kn(5,nxmin+1,j,k) ; jzp_p=kn(6,nxmin+1,j,k)
        endif
        A(:ng,:ng+1)=0.0
        LLR(:ng,:14*ng)=0.0
        if((iqf1 > 0).or.(iqf1 == -1)) then
          ! physical albedo
          Lambda(:ng,:ng)=0.0
          do ig=1,ng
            Lambda(ig,ig)=qfr(1,nxmin,j,k,ig)
          enddo
          LLR(:ng,:14*ng)=matmul(Lambda(:ng,:ng),Lx(:ng,:14*ng,nxmin,j,k))
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
          LLR(:ng,:14*ng)=Lx(:ng,:14*ng,nxmin,j,k)
          A(:ng,:ng)=real(-LLR(:ng,ng+1:2*ng),4)
        else if(iqf1 == -4) then
          call XABORT('NSSANM3: SYME boundary condition is not supported.(1)')
        else
          call XABORT('NSSANM3: illegal left X-boundary condition.')
        endif
        if(iqf1 /= -2) then
          A(:ng,ng+1)=real(matmul(LLR(:ng,:ng),savg(ind1,:ng)),4)
          do ig=1,ng
            do jg=1,ng
              if(jym /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,3*ng+jg)*savg(7*ll4f+ll4x+jym,jg),4)
              if(jym_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,4*ng+jg)*savg(7*ll4f+ll4x+jym_p,jg),4)
              if(jyp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,6*ng+jg)*savg(7*ll4f+ll4x+jyp,jg),4)
              if(jyp_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,7*ng+jg)*savg(7*ll4f+ll4x+jyp_p,jg),4)
              !
              if(jzm /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,9*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzm,jg),4)
              if(jzm_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,10*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzm_p,jg),4)
              if(jzp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,12*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzp,jg),4)
              if(jzp_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,13*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzp_p,jg),4)
            enddo
          enddo
        endif
        call ALSB(ng,1,A,ier,ng)
        if(ier /= 0) call XABORT('NSSANM3: singular matrix.(1)')
        if(jxm /= 0) savg(7*ll4f+jxm,:ng)=A(:ng,ng+1)
        !
        ! two-node relations
        do i=nxmin,nxmax-1
          ind1=idl(i,j,k)
          if(ind1 == 0) call XABORT('NSSANM3: invalid idl index.(2)')
          ind2=idl(i+1,j,k)
          if(ind2 == 0) call XABORT('NSSANM3: invalid idl index.(3)')
          if(kn(1,i+1,j,k) /= kn(2,i,j,k)) call XABORT('NSSANM3: invalid kn index.(1)')
          if(iqfr(2,i,j,k) /= 0) call XABORT('NSSANM3: invalid iqfr index.(1)')
          if(iqfr(1,i+1,j,k) /= 0) call XABORT('NSSANM3: invalid iqfr index.(2)')
          jxm=kn(1,i,j,k) ; jxp=kn(2,i,j,k) ; jym=kn(3,i,j,k) ; jyp=kn(4,i,j,k)
          jzm=kn(5,i,j,k) ; jzp=kn(6,i,j,k)
          jym_m=0 ; jyp_m=0 ; jym_pp=0 ; jyp_pp=0
          jzm_m=0 ; jzp_m=0 ; jzm_pp=0 ; jzp_pp=0
          if((i == 1).and.(iqfr(1,1,j,k) == -2)) then
            jym_m=kn(3,1,j,k) ; jyp_m=kn(4,1,j,k)
            jzm_m=kn(5,1,j,k) ; jzp_m=kn(6,1,j,k)
          else if(i > 1) then
            jym_m=kn(3,i-1,j,k) ; jyp_m=kn(4,i-1,j,k)
            jzm_m=kn(5,i-1,j,k) ; jzp_m=kn(6,i-1,j,k)
          endif
          jym_p=kn(3,i+1,j,k) ; jyp_p=kn(4,i+1,j,k)
          jzm_p=kn(5,i+1,j,k) ; jzp_p=kn(6,i+1,j,k)
          if((i == nx-1).and.(iqfr(2,nx,j,k) == -2)) then
            jym_pp=kn(3,nx,j,k) ; jyp_pp=kn(4,nx,j,k)
            jzm_pp=kn(5,nx,j,k) ; jzp_pp=kn(6,nx,j,k)
          else if(i < nx-1) then
            jym_pp=kn(3,i+2,j,k) ; jyp_pp=kn(4,i+2,j,k)
            jzm_pp=kn(5,i+2,j,k) ; jzp_pp=kn(6,i+2,j,k)
          endif
          !
          A(:ng,:ng+1)=0.0
          ! node i
          LLR(:ng,:14*ng)=matmul(fd(mat(i,j,k),2,:ng,:ng),Rx(:ng,:14*ng,i,j,k))
          do ig=1,ng
            A(:ng,ig)=A(:ng,ig)+real(LLR(:ng,ng+ig),4)
            do jg=1,ng
              A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,jg)*savg(ind1,jg),4)
              if(jym_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,2*ng+jg)*savg(7*ll4f+ll4x+jym_m,jg),4)
              if(jym /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,3*ng+jg)*savg(7*ll4f+ll4x+jym,jg),4)
              if(jym_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,4*ng+jg)*savg(7*ll4f+ll4x+jym_p,jg),4)
              if(jyp_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,5*ng+jg)*savg(7*ll4f+ll4x+jyp_m,jg),4)
              if(jyp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,6*ng+jg)*savg(7*ll4f+ll4x+jyp,jg),4)
              if(jyp_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,7*ng+jg)*savg(7*ll4f+ll4x+jyp_p,jg),4)
              !
              if(jzm_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,8*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzm_m,jg),4)
              if(jzm /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,9*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzm,jg),4)
              if(jzm_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,10*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzm_p,jg),4)
              if(jzp_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,11*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzp_m,jg),4)
              if(jzp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,12*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzp,jg),4)
              if(jzp_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,13*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzp_p,jg),4)
          enddo
        enddo
        ! node i+1
        LLR(:ng,:14*ng)=matmul(fd(mat(i+1,j,k),1,:ng,:ng),Lx(:ng,:14*ng,i+1,j,k))
        do ig=1,ng
          A(:ng,ig)=A(:ng,ig)+real(-LLR(:ng,ng+ig),4)
          do jg=1,ng
              A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,jg)*savg(ind2,jg),4)
              if(jym /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,2*ng+jg)*savg(7*ll4f+ll4x+jym,jg),4)
              if(jym_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,3*ng+jg)*savg(7*ll4f+ll4x+jym_p,jg),4)
              if(jym_pp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,4*ng+jg)*savg(7*ll4f+ll4x+jym_pp,jg),4)
              if(jyp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,5*ng+jg)*savg(7*ll4f+ll4x+jyp,jg),4)
              if(jyp_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,6*ng+jg)*savg(7*ll4f+ll4x+jyp_p,jg),4)
              if(jyp_pp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,7*ng+jg)*savg(7*ll4f+ll4x+jyp_pp,jg),4)
              !
              if(jzm /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,8*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzm,jg),4)
              if(jzm_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,9*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzm_p,jg),4)
              if(jzm_pp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,10*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzm_pp,jg),4)
              if(jzp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,11*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzp,jg),4)
              if(jzp_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,12*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzp_p,jg),4)
              if(jzp_pp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,13*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzp_pp,jg),4)
            enddo
          enddo
          call ALSB(ng,1,A,ier,ng)
          if(ier /= 0) call XABORT('NSSANM3: singular matrix.(2)')
          if(jxp /= 0) savg(7*ll4f+jxp,:ng)=A(:ng,ng+1)
        enddo
        !
        ! one-node relation at right
        ind1=idl(nxmax,j,k)
        if(ind1 == 0) call XABORT('NSSANM3: invalid idl index.(4)')
        iqf2=iqfr(2,nxmax,j,k)
        jxm=kn(1,nxmax,j,k) ; jxp=kn(2,nxmax,j,k) ; jym=kn(3,nxmax,j,k) ; jyp=kn(4,nxmax,j,k)
        jzm=kn(5,nxmax,j,k) ; jzp=kn(6,nxmax,j,k)
        jym_m=0 ; jyp_m=0 ; jzm_m=0 ; jzp_m=0
        if(nxmax > 1) then
          jym_m=kn(3,nxmax-1,j,k) ; jyp_m=kn(4,nxmax-1,j,k)
          jzm_m=kn(5,nxmax-1,j,k) ; jzp_m=kn(6,nxmax-1,j,k)
        endif
        A(:ng,:ng+1)=0.0
        LLR(:ng,:14*ng)=0.0
        if((iqf2 > 0).or.(iqf2 == -1)) then
          ! physical albedo
          Lambda(:ng,:ng)=0.0
          do ig=1,ng
            Lambda(ig,ig)=qfr(2,nxmax,j,k,ig)
          enddo
          LLR(:ng,:14*ng)=matmul(Lambda(:ng,:ng),Rx(:ng,:14*ng,nxmax,j,k))
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
          LLR(:ng,:14*ng)=Rx(:ng,:14*ng,nxmax,j,k)
          A(:ng,:ng)=real(LLR(:ng,ng+1:2*ng),4)
        else if(iqf2 == -4) then
          call XABORT('NSSANM3: SYME boundary condition is not supported.(2)')
        else
          call XABORT('NSSANM3: illegal right X-boundary condition.')
        endif
        if(iqf2 /= -2) then
          A(:ng,ng+1)=real(matmul(-LLR(:ng,:ng),savg(ind1,:ng)),4)
          do ig=1,ng
            do jg=1,ng
              if(jym_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,2*ng+jg)*savg(7*ll4f+ll4x+jym_m,jg),4)
              if(jym /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,3*ng+jg)*savg(7*ll4f+ll4x+jym,jg),4)
              if(jyp_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,5*ng+jg)*savg(7*ll4f+ll4x+jyp_m,jg),4)
              if(jyp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,6*ng+jg)*savg(7*ll4f+ll4x+jyp,jg),4)
              !
              if(jzm_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,8*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzm_m,jg),4)
              if(jzm /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,9*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzm,jg),4)
              if(jzp_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,11*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzp_m,jg),4)
              if(jzp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,12*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzp,jg),4)
            enddo
          enddo
        endif
        call ALSB(ng,1,A,ier,ng)
        if(ier /= 0) call XABORT('NSSANM3: singular matrix.(3)')
        if(jxp /= 0) savg(7*ll4f+jxp,:ng)=A(:ng,ng+1)
      enddo
    enddo
    !----
    !  one- and two-node relations along Y axis
    !----
    do i=1,nx
      do k=1,nz
        nymin=1
        do j=1,ny
          if(mat(i,j,k) > 0) exit
          nymin=j+1
        enddo
        if(nymin > ny) cycle
        nymax=ny
        do j=ny,1,-1
          if(mat(i,j,k) > 0) exit
          nymax=j-1
        enddo
        ! one-node relation at left
        ind1=idl(i,nymin,k)
        if(ind1 == 0) call XABORT('NSSANM3: invalid idl index.(5)')
        iqf3=iqfr(3,i,nymin,k)
        jxm=kn(1,i,nymin,k) ; jxp=kn(2,i,nymin,k) ; jym=kn(3,i,nymin,k) ; jyp=kn(4,i,nymin,k)
        jzm=kn(5,i,nymin,k) ; jzp=kn(6,i,nymin,k)
        jxm_p=0 ; jxp_p=0 ; jzm_p=0 ; jzp_p=0
        if(nymin < ny) then
          jxm_p=kn(1,i,nymin+1,k) ; jxp_p=kn(2,i,nymin+1,k)
          jzm_p=kn(5,i,nymin+1,k) ; jzp_p=kn(6,i,nymin+1,k)
        endif
        A(:ng,:ng+1)=0.0
        LLR(:ng,:14*ng)=0.0
        if((iqf3 > 0).or.(iqf3 == -1)) then
          ! physical albedo
          Lambda(:ng,:ng)=0.0
          do ig=1,ng
            Lambda(ig,ig)=qfr(3,i,nymin,k,ig)
          enddo
          LLR(:ng,:14*ng)=matmul(Lambda(:ng,:ng),Ly(:ng,:14*ng,i,nymin,k))
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
          LLR(:ng,:14*ng)=Ly(:ng,:14*ng,i,nymin,k)
          A(:ng,:ng)=real(-LLR(:ng,ng+1:2*ng),4)
        else if(iqf3 == -4) then
          call XABORT('NSSANM3: SYME boundary condition is not supported.(3)')
        else
          call XABORT('NSSANM3: illegal left Y-boundary condition.')
        endif
        if(iqf3 /= -2) then
          A(:ng,ng+1)=real(matmul(LLR(:ng,:ng),savg(ind1,:ng)),4)
          do ig=1,ng
            do jg=1,ng
              if(jzm /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,3*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzm,jg),4)
              if(jzm_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,4*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzm_p,jg),4)
              if(jzp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,6*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzp,jg),4)
              if(jzp_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,7*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzp_p,jg),4)
              !
              if(jxm /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,9*ng+jg)*savg(7*ll4f+jxm,jg),4)
              if(jxm_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,10*ng+jg)*savg(7*ll4f+jxm_p,jg),4)
              if(jxp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,12*ng+jg)*savg(7*ll4f+jxp,jg),4)
              if(jxp_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,13*ng+jg)*savg(7*ll4f+jxp_p,jg),4)
            enddo
          enddo
        endif
        call ALSB(ng,1,A,ier,ng)
        if(ier /= 0) call XABORT('NSSANM3: singular matrix.(4)')
        if(jym /= 0) savg(7*ll4f+ll4x+jym,:ng)=A(:ng,ng+1)
        !
        ! two-node relations
        do j=nymin,nymax-1
          ind1=idl(i,j,k)
          if(ind1 == 0) call XABORT('NSSANM3: invalid idl index.(6)')
          ind2=idl(i,j+1,k)
          if(ind2 == 0) call XABORT('NSSANM3: invalid idl index.(7)')
          if(kn(3,i,j+1,k) /= kn(4,i,j,k)) call XABORT('NSSANM3: invalid kn index.(2)')
          if(iqfr(4,i,j,k) /= 0) call XABORT('NSSANM3: invalid iqfr index.(3)')
          if(iqfr(3,i,j+1,k) /= 0) call XABORT('NSSANM3: invalid iqfr index.(4)')
          jxm=kn(1,i,j,k) ; jxp=kn(2,i,j,k) ; jym=kn(3,i,j,k) ; jyp=kn(4,i,j,k)
          jzm=kn(5,i,j,k) ; jzp=kn(6,i,j,k)
          jxm_m=0 ; jxp_m=0 ; jxm_pp=0 ; jxp_pp=0
          jzm_m=0 ; jzp_m=0 ; jzm_pp=0 ; jzp_pp=0
          if((j == 1).and.(iqfr(3,i,1,k) == -2)) then
            jxm_m=kn(1,i,1,k) ; jxp_m=kn(2,i,1,k)
            jzm_m=kn(5,i,1,k) ; jzp_m=kn(6,i,1,k)
          else if(j > 1) then
            jxm_m=kn(1,i,j-1,k) ; jxp_m=kn(2,i,j-1,k)
            jzm_m=kn(5,i,j-1,k) ; jzp_m=kn(6,i,j-1,k)
          endif
          jxm_p=kn(1,i,j+1,k) ; jxp_p=kn(2,i,j+1,k)
          jzm_p=kn(5,i,j+1,k) ; jzp_p=kn(6,i,j+1,k)
          if((j == ny-1).and.(iqfr(4,i,ny,k) == -2)) then
            jxm_pp=kn(1,i,ny,k) ; jxp_pp=kn(2,i,ny,k)
            jzm_pp=kn(5,i,ny,k) ; jzp_pp=kn(6,i,ny,k)
          else if(j < ny-1) then
            jxm_pp=kn(1,i,j+2,k) ; jxp_pp=kn(2,i,j+2,k)
            jzm_pp=kn(5,i,j+2,k) ; jzp_pp=kn(6,i,j+2,k)
          endif
          !
          A(:ng,:ng+1)=0.0
          ! node j
          LLR(:ng,:14*ng)=matmul(fd(mat(i,j,k),4,:ng,:ng),Ry(:ng,:14*ng,i,j,k))
          do ig=1,ng
            A(:ng,ig)=A(:ng,ig)+real(LLR(:ng,ng+ig),4)
            do jg=1,ng
              A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,jg)*savg(ind1,jg),4)
              if(jzm_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,2*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzm_m,jg),4)
              if(jzm /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,3*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzm,jg),4)
              if(jzm_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,4*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzm_p,jg),4)
              if(jzp_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,5*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzp_m,jg),4)
              if(jzp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,6*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzp,jg),4)
              if(jzp_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,7*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzp_p,jg),4)
              !
              if(jxm_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,8*ng+jg)*savg(7*ll4f+jxm_m,jg),4)
              if(jxm /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,9*ng+jg)*savg(7*ll4f+jxm,jg),4)
              if(jxm_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,10*ng+jg)*savg(7*ll4f+jxm_p,jg),4)
              if(jxp_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,11*ng+jg)*savg(7*ll4f+jxp_m,jg),4)
              if(jxp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,12*ng+jg)*savg(7*ll4f+jxp,jg),4)
              if(jxp_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,13*ng+jg)*savg(7*ll4f+jxp_p,jg),4)
            enddo
          enddo
          ! node j+1
          LLR(:ng,:14*ng)=matmul(fd(mat(i,j+1,k),3,:ng,:ng),Ly(:ng,:14*ng,i,j+1,k))
          do ig=1,ng
            A(:ng,ig)=A(:ng,ig)+real(-LLR(:ng,ng+ig),4)
            do jg=1,ng
              A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,jg)*savg(ind2,jg),4)
              if(jzm /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,2*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzm,jg),4)
              if(jzm_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,3*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzm_p,jg),4)
              if(jzm_pp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,4*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzm_pp,jg),4)
              if(jzp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,5*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzp,jg),4)
              if(jzp_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,6*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzp_p,jg),4)
              if(jzp_pp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,7*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzp_pp,jg),4)
              !
              if(jxm /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,8*ng+jg)*savg(7*ll4f+jxm,jg),4)
              if(jxm_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,9*ng+jg)*savg(7*ll4f+jxm_p,jg),4)
              if(jxm_pp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,10*ng+jg)*savg(7*ll4f+jxm_pp,jg),4)
              if(jxp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,11*ng+jg)*savg(7*ll4f+jxp,jg),4)
              if(jxp_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,12*ng+jg)*savg(7*ll4f+jxp_p,jg),4)
              if(jxp_pp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,13*ng+jg)*savg(7*ll4f+jxp_pp,jg),4)
            enddo
          enddo
          call ALSB(ng,1,A,ier,ng)
          if(ier /= 0) call XABORT('NSSANM3: singular matrix.(5)')
          if(jyp /= 0) savg(7*ll4f+ll4x+jyp,:ng)=A(:ng,ng+1)
        enddo
        !
        ! one-node relation at right
        ind1=idl(i,nymax,k)
        if(ind1 == 0) call XABORT('NSSANM3: invalid idl index.(8)')
        iqf4=iqfr(4,i,nymax,k)
        jxm=kn(1,i,nymax,k) ; jxp=kn(2,i,nymax,k) ; jym=kn(3,i,nymax,k) ; jyp=kn(4,i,nymax,k)
        jzm=kn(5,i,nymax,k) ; jzp=kn(6,i,nymax,k)
        jxm_m=0 ; jxp_m=0 ; jzm_m=0 ; jzp_m=0
        if(nymax > 1) then
          jxm_m=kn(1,i,nymax-1,k) ; jxp_m=kn(2,i,nymax-1,k)
          jzm_m=kn(5,i,nymax-1,k) ; jzp_m=kn(6,i,nymax-1,k)
        endif
        A(:ng,:ng+1)=0.0
        LLR(:ng,:14*ng)=0.0
        if((iqf4 > 0).or.(iqf4 == -1)) then
          ! physical albedo
          Lambda(:ng,:ng)=0.0
          do ig=1,ng
            Lambda(ig,ig)=qfr(4,i,nymax,k,ig)
          enddo
          LLR(:ng,:14*ng)=matmul(Lambda(:ng,:ng),Ry(:ng,:14*ng,i,nymax,k))
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
          LLR(:ng,:14*ng)=Ry(:ng,:14*ng,i,nymax,k)
          A(:ng,:ng)=real(LLR(:ng,ng+1:2*ng),4)
        else if(iqf4 == -4) then
          call XABORT('NSSANM3: SYME boundary condition is not supported.(4)')
        else
          call XABORT('NSSANM3: illegal right Y-boundary condition.')
        endif
        if(iqf4 /= -2) then
          A(:ng,ng+1)=real(matmul(-LLR(:ng,:ng),savg(ind1,:ng)),4)
          do ig=1,ng
            do jg=1,ng
              if(jzm_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,2*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzm_m,jg),4)
              if(jzm /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,3*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzm,jg),4)
              if(jzp_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,5*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzp_m,jg),4)
              if(jzp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,6*ng+jg)*savg(7*ll4f+ll4x+ll4y+jzp,jg),4)
              !
              if(jxm_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,8*ng+jg)*savg(7*ll4f+jxm_m,jg),4)
              if(jxm /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,9*ng+jg)*savg(7*ll4f+jxm,jg),4)
              if(jxp_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,11*ng+jg)*savg(7*ll4f+jxp_m,jg),4)
              if(jxp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,12*ng+jg)*savg(7*ll4f+jxp,jg),4)
            enddo
          enddo
        endif
        call ALSB(ng,1,A,ier,ng)
        if(ier /= 0) call XABORT('NSSANM3: singular matrix.(6)')
        if(jyp /= 0) savg(7*ll4f+ll4x+jyp,:ng)=A(:ng,ng+1)
      enddo
    enddo
    !----
    !  one- and two-node relations along Z axis
    !----
    do j=1,ny
      do i=1,nx
        nzmin=1
        do k=1,nz
          if(mat(i,j,k) > 0) exit
          nzmin=k+1
        enddo
        if(nzmin > nz) cycle
        nzmax=nz
        do k=nz,1,-1
          if(mat(i,j,k) > 0) exit
          nzmax=k-1
        enddo
        ! one-node relation at left
        ind1=idl(i,j,nzmin)
        if(ind1 == 0) call XABORT('NSSANM3: invalid idl index.(9)')
        iqf5=iqfr(5,i,j,nzmin)
        jxm=kn(1,i,j,nzmin) ; jxp=kn(2,i,j,nzmin) ; jym=kn(3,i,j,nzmin) ; jyp=kn(4,i,j,nzmin)
        jzm=kn(5,i,j,nzmin) ; jzp=kn(6,i,j,nzmin)
        jxm_p=0 ; jxp_p=0 ; jym_p=0 ; jyp_p=0
        if(nzmin < nz) then
          jxm_p=kn(1,i,j,nzmin+1) ; jxp_p=kn(2,i,j,nzmin+1)
          jym_p=kn(3,i,j,nzmin+1) ; jyp_p=kn(4,i,j,nzmin+1)
        endif
        A(:ng,:ng+1)=0.0
        LLR(:ng,:14*ng)=0.0
        if((iqf5 > 0).or.(iqf5 == -1)) then
          ! physical albedo
          Lambda(:ng,:ng)=0.0
          do ig=1,ng
            Lambda(ig,ig)=qfr(5,i,j,nzmin,ig)
          enddo
          LLR(:ng,:14*ng)=matmul(Lambda(:ng,:ng),Lz(:ng,:14*ng,i,j,nzmin))
          A(:ng,:ng)=-real(LLR(:ng,ng+1:2*ng),4)
          do ig=1,ng
            A(ig,ig)=-1.0+A(ig,ig)
          enddo
        else if(iqf5 == -2) then
          ! zero net current
          do ig=1,ng
            A(ig,ig)=1.0
          enddo
        else if(iqf5 == -3) then
          ! zero flux
          LLR(:ng,:14*ng)=Lz(:ng,:14*ng,i,j,nzmin)
          A(:ng,:ng)=real(-LLR(:ng,ng+1:2*ng),4)
        else if(iqf5 == -4) then
          call XABORT('NSSANM3: SYME boundary condition is not supported.(5)')
        else
          call XABORT('NSSANM3: illegal left Z-boundary condition.')
        endif
        if(iqf5 /= -2) then
          A(:ng,ng+1)=real(matmul(LLR(:ng,:ng),savg(ind1,:ng)),4)
          do ig=1,ng
            do jg=1,ng
              if(jxm /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,3*ng+jg)*savg(7*ll4f+jxm,jg),4)
              if(jxm_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,4*ng+jg)*savg(7*ll4f+jxm_p,jg),4)
              if(jxp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,6*ng+jg)*savg(7*ll4f+jxp,jg),4)
              if(jxp_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,7*ng+jg)*savg(7*ll4f+jxp_p,jg),4)
              !
              if(jym /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,9*ng+jg)*savg(7*ll4f+ll4x+jym,jg),4)
              if(jym_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,10*ng+jg)*savg(7*ll4f+ll4x+jym_p,jg),4)
              if(jyp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,12*ng+jg)*savg(7*ll4f+ll4x+jyp,jg),4)
              if(jyp_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,13*ng+jg)*savg(7*ll4f+ll4x+jyp_p,jg),4)
            enddo
          enddo
        endif
        call ALSB(ng,1,A,ier,ng)
        if(ier /= 0) call XABORT('NSSANM3: singular matrix.(7)')
        if(jzm /= 0) savg(7*ll4f+ll4x+ll4y+jzm,:ng)=A(:ng,ng+1)
        !
        ! two-node relations
        do k=nzmin,nzmax-1
          ind1=idl(i,j,k)
          if(ind1 == 0) call XABORT('NSSANM3: invalid idl index.(10)')
          ind2=idl(i,j,k+1)
          if(ind2 == 0) call XABORT('NSSANM3: invalid idl index.(11)')
          if(kn(5,i,j,k+1) /= kn(6,i,j,k)) call XABORT('NSSANM3: invalid kn index.(3)')
          if(iqfr(6,i,j,k) /= 0) call XABORT('NSSANM3: invalid iqfr index.(5)')
          if(iqfr(5,i,j,k+1) /= 0) call XABORT('NSSANM3: invalid iqfr index.(6)')
          jxm=kn(1,i,j,k) ; jxp=kn(2,i,j,k) ; jym=kn(3,i,j,k) ; jyp=kn(4,i,j,k)
          jzm=kn(5,i,j,k) ; jzp=kn(6,i,j,k)
          jxm_m=0 ; jxp_m=0 ; jxm_pp=0 ; jxp_pp=0
          jym_m=0 ; jyp_m=0 ; jym_pp=0 ; jyp_pp=0
          if((k == 1).and.(iqfr(5,i,j,1) == -2)) then
            jxm_m=kn(1,i,j,1) ; jxp_m=kn(2,i,j,1)
            jym_m=kn(3,i,j,1) ; jyp_m=kn(4,i,j,1)
          else if(k > 1) then
            jxm_m=kn(1,i,j,k-1) ; jxp_m=kn(2,i,j,k-1)
            jym_m=kn(3,i,j,k-1) ; jyp_m=kn(4,i,j,k-1)
          endif
          jxm_p=kn(1,i,j,k+1) ; jxp_p=kn(2,i,j,k+1)
          jym_p=kn(3,i,j,k+1) ; jyp_p=kn(4,i,j,k+1)
          if((k == nz-1).and.(iqfr(6,i,j,nz) == -2)) then
            jxm_pp=kn(1,i,j,nz) ; jxp_pp=kn(2,i,j,nz)
            jym_pp=kn(3,i,j,nz) ; jyp_pp=kn(4,i,j,nz)
          else if(k < nz-1) then
            jxm_pp=kn(1,i,j,k+2) ; jxp_pp=kn(2,i,j,k+2)
            jym_pp=kn(3,i,j,k+2) ; jyp_pp=kn(4,i,j,k+2)
          endif
          !
          A(:ng,:ng+1)=0.0
          ! node i
          LLR(:ng,:14*ng)=matmul(fd(mat(i,j,k),6,:ng,:ng),Rz(:ng,:14*ng,i,j,k))
          do ig=1,ng
            A(:ng,ig)=A(:ng,ig)+real(LLR(:ng,ng+ig),4)
            do jg=1,ng
              A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,jg)*savg(ind1,jg),4)
              if(jxm_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,2*ng+jg)*savg(7*ll4f+jxm_m,jg),4)
              if(jxm /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,3*ng+jg)*savg(7*ll4f+jxm,jg),4)
              if(jxm_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,4*ng+jg)*savg(7*ll4f+jxm_p,jg),4)
              if(jxp_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,5*ng+jg)*savg(7*ll4f+jxp_m,jg),4)
              if(jxp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,6*ng+jg)*savg(7*ll4f+jxp,jg),4)
              if(jxp_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,7*ng+jg)*savg(7*ll4f+jxp_p,jg),4)
              !
              if(jym_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,8*ng+jg)*savg(7*ll4f+ll4x+jym_m,jg),4)
              if(jym /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,9*ng+jg)*savg(7*ll4f+ll4x+jym,jg),4)
              if(jym_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,10*ng+jg)*savg(7*ll4f+ll4x+jym_p,jg),4)
              if(jyp_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,11*ng+jg)*savg(7*ll4f+ll4x+jyp_m,jg),4)
              if(jyp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,12*ng+jg)*savg(7*ll4f+ll4x+jyp,jg),4)
              if(jyp_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,13*ng+jg)*savg(7*ll4f+ll4x+jyp_p,jg),4)
            enddo
          enddo
          ! node i+1
          LLR(:ng,:14*ng)=matmul(fd(mat(i,j,k+1),5,:ng,:ng),Lz(:ng,:14*ng,i,j,k+1))
          do ig=1,ng
            A(:ng,ig)=A(:ng,ig)+real(-LLR(:ng,ng+ig),4)
            do jg=1,ng
              A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,jg)*savg(ind2,jg),4)
              if(jxm /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,2*ng+jg)*savg(7*ll4f+jxm,jg),4)
              if(jxm_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,3*ng+jg)*savg(7*ll4f+jxm_p,jg),4)
              if(jxm_pp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,4*ng+jg)*savg(7*ll4f+jxm_pp,jg),4)
              if(jxp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,5*ng+jg)*savg(7*ll4f+jxp,jg),4)
              if(jxp_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,6*ng+jg)*savg(7*ll4f+jxp_p,jg),4)
              if(jxp_pp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,7*ng+jg)*savg(7*ll4f+jxp_pp,jg),4)
              !
              if(jym /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,8*ng+jg)*savg(7*ll4f+ll4x+jym,jg),4)
              if(jym_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,9*ng+jg)*savg(7*ll4f+ll4x+jym_p,jg),4)
              if(jym_pp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,10*ng+jg)*savg(7*ll4f+ll4x+jym_pp,jg),4)
              if(jyp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,11*ng+jg)*savg(7*ll4f+ll4x+jyp,jg),4)
              if(jyp_p /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,12*ng+jg)*savg(7*ll4f+ll4x+jyp_p,jg),4)
              if(jyp_pp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(LLR(ig,13*ng+jg)*savg(7*ll4f+ll4x+jyp_pp,jg),4)
            enddo
          enddo
          call ALSB(ng,1,A,ier,ng)
          if(ier /= 0) call XABORT('NSSANM3: singular matrix.(8)')
          if(jzp /= 0) savg(7*ll4f+ll4x+ll4y+jzp,:ng)=A(:ng,ng+1)
        enddo
        !
        ! one-node relation at right
        ind1=idl(i,j,nzmax)
        if(ind1 == 0) call XABORT('NSSANM3: invalid idl index.(12)')
        iqf6=iqfr(6,i,j,nzmax)
        jxm=kn(1,i,j,nzmax) ; jxp=kn(2,i,j,nzmax) ; jym=kn(3,i,j,nzmax) ; jyp=kn(4,i,j,nzmax)
        jzm=kn(5,i,j,nzmax) ; jzp=kn(6,i,j,nzmax)
        jxm_m=0 ; jxp_m=0 ; jym_m=0 ; jyp_m=0
        if(nzmax > 1) then
          jxm_m=kn(1,i,j,nzmax-1) ; jxp_m=kn(2,i,j,nzmax-1)
          jym_m=kn(3,i,j,nzmax-1) ; jyp_m=kn(4,i,j,nzmax-1)
        endif
        A(:ng,:ng+1)=0.0
        LLR(:ng,:14*ng)=0.0
        if((iqf6 > 0).or.(iqf6 == -1)) then
          ! physical albedo
          Lambda(:ng,:ng)=0.0
          do ig=1,ng
            Lambda(ig,ig)=qfr(6,i,j,nzmax,ig)
          enddo
          LLR(:ng,:14*ng)=matmul(Lambda(:ng,:ng),Rz(:ng,:14*ng,i,j,nzmax))
          A(:ng,:ng)=real(LLR(:ng,ng+1:2*ng),4)
          do ig=1,ng
            A(ig,ig)=-1.0+A(ig,ig)
          enddo
        else if(iqf6 == -2) then
          ! zero net current
          do ig=1,ng
            A(ig,ig)=1.0
          enddo
        else if(iqf6 == -3) then
          ! zero flux
          LLR(:ng,:14*ng)=Rz(:ng,:14*ng,i,j,nzmax)
          A(:ng,:ng)=real(LLR(:ng,ng+1:2*ng),4)
        else if(iqf6 == -4) then
          call XABORT('NSSANM3: SYME boundary condition is not supported.(6)')
        else
          call XABORT('NSSANM3: illegal right Z-boundary condition.')
        endif
        if(iqf6 /= -2) then
          A(:ng,ng+1)=real(matmul(-LLR(:ng,:ng),savg(ind1,:ng)),4)
          do ig=1,ng
            do jg=1,ng
              if(jxm_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,2*ng+jg)*savg(7*ll4f+jxm_m,jg),4)
              if(jxm /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,3*ng+jg)*savg(7*ll4f+jxm,jg),4)
              if(jxp_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,5*ng+jg)*savg(7*ll4f+jxp_m,jg),4)
              if(jxp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,6*ng+jg)*savg(7*ll4f+jxp,jg),4)
              !
              if(jym_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,8*ng+jg)*savg(7*ll4f+ll4x+jym_m,jg),4)
              if(jym /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,9*ng+jg)*savg(7*ll4f+ll4x+jym,jg),4)
              if(jyp_m /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,11*ng+jg)*savg(7*ll4f+ll4x+jyp_m,jg),4)
              if(jyp /= 0) A(ig,ng+1)=A(ig,ng+1)+real(-LLR(ig,12*ng+jg)*savg(7*ll4f+ll4x+jyp,jg),4)
            enddo
          enddo
        endif
        call ALSB(ng,1,A,ier,ng)
        if(ier /= 0) call XABORT('NSSANM3: singular matrix.(9)')
        if(jzp /= 0) savg(7*ll4f+ll4x+ll4y+jzp,:ng)=A(:ng,ng+1)
      enddo
    enddo
    !----
    !  end of transverse current iterations
    !----
  enddo
  !----
  ! compute boundary fluxes
  !----
  do k=1,nz
    do j=1,ny
      do i=1,nx
        ind1=idl(i,j,k)
        if(ind1 == 0) cycle
        jxm=kn(1,i,j,k) ; jxp=kn(2,i,j,k) ; jym=kn(3,i,j,k) ; jyp=kn(4,i,j,k)
        jzm=kn(5,i,j,k) ; jzp=kn(6,i,j,k)
        !
        jym_m=0 ; jyp_m=0 ; jym_p=0 ; jyp_p=0
        jzm_m=0 ; jzp_m=0 ; jzm_p=0 ; jzp_p=0
        if((i == 1).and.(iqfr(1,1,j,k) == -2)) then
          jym_m=kn(3,1,j,k) ; jyp_m=kn(4,1,j,k)
          jzm_m=kn(5,1,j,k) ; jzp_m=kn(6,1,j,k)
        else if(i > 1) then
          jym_m=kn(3,i-1,j,k) ; jyp_m=kn(4,i-1,j,k)
          jzm_m=kn(5,i-1,j,k) ; jzp_m=kn(6,i-1,j,k)
        endif
        if((i == nx).and.(iqfr(2,nx,j,k) == -2)) then
          jym_p=kn(3,nx,j,k) ; jyp_p=kn(4,nx,j,k)
          jzm_p=kn(5,nx,j,k) ; jzp_p=kn(6,nx,j,k)
        else if(i < nx) then
          jym_p=kn(3,i+1,j,k) ; jyp_p=kn(4,i+1,j,k)
          jzm_p=kn(5,i+1,j,k) ; jzp_p=kn(6,i+1,j,k)
        endif
        ! x- relations
        savg(ll4f+ind1,:ng)=real(matmul(Lx(:ng,:ng,i,j,k),savg(ind1,:ng)),4)
        if(jxm /= 0) savg(ll4f+ind1,:ng)=savg(ll4f+ind1,:ng)+ &
          & real(matmul(Lx(:ng,ng+1:2*ng,i,j,k),savg(7*ll4f+jxm,:ng)),4)
        !
        if(jym_m /= 0) savg(ll4f+ind1,:ng)=savg(ll4f+ind1,:ng)+ &
          & real(matmul(Lx(:ng,2*ng+1:3*ng,i,j,k),savg(7*ll4f+ll4x+jym_m,:ng)),4)
        if(jym /= 0) savg(ll4f+ind1,:ng)=savg(ll4f+ind1,:ng)+ &
          & real(matmul(Lx(:ng,3*ng+1:4*ng,i,j,k),savg(7*ll4f+ll4x+jym,:ng)),4)
        if(jym_p /= 0) savg(ll4f+ind1,:ng)=savg(ll4f+ind1,:ng)+ &
          & real(matmul(Lx(:ng,4*ng+1:5*ng,i,j,k),savg(7*ll4f+ll4x+jym_p,:ng)),4)
        if(jyp_m /= 0) savg(ll4f+ind1,:ng)=savg(ll4f+ind1,:ng)+ &
          & real(matmul(Lx(:ng,5*ng+1:6*ng,i,j,k),savg(7*ll4f+ll4x+jyp_m,:ng)),4)
        if(jyp /= 0) savg(ll4f+ind1,:ng)=savg(ll4f+ind1,:ng)+ &
          & real(matmul(Lx(:ng,6*ng+1:7*ng,i,j,k),savg(7*ll4f+ll4x+jyp,:ng)),4)
        if(jyp_p /= 0) savg(ll4f+ind1,:ng)=savg(ll4f+ind1,:ng)+ &
          & real(matmul(Lx(:ng,7*ng+1:8*ng,i,j,k),savg(7*ll4f+ll4x+jyp_p,:ng)),4)
        !
        if(jzm_m /= 0) savg(ll4f+ind1,:ng)=savg(ll4f+ind1,:ng)+ &
          & real(matmul(Lx(:ng,8*ng+1:9*ng,i,j,k),savg(7*ll4f+ll4x+ll4y+jzm_m,:ng)),4)
        if(jzm /= 0) savg(ll4f+ind1,:ng)=savg(ll4f+ind1,:ng)+ &
          & real(matmul(Lx(:ng,9*ng+1:10*ng,i,j,k),savg(7*ll4f+ll4x+ll4y+jzm,:ng)),4)
        if(jzm_p /= 0) savg(ll4f+ind1,:ng)=savg(ll4f+ind1,:ng)+ &
          & real(matmul(Lx(:ng,10*ng+1:11*ng,i,j,k),savg(7*ll4f+ll4x+ll4y+jzm_p,:ng)),4)
        if(jzp_m /= 0) savg(ll4f+ind1,:ng)=savg(ll4f+ind1,:ng)+ &
          & real(matmul(Lx(:ng,11*ng+1:12*ng,i,j,k),savg(7*ll4f+ll4x+ll4y+jzp_m,:ng)),4)
        if(jzp /= 0) savg(ll4f+ind1,:ng)=savg(ll4f+ind1,:ng)+ &
          & real(matmul(Lx(:ng,12*ng+1:13*ng,i,j,k),savg(7*ll4f+ll4x+ll4y+jzp,:ng)),4)
        if(jzp_p /= 0) savg(ll4f+ind1,:ng)=savg(ll4f+ind1,:ng)+ &
          & real(matmul(Lx(:ng,13*ng+1:14*ng,i,j,k),savg(7*ll4f+ll4x+ll4y+jzp_p,:ng)),4)
        !
        ! x+ relations
        savg(2*ll4f+ind1,:ng)=real(matmul(Rx(:ng,:ng,i,j,k),savg(ind1,:ng)),4)
        if(jxp /= 0) savg(2*ll4f+ind1,:ng)=savg(2*ll4f+ind1,:ng)+ &
          & real(matmul(Rx(:ng,ng+1:2*ng,i,j,k),savg(7*ll4f+jxp,:ng)),4)
        !
        if(jym_m /= 0) savg(2*ll4f+ind1,:ng)=savg(2*ll4f+ind1,:ng)+ &
          & real(matmul(Rx(:ng,2*ng+1:3*ng,i,j,k),savg(7*ll4f+ll4x+jym_m,:ng)),4)
        if(jym /= 0) savg(2*ll4f+ind1,:ng)=savg(2*ll4f+ind1,:ng)+ &
          & real(matmul(Rx(:ng,3*ng+1:4*ng,i,j,k),savg(7*ll4f+ll4x+jym,:ng)),4)
        if(jym_p /= 0) savg(2*ll4f+ind1,:ng)=savg(2*ll4f+ind1,:ng)+ &
          & real(matmul(Rx(:ng,4*ng+1:5*ng,i,j,k),savg(7*ll4f+ll4x+jym_p,:ng)),4)
        if(jyp_m /= 0) savg(2*ll4f+ind1,:ng)=savg(2*ll4f+ind1,:ng)+ &
          & real(matmul(Rx(:ng,5*ng+1:6*ng,i,j,k),savg(7*ll4f+ll4x+jyp_m,:ng)),4)
        if(jyp /= 0) savg(2*ll4f+ind1,:ng)=savg(2*ll4f+ind1,:ng)+ &
          & real(matmul(Rx(:ng,6*ng+1:7*ng,i,j,k),savg(7*ll4f+ll4x+jyp,:ng)),4)
        if(jyp_p /= 0) savg(2*ll4f+ind1,:ng)=savg(2*ll4f+ind1,:ng)+ &
          & real(matmul(Rx(:ng,7*ng+1:8*ng,i,j,k),savg(7*ll4f+ll4x+jyp_p,:ng)),4)
        !
        if(jzm_m /= 0) savg(2*ll4f+ind1,:ng)=savg(2*ll4f+ind1,:ng)+ &
          & real(matmul(Rx(:ng,8*ng+1:9*ng,i,j,k),savg(7*ll4f+ll4x+ll4y+jzm_m,:ng)),4)
        if(jzm /= 0) savg(2*ll4f+ind1,:ng)=savg(2*ll4f+ind1,:ng)+ &
          & real(matmul(Rx(:ng,9*ng+1:10*ng,i,j,k),savg(7*ll4f+ll4x+ll4y+jzm,:ng)),4)
        if(jzm_p /= 0) savg(2*ll4f+ind1,:ng)=savg(2*ll4f+ind1,:ng)+ &
          & real(matmul(Rx(:ng,10*ng+1:11*ng,i,j,k),savg(7*ll4f+ll4x+ll4y+jzm_p,:ng)),4)
        if(jzp_m /= 0) savg(2*ll4f+ind1,:ng)=savg(2*ll4f+ind1,:ng)+ &
          & real(matmul(Rx(:ng,11*ng+1:12*ng,i,j,k),savg(7*ll4f+ll4x+ll4y+jzp_m,:ng)),4)
        if(jzp /= 0) savg(2*ll4f+ind1,:ng)=savg(2*ll4f+ind1,:ng)+ &
          & real(matmul(Rx(:ng,12*ng+1:13*ng,i,j,k),savg(7*ll4f+ll4x+ll4y+jzp,:ng)),4)
        if(jzp_p /= 0) savg(2*ll4f+ind1,:ng)=savg(2*ll4f+ind1,:ng)+ &
          & real(matmul(Rx(:ng,13*ng+1:14*ng,i,j,k),savg(7*ll4f+ll4x+ll4y+jzp_p,:ng)),4)
        !
        jxm_m=0 ; jxp_m=0 ; jxm_p=0 ; jxp_p=0
        jzm_m=0 ; jzp_m=0 ; jzm_p=0 ; jzp_p=0
        jxm=kn(1,i,j,k) ; jxp=kn(2,i,j,k) ; jym=kn(3,i,j,k) ; jyp=kn(4,i,j,k)
        jzm=kn(5,i,j,k) ; jzp=kn(6,i,j,k)
        if((j == 1).and.(iqfr(3,i,1,k) == -2)) then
          jxm_m=kn(1,i,1,k) ; jxp_m=kn(2,i,1,k)
          jzm_m=kn(5,i,1,k) ; jzp_m=kn(6,i,1,k)
        else if(j > 1) then
          jxm_m=kn(1,i,j-1,k) ; jxp_m=kn(2,i,j-1,k)
          jzm_m=kn(5,i,j-1,k) ; jzp_m=kn(6,i,j-1,k)
        endif
        if((j == ny).and.(iqfr(4,i,ny,k) == -2)) then
          jxm_p=kn(1,i,ny,k) ; jxp_p=kn(2,i,ny,k)
          jzm_p=kn(5,i,ny,k) ; jzp_p=kn(6,i,ny,k)
        else if(j < ny) then
          jxm_p=kn(1,i,j+1,k) ; jxp_p=kn(2,i,j+1,k)
          jzm_p=kn(5,i,j+1,k) ; jzp_p=kn(6,i,j+1,k)
        endif
        ! y- relations
        savg(3*ll4f+ind1,:ng)=real(matmul(Ly(:ng,:ng,i,j,k),savg(ind1,:ng)),4)
        if(jym /= 0) savg(3*ll4f+ind1,:ng)=savg(3*ll4f+ind1,:ng)+ &
          & real(matmul(Ly(:ng,ng+1:2*ng,i,j,k),savg(7*ll4f+ll4x+jym,:ng)),4)
        !
        if(jzm_m /= 0) savg(3*ll4f+ind1,:ng)=savg(3*ll4f+ind1,:ng)+ &
          & real(matmul(Ly(:ng,2*ng+1:3*ng,i,j,k),savg(7*ll4f+ll4x+ll4y+jzm_m,:ng)),4)
        if(jzm /= 0) savg(3*ll4f+ind1,:ng)=savg(3*ll4f+ind1,:ng)+ &
          & real(matmul(Ly(:ng,3*ng+1:4*ng,i,j,k),savg(7*ll4f+ll4x+ll4y+jzm,:ng)),4)
        if(jzm_p /= 0) savg(3*ll4f+ind1,:ng)=savg(3*ll4f+ind1,:ng)+ &
          & real(matmul(Ly(:ng,4*ng+1:5*ng,i,j,k),savg(7*ll4f+ll4x+ll4y+jzm_p,:ng)),4)
        if(jzp_m /= 0) savg(3*ll4f+ind1,:ng)=savg(3*ll4f+ind1,:ng)+ &
          & real(matmul(Ly(:ng,5*ng+1:6*ng,i,j,k),savg(7*ll4f+ll4x+ll4y+jzp_m,:ng)),4)
        if(jzp /= 0) savg(3*ll4f+ind1,:ng)=savg(3*ll4f+ind1,:ng)+ &
          & real(matmul(Ly(:ng,6*ng+1:7*ng,i,j,k),savg(7*ll4f+ll4x+ll4y+jzp,:ng)),4)
        if(jzp_p /= 0) savg(3*ll4f+ind1,:ng)=savg(3*ll4f+ind1,:ng)+ &
          & real(matmul(Ly(:ng,7*ng+1:8*ng,i,j,k),savg(7*ll4f+ll4x+ll4y+jzp_p,:ng)),4)
        !
        if(jxm_m /= 0) savg(3*ll4f+ind1,:ng)=savg(3*ll4f+ind1,:ng)+ &
          & real(matmul(Ly(:ng,8*ng+1:9*ng,i,j,k),savg(7*ll4f+jxm_m,:ng)),4)
        if(jxm /= 0) savg(3*ll4f+ind1,:ng)=savg(3*ll4f+ind1,:ng)+ &
          & real(matmul(Ly(:ng,9*ng+1:10*ng,i,j,k),savg(7*ll4f+jxm,:ng)),4)
        if(jxm_p /= 0) savg(3*ll4f+ind1,:ng)=savg(3*ll4f+ind1,:ng)+ &
          & real(matmul(Ly(:ng,10*ng+1:11*ng,i,j,k),savg(7*ll4f+jxm_p,:ng)),4)
        if(jxp_m /= 0) savg(3*ll4f+ind1,:ng)=savg(3*ll4f+ind1,:ng)+ &
          & real(matmul(Ly(:ng,11*ng+1:12*ng,i,j,k),savg(7*ll4f+jxp_m,:ng)),4)
        if(jxp /= 0) savg(3*ll4f+ind1,:ng)=savg(3*ll4f+ind1,:ng)+ &
          & real(matmul(Ly(:ng,12*ng+1:13*ng,i,j,k),savg(7*ll4f+jxp,:ng)),4)
        if(jxp_p /= 0) savg(3*ll4f+ind1,:ng)=savg(3*ll4f+ind1,:ng)+ &
          & real(matmul(Ly(:ng,13*ng+1:14*ng,i,j,k),savg(7*ll4f+jxp_p,:ng)),4)
        !
        ! y+ relations
        savg(4*ll4f+ind1,:ng)=real(matmul(Ry(:ng,:ng,i,j,k),savg(ind1,:ng)),4)
        if(jyp /= 0) savg(4*ll4f+ind1,:ng)=savg(4*ll4f+ind1,:ng)+ &
          & real(matmul(Ry(:ng,ng+1:2*ng,i,j,k),savg(7*ll4f+ll4x+jyp,:ng)),4)
        !
        if(jzm_m /= 0) savg(4*ll4f+ind1,:ng)=savg(4*ll4f+ind1,:ng)+ &
          & real(matmul(Ry(:ng,2*ng+1:3*ng,i,j,k),savg(7*ll4f+ll4x+ll4y+jzm_m,:ng)),4)
        if(jzm /= 0) savg(4*ll4f+ind1,:ng)=savg(4*ll4f+ind1,:ng)+ &
          & real(matmul(Ry(:ng,3*ng+1:4*ng,i,j,k),savg(7*ll4f+ll4x+ll4y+jzm,:ng)),4)
        if(jzm_p /= 0) savg(4*ll4f+ind1,:ng)=savg(4*ll4f+ind1,:ng)+ &
          & real(matmul(Ry(:ng,4*ng+1:5*ng,i,j,k),savg(7*ll4f+ll4x+ll4y+jzm_p,:ng)),4)
        if(jzp_m /= 0) savg(4*ll4f+ind1,:ng)=savg(4*ll4f+ind1,:ng)+ &
          & real(matmul(Ry(:ng,5*ng+1:6*ng,i,j,k),savg(7*ll4f+ll4x+ll4y+jzp_m,:ng)),4)
        if(jzp /= 0) savg(4*ll4f+ind1,:ng)=savg(4*ll4f+ind1,:ng)+ &
          & real(matmul(Ry(:ng,6*ng+1:7*ng,i,j,k),savg(7*ll4f+ll4x+ll4y+jzp,:ng)),4)
        if(jzp_p /= 0) savg(4*ll4f+ind1,:ng)=savg(4*ll4f+ind1,:ng)+ &
          & real(matmul(Ry(:ng,7*ng+1:8*ng,i,j,k),savg(7*ll4f+ll4x+ll4y+jzp_p,:ng)),4)
        !
        if(jxm_m /= 0) savg(4*ll4f+ind1,:ng)=savg(4*ll4f+ind1,:ng)+ &
          & real(matmul(Ry(:ng,8*ng+1:9*ng,i,j,k),savg(7*ll4f+jxm_m,:ng)),4)
        if(jxm /= 0) savg(4*ll4f+ind1,:ng)=savg(4*ll4f+ind1,:ng)+ &
          & real(matmul(Ry(:ng,9*ng+1:10*ng,i,j,k),savg(7*ll4f+jxm,:ng)),4)
        if(jxm_p /= 0) savg(4*ll4f+ind1,:ng)=savg(4*ll4f+ind1,:ng)+ &
          & real(matmul(Ry(:ng,10*ng+1:11*ng,i,j,k),savg(7*ll4f+jxm_p,:ng)),4)
        if(jxp_m /= 0) savg(4*ll4f+ind1,:ng)=savg(4*ll4f+ind1,:ng)+ &
          & real(matmul(Ry(:ng,11*ng+1:12*ng,i,j,k),savg(7*ll4f+jxp_m,:ng)),4)
        if(jxp /= 0) savg(4*ll4f+ind1,:ng)=savg(4*ll4f+ind1,:ng)+ &
          & real(matmul(Ry(:ng,12*ng+1:13*ng,i,j,k),savg(7*ll4f+jxp,:ng)),4)
        if(jxp_p /= 0) savg(4*ll4f+ind1,:ng)=savg(4*ll4f+ind1,:ng)+ &
          & real(matmul(Ry(:ng,13*ng+1:14*ng,i,j,k),savg(7*ll4f+jxp_p,:ng)),4)
        !
        jxm_m=0 ; jxp_m=0 ; jxm_p=0 ; jxp_p=0
        jym_m=0 ; jyp_m=0 ; jym_p=0 ; jyp_p=0
        if((k == 1).and.(iqfr(5,i,j,1) == -2)) then
          jxm_m=kn(1,i,j,1) ; jxp_m=kn(2,i,j,1)
          jym_m=kn(3,i,j,1) ; jyp_m=kn(4,i,j,1)
        else if(k > 1) then
          jxm_m=kn(1,i,j,k-1) ; jxp_m=kn(2,i,j,k-1)
          jym_m=kn(3,i,j,k-1) ; jyp_m=kn(4,i,j,k-1)
        endif
        if((k == nz).and.(iqfr(6,i,j,nz) == -2)) then
          jxm_p=kn(1,i,j,nz) ; jxp_p=kn(2,i,j,nz)
          jym_p=kn(3,i,j,nz) ; jyp_p=kn(4,i,j,nz)
        else if(k < nz) then
          jxm_p=kn(1,i,j,k+1) ; jxp_p=kn(2,i,j,k+1)
          jym_p=kn(3,i,j,k+1) ; jyp_p=kn(4,i,j,k+1)
        endif
        ! z- relations
        savg(5*ll4f+ind1,:ng)=real(matmul(Lz(:ng,:ng,i,j,k),savg(ind1,:ng)),4)
        if(jzm /= 0) savg(5*ll4f+ind1,:ng)=savg(5*ll4f+ind1,:ng)+ &
          & real(matmul(Lz(:ng,ng+1:2*ng,i,j,k),savg(7*ll4f+ll4x+ll4y+jzm,:ng)),4)
        !
        if(jxm_m /= 0) savg(5*ll4f+ind1,:ng)=savg(5*ll4f+ind1,:ng)+ &
          & real(matmul(Lz(:ng,2*ng+1:3*ng,i,j,k),savg(7*ll4f+jxm_m,:ng)),4)
        if(jxm /= 0) savg(5*ll4f+ind1,:ng)=savg(5*ll4f+ind1,:ng)+ &
          & real(matmul(Lz(:ng,3*ng+1:4*ng,i,j,k),savg(7*ll4f+jxm,:ng)),4)
        if(jxm_p /= 0) savg(5*ll4f+ind1,:ng)=savg(5*ll4f+ind1,:ng)+ &
          & real(matmul(Lz(:ng,4*ng+1:5*ng,i,j,k),savg(7*ll4f+jxm_p,:ng)),4)
        if(jxp_m /= 0) savg(5*ll4f+ind1,:ng)=savg(5*ll4f+ind1,:ng)+ &
          & real(matmul(Lz(:ng,5*ng+1:6*ng,i,j,k),savg(7*ll4f+jxp_m,:ng)),4)
        if(jxp /= 0) savg(5*ll4f+ind1,:ng)=savg(5*ll4f+ind1,:ng)+ &
          & real(matmul(Lz(:ng,6*ng+1:7*ng,i,j,k),savg(7*ll4f+jxp,:ng)),4)
        if(jxp_p /= 0) savg(5*ll4f+ind1,:ng)=savg(5*ll4f+ind1,:ng)+ &
          & real(matmul(Lz(:ng,7*ng+1:8*ng,i,j,k),savg(7*ll4f+jxp_p,:ng)),4)
        !
        if(jym_m /= 0) savg(5*ll4f+ind1,:ng)=savg(5*ll4f+ind1,:ng)+ &
          & real(matmul(Lz(:ng,8*ng+1:9*ng,i,j,k),savg(7*ll4f+ll4x+jym_m,:ng)),4)
        if(jym /= 0) savg(5*ll4f+ind1,:ng)=savg(5*ll4f+ind1,:ng)+ &
          & real(matmul(Lz(:ng,9*ng+1:10*ng,i,j,k),savg(7*ll4f+ll4x+jym,:ng)),4)
        if(jym_p /= 0) savg(5*ll4f+ind1,:ng)=savg(5*ll4f+ind1,:ng)+ &
          & real(matmul(Lz(:ng,10*ng+1:11*ng,i,j,k),savg(7*ll4f+ll4x+jym_p,:ng)),4)
        if(jyp_m /= 0) savg(5*ll4f+ind1,:ng)=savg(5*ll4f+ind1,:ng)+ &
          & real(matmul(Lz(:ng,11*ng+1:12*ng,i,j,k),savg(7*ll4f+ll4x+jyp_m,:ng)),4)
        if(jyp /= 0) savg(5*ll4f+ind1,:ng)=savg(5*ll4f+ind1,:ng)+ &
          & real(matmul(Lz(:ng,12*ng+1:13*ng,i,j,k),savg(7*ll4f+ll4x+jyp,:ng)),4)
        if(jyp_p /= 0) savg(5*ll4f+ind1,:ng)=savg(5*ll4f+ind1,:ng)+ &
          & real(matmul(Lz(:ng,13*ng+1:14*ng,i,j,k),savg(7*ll4f+ll4x+jyp_p,:ng)),4)
        !
        ! z+ relations
        savg(6*ll4f+ind1,:ng)=real(matmul(Rz(:ng,:ng,i,j,k),savg(ind1,:ng)),4)
        if(jzp /= 0) savg(6*ll4f+ind1,:ng)=savg(6*ll4f+ind1,:ng)+ &
          & real(matmul(Rz(:ng,ng+1:2*ng,i,j,k),savg(7*ll4f+ll4x+ll4y+jzp,:ng)),4)
        !
        if(jxm_m /= 0) savg(6*ll4f+ind1,:ng)=savg(6*ll4f+ind1,:ng)+ &
          & real(matmul(Rz(:ng,2*ng+1:3*ng,i,j,k),savg(7*ll4f+jxm_m,:ng)),4)
        if(jxm /= 0) savg(6*ll4f+ind1,:ng)=savg(6*ll4f+ind1,:ng)+ &
          & real(matmul(Rz(:ng,3*ng+1:4*ng,i,j,k),savg(7*ll4f+jxm,:ng)),4)
        if(jxm_p /= 0) savg(6*ll4f+ind1,:ng)=savg(6*ll4f+ind1,:ng)+ &
          & real(matmul(Rz(:ng,4*ng+1:5*ng,i,j,k),savg(7*ll4f+jxm_p,:ng)),4)
        if(jxp_m /= 0) savg(6*ll4f+ind1,:ng)=savg(6*ll4f+ind1,:ng)+ &
          & real(matmul(Rz(:ng,5*ng+1:6*ng,i,j,k),savg(7*ll4f+jxp_m,:ng)),4)
        if(jxp /= 0) savg(6*ll4f+ind1,:ng)=savg(6*ll4f+ind1,:ng)+ &
          & real(matmul(Rz(:ng,6*ng+1:7*ng,i,j,k),savg(7*ll4f+jxp,:ng)),4)
        if(jxp_p /= 0) savg(6*ll4f+ind1,:ng)=savg(6*ll4f+ind1,:ng)+ &
          & real(matmul(Rz(:ng,7*ng+1:8*ng,i,j,k),savg(7*ll4f+jxp_p,:ng)),4)
        !
        if(jym_m /= 0) savg(6*ll4f+ind1,:ng)=savg(6*ll4f+ind1,:ng)+ &
          & real(matmul(Rz(:ng,8*ng+1:9*ng,i,j,k),savg(7*ll4f+ll4x+jym_m,:ng)),4)
        if(jym /= 0) savg(6*ll4f+ind1,:ng)=savg(6*ll4f+ind1,:ng)+ &
          & real(matmul(Rz(:ng,9*ng+1:10*ng,i,j,k),savg(7*ll4f+ll4x+jym,:ng)),4)
        if(jym_p /= 0) savg(6*ll4f+ind1,:ng)=savg(6*ll4f+ind1,:ng)+ &
          & real(matmul(Rz(:ng,10*ng+1:11*ng,i,j,k),savg(7*ll4f+ll4x+jym_p,:ng)),4)
        if(jyp_m /= 0) savg(6*ll4f+ind1,:ng)=savg(6*ll4f+ind1,:ng)+ &
          & real(matmul(Rz(:ng,11*ng+1:12*ng,i,j,k),savg(7*ll4f+ll4x+jyp_m,:ng)),4)
        if(jyp /= 0) savg(6*ll4f+ind1,:ng)=savg(6*ll4f+ind1,:ng)+ &
          & real(matmul(Rz(:ng,12*ng+1:13*ng,i,j,k),savg(7*ll4f+ll4x+jyp,:ng)),4)
        if(jyp_p /= 0) savg(6*ll4f+ind1,:ng)=savg(6*ll4f+ind1,:ng)+ &
          & real(matmul(Rz(:ng,13*ng+1:14*ng,i,j,k),savg(7*ll4f+ll4x+jyp_p,:ng)),4)
      enddo
    enddo
  enddo
  !----
  ! scratch storage deallocation
  !----
  deallocate(Rz,Lz,Ry,Ly,Rx,Lx)
  deallocate(work5,work4,work3,work2,work1)
  deallocate(LLR,Lambda,A)
end subroutine NSSANM3
