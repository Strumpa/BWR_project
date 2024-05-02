subroutine NSSANM1(nel,ng,nmix,iqfr,qfr,mat,xxx,keff,diff,sigr,chi,sigf,scat,fd,savg)
!
!-----------------------------------------------------------------------
!
!Purpose:
! Compute the ANM volume fluxes and boundary fluxes and currents using
! a solution of one- and two-node relations in Cartesian 1D geometry.
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
! nel     number of nodes in the nodal calculation.
! ng      number of energy groups.
! nmix    number of material mixtures in the nodal calculation.
! iqfr    node-ordered physical albedo indices.
! qfr     albedo function information.
! mat     material mixture index in eacn node.
! xxx     Cartesian coordinates along the X axis.
! keff    effective multiplication facctor.
! diff    diffusion coefficients
! sigr    removal cross sections.
! chi     fission spectra.
! sigf    nu times fission cross section.
! scat    scattering cross section.
! fd      discontinuity factors
! savg    nodal fluxes.
!
!Parameters: output
! savg    boundary fluxes and currents.
!
!-----------------------------------------------------------------------
  !
  !----
  !  subroutine arguments
  !----
  integer,intent(in) :: nel,ng,nmix,iqfr(6,nel),mat(nel)
  real,intent(in) :: qfr(6,nel,ng),xxx(nel+1),keff,diff(nmix,ng),sigr(nmix,ng), &
  & chi(nmix,ng),sigf(nmix,ng),scat(nmix,ng,ng),fd(nmix,2,ng,ng)
  real, dimension(4*nel+1,ng),intent(inout) :: savg
  !----
  ! allocatable arrays
  !----
  real, allocatable, dimension(:) :: work1,work2,work4,work5
  real, allocatable, dimension(:,:) :: A,B,Lambda,work3
  real(kind=8), allocatable, dimension(:,:,:) :: Lx,Rx
  !----
  ! scratch storage allocation
  !----
  allocate(A(ng,ng+1),B(ng,ng),Lambda(ng,ng))
  allocate(work1(ng),work2(ng),work3(ng,ng),work4(ng),work5(ng))
  allocate(Lx(ng,2*ng,nel),Rx(ng,2*ng,nel))
  !
  ! compute nodal coefficients
  do iel=1,nel
    ibm=mat(iel)
    if(ibm == 0) cycle
    work1(:ng)=diff(ibm,:ng)
    work2(:ng)=sigr(ibm,:ng)
    work3(:ng,:ng)=scat(ibm,:ng,:ng)
    work4(:ng)=chi(ibm,:ng)
    work5(:ng)=sigf(ibm,:ng)
    delx=xxx(iel+1)-xxx(iel)
    call NSSLR1(keff,ng,delx,work1,work2,work3,work4,work5, &
    & Lx(1,1,iel),Rx(1,1,iel))
  enddo
  !----
  ! compute boundary currents
  ! left one-node relation
  !----
  A(:ng,:ng+1)=0.0
  if((iqfr(1,1) > 0).or.(iqfr(1,1) == -1)) then
    ! physical albedo
    Lambda(:ng,:ng)=0.0
    do ig=1,ng
      Lambda(ig,ig)=qfr(1,1,ig)
    enddo
    A(:ng,:ng)=real(matmul(Lambda(:ng,:ng),Lx(:ng,ng+1:2*ng,1)),4)
    B(:ng,:ng)=real(matmul(Lambda(:ng,:ng),Lx(:ng,:ng,1)),4)
    do ig=1,ng
      A(ig,ig)=1.0+A(ig,ig)
    enddo
    A(:ng,ng+1)=-matmul(B(:ng,:ng),savg(1,:ng))
  else if(iqfr(1,1) == -2) then
    ! zero net current
    do ig=1,ng
      A(ig,ig)=1.0
    enddo
  else if(iqfr(1,1) == -3) then
    ! zero flux
    A(:ng,:ng)=real(Lx(:ng,ng+1:2*ng,1),4)
    A(:ng,ng+1)=real(-matmul(Lx(:ng,:ng,1),savg(1,:ng)),4)
  else
    call XABORT('NSSANM1: illegal left boundary condition.')
  endif
  call ALSB(ng,1,A,ier,ng)
  if(ier /= 0) call XABORT('NSSANM1: singular matrix.(1)')
  savg(3*nel+1,:ng)=A(:ng,ng+1)
  ! two-node relations
  do i=2,nel
    A(:ng,:ng)=real(matmul(fd(mat(i-1),2,:ng,:ng),Rx(:ng,ng+1:2*ng,i-1))- &
                  & matmul(fd(mat(i),1,:ng,:ng),Lx(:ng,ng+1:2*ng,i)),4)
    A(:ng,ng+1)=-real(matmul(matmul(fd(mat(i-1),2,:ng,:ng),Rx(:ng,:ng,i-1)),savg(i-1,:ng))- &
                  & matmul(matmul(fd(mat(i),1,:ng,:ng),Lx(:ng,:ng,i)),savg(i,:ng)),4)
    call ALSB(ng,1,A,ier,ng)
    if(ier /= 0) call XABORT('NSSANM1: singular matrix.(2)')
    savg(3*nel+i,:ng)=A(:ng,ng+1)
  enddo
  ! right one-node relation
  if((iqfr(2,nel) > 0).or.(iqfr(2,nel) == -1)) then
    ! physical albedo
    Lambda(:ng,:ng)=0.0
    do ig=1,ng
      Lambda(ig,ig)=qfr(2,nel,ig)
    enddo
    A(:ng,:ng)=real(matmul(Lambda(:ng,:ng),Rx(:ng,ng+1:2*ng,nel)),4)
    B(:ng,:ng)=real(matmul(Lambda(:ng,:ng),Rx(:ng,:ng,nel)),4)
    do ig=1,ng
      A(ig,ig)=-1.0+A(ig,ig)
    enddo
    A(:ng,ng+1)=-matmul(B(:ng,:ng),savg(nel,:ng))
  else if(iqfr(2,nel) == -2) then
    ! zero net current
    do ig=1,ng
      A(2*nel*ng+ig,2*nel*ng+ig)=1.0
    enddo
  else if(iqfr(2,nel) == -3) then
    ! zero flux
    A(:ng,:ng)=real(Rx(:ng,ng+1:2*ng,nel),4)
    A(:ng,ng+1)=real(-matmul(Rx(:ng,:ng,nel),savg(nel,:ng)),4)
  else
    call XABORT('NSSANM1: illegal right boundary condition.')
  endif
  call ALSB(ng,1,A,ier,ng)
  if(ier /= 0) call XABORT('NSSANM1: singular matrix.(3)')
  savg(4*nel+1,:ng)=A(:ng,ng+1)
  !----
  ! compute boundary fluxes
  !----
  do i=1,nel
    savg(nel+i,:ng)=real(matmul(Lx(:ng,:ng,i),savg(i,:ng))+ &
    & matmul(Lx(:ng,ng+1:2*ng,i),savg(3*nel+i,:ng)),4)
    savg(2*nel+i,:ng)=real(matmul(Rx(:ng,:ng,i),savg(i,:ng))+ &
    & matmul(Rx(:ng,ng+1:2*ng,i),savg(3*nel+i+1,:ng)),4)
  enddo
  !----
  ! scratch storage deallocation
  !----
  deallocate(Rx,Lx)
  deallocate(work5,work4,work3,work2,work1)
  deallocate(Lambda,B,A)
end subroutine NSSANM1
