subroutine NSSLR1(keff, ng, delx, diff, sigr, scat, chi, nusigf, L, R)
!
!-----------------------------------------------------------------------
!
!Purpose:
! Compute the 1D ANM coupling matrices for a single node.
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
! keff    effective multiplication factor.
! ng      number of energy groups.
! delx    node width along X-axis.
! diff    diffusion coefficient array (cm).
! sigr    removal cross section array (cm^-1).
! scat    P0 scattering cross section matrix (cm^-1).
! chi     fission spectrum array.
! nusigf  nu*fission cross section array (cm^-1).
!
!Parameters: output
! L       left nodal coupling matrix.
! R       right nodal coupling matrix.
!
!-----------------------------------------------------------------------
  !
  !----
  !  subroutine arguments
  !----
  integer, intent(in) :: ng
  real, intent(in) :: keff,delx
  real, dimension(ng), intent(in) :: diff, sigr, chi, nusigf
  real, dimension(ng,ng), intent(in) :: scat
  real(kind=8), dimension(ng,2*ng), intent(out) :: L,R
  !----
  !  local variables
  !----
  real(kind=8) :: Lambda_r,sqla,mmax2
  !----
  ! allocatable arrays
  !----
  complex(kind=8), allocatable, dimension(:,:) :: T,Lambda
  real(kind=8), allocatable, dimension(:,:) :: F,DI,T_r,TI,S,Mm, &
    & Mp,Nm,Np,GAR1,GAR2
  !----
  ! scratch storage allocation
  !----
  allocate(F(ng,ng),T_r(ng,ng),T(ng,ng),TI(ng,ng),DI(ng,ng), &
  & Lambda(ng,ng),S(ng,ng),Mm(2*ng,2*ng),Mp(2*ng,2*ng), &
  & Nm(ng,2*ng),Np(ng,2*ng),GAR1(ng,2*ng),GAR2(ng,2*ng))
  !----
  ! compute matrices L and R
  !----
  Mm(:,:)=0.0d0
  Mp(:,:)=0.0d0
  Nm(:,:)=0.0d0
  Np(:,:)=0.0d0
  DI(:,:)=0.0d0
  xm=0.0 ; xp=delx
  do ig=1,ng
    do jg=1,ng
      if(ig == jg) then
        F(ig,ig)=(chi(ig)*nusigf(ig)/keff-sigr(ig))/diff(ig)
      else
        F(ig,jg)=(chi(ig)*nusigf(jg)/keff+scat(ig,jg))/diff(ig)
      endif
    enddo
    DI(ig,ig)=1.d0/diff(ig)
  enddo
  maxiter=40
  call ALHQR(ng,ng,F,maxiter,iter,T,Lambda)
  mmax2=0.0d0
  do ig=1,ng
    do jg=1,ng
      mmax2=max(mmax2,abs(aimag(T(ig,jg))))
    enddo
  enddo
  if(mmax2 > 1.0e-6) then
    write(6,'(3h T=)')
    do ig=1,ng
      write(6,'(1p,12e12.4)') T(ig,:)
    enddo
    call XABORT('NSSLR1: complex eigenvalues.')
  endif
  T_r(:,:)=real(T(:,:),8)
  do ig=1,ng
    Lambda_r=real(Lambda(ig,ig),8)
    sqla=sqrt(abs(Lambda_r))
    if(delx*sqla < 1.e-6) then
      if(Lambda_r >= 0) then
        Mm(ig,ig)=-(delx*sqla)**6/5040.+(delx*sqla)**4/120.-(delx*sqla)**2/6.+1.
        Mm(ig,ng+ig)=(delx*sqla)**5/720.-(delx*sqla)**3/24.+(delx*sqla)/2.
        Mm(ng+ig,ng+ig)=-sqla
        Mp(ng+ig,ig)=((delx*sqla)**6/120.-(delx*sqla)**4/6.+(delx*sqla)**2)/delx
        Mp(ng+ig,ng+ig)=(-(delx*sqla)**5/24.+(delx*sqla)**3/2.-(delx*sqla))/delx
        Nm(ig,ig)=1.
        Np(ig,ig)=-(delx*sqla)**6/720.+(delx*sqla)**4/24.-(delx*sqla)**2/2.+1.
        Np(ig,ng+ig)=(delx*sqla)**5/120.-(delx*sqla)**3/6.+(delx*sqla)
      else
        Mm(ig,ig)=(delx*sqla)**4/120.+(delx*sqla)**3/24.+(delx*sqla)**2/6.+(delx*sqla)/2. + 1.
        Mm(ig,ng+ig)=-(delx*sqla)**3/24.+(delx*sqla)**2/6.-(delx*sqla)/2. + 1.
        Mm(ng+ig,ig)=-sqla ; Mm(ng+ig,ng+ig)=sqla ;
        Mp(ng+ig,ig)=(-(delx*sqla)**4/6.-(delx*sqla)**3/2.-(delx*sqla)**2-(delx*sqla))/delx
        Mp(ng+ig,ng+ig)=(-(delx*sqla)**4/6+(delx*sqla)**3/2.-(delx*sqla)**2+(delx*sqla))/delx
        Nm(ig,ig)=1. ; Nm(ig,ng+ig)=1. ;
        Np(ig,ig)=(delx*sqla)**4/24.+(delx*sqla)**3/6.+(delx*sqla)**2/2.+(delx*sqla)+1.
        Np(ig,ng+ig)=(delx*sqla)**4/24.-(delx*sqla)**3/6.+(delx*sqla)**2/2.-(delx*sqla)+1.
      endif
    else if(Lambda_r >= 0) then
      Mm(ig,ig)=(sin(sqla*xp)-sin(sqla*xm))/(delx*sqla)
      Mm(ig,ng+ig)=-(cos(sqla*xp)-cos(sqla*xm))/(delx*sqla)
      Mm(ng+ig,ig)=sqla*sin(sqla*xm)
      Mm(ng+ig,ng+ig)=-sqla*cos(sqla*xm)
      Mp(ng+ig,ig)=sqla*sin(sqla*xp)
      Mp(ng+ig,ng+ig)=-sqla*cos(sqla*xp)
      Nm(ig,ig)=cos(sqla*xm)
      Nm(ig,ng+ig)=sin(sqla*xm)
      Np(ig,ig)=cos(sqla*xp)
      Np(ig,ng+ig)=sin(sqla*xp)
    else
      Mm(ig,ig)=exp(sqla*xm)*(exp(sqla*(xp-xm))-1.0d0)/(delx*sqla)
      Mm(ig,ng+ig)=-exp(-sqla*xm)*(exp(-sqla*(xp-xm))-1.0d0)/(delx*sqla)
      Mm(ng+ig,ig)=-sqla*exp(sqla*xm)
      Mm(ng+ig,ng+ig)=sqla*exp(-sqla*xm)
      Mp(ng+ig,ig)=-sqla*exp(sqla*xp)
      Mp(ng+ig,ng+ig)=sqla*exp(-sqla*xp)
      Nm(ig,ig)=exp(sqla*xm)
      Nm(ig,ng+ig)=exp(-sqla*xm)
      Np(ig,ig)=exp(sqla*xp)
      Np(ig,ng+ig)=exp(-sqla*xp)
    endif
    Mp(ig,ig)=Mm(ig,ig)
    Mp(ig,ng+ig)=Mm(ig,ng+ig)
  enddo
  !
  TI(:,:)=T_r(:,:)
  call ALINVD(2*ng,Mm,2*ng,ier)
  if(ier /= 0) call XABORT('NSSLR1: singular matrix.(1)')
  call ALINVD(2*ng,Mp,2*ng,ier)
  if(ier /= 0) call XABORT('NSSLR1: singular matrix.(2)')
  call ALINVD(ng,TI,ng,ier)
  if(ier /= 0) call XABORT('NSSLR1: singular matrix.(3)')
  !
  GAR1=matmul(Nm,Mm)  ! ng,2*ng
  GAR2=matmul(Np,Mp)  ! ng,2*ng
  S=matmul(TI,DI)     ! ng,ng
  GAR1=matmul(T_r,GAR1) ! ng,2*ng
  GAR2=matmul(T_r,GAR2) ! ng,2*ng
  !
  L(:ng,:ng)=matmul(GAR1(:ng,:ng),TI(:ng,:ng))
  L(:ng,ng+1:2*ng)=matmul(GAR1(:ng,ng+1:2*ng),S(:ng,:ng))
  R(:ng,:ng)=matmul(GAR2(:ng,:ng),TI(:ng,:ng))
  R(:ng,ng+1:2*ng)=matmul(GAR2(:ng,ng+1:2*ng),S(:ng,:ng))
  !----
  ! scratch storage deallocation
  !----
  deallocate(GAR2,GAR1,Np,Nm,Mp,Mm,S,Lambda,DI,TI,T,T_r,F)
end subroutine NSSLR1
