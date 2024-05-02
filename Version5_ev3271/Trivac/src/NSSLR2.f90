subroutine NSSLR2(keff, ng, bndtl, xxx, dely, diff, sigr, scat, chi, nusigf, L, R)
!
!-----------------------------------------------------------------------
!
!Purpose:
! Compute the 2D ANM coupling matrices for a single node.
!
!Copyright:
! Copyright (C) 2023 Ecole Polytechnique de Montreal
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
! bndtl   set to 'flat' or 'quadratic'.
! xxx     node support along X-axis.
! dely    node width along Y-axis.
! diff    diffusion coefficient array (cm).
! sigr    removal cross section array (cm-1).
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
  real, intent(in) :: keff, xxx(4), dely
  real, dimension(ng), intent(in) :: diff, sigr, chi, nusigf
  real, dimension(ng,ng), intent(in) :: scat
  character(len=12), intent(in) :: bndtl
  real(kind=8), dimension(ng,8*ng), intent(out) :: L, R
  !----
  !  local variables
  !----
  real(kind=8) :: m0(3,3),m2(3,3),m3(2,3),m4(1,3),Lambda_r,sqla,mmax2
  !----
  ! allocatable arrays
  !----
  complex(kind=8), allocatable, dimension(:,:) :: T,Lambda
  real(kind=8), allocatable, dimension(:,:) :: F,DI,T_r,TI,S,Mm,Mp,Nm,Np, &
  & GAR1,GAR2,GAR3,GAR4,Vm,Vp,Um,Up,MAT1,MAT2,S7
  !----
  ! scratch storage allocation
  !----
  allocate(F(ng,ng),T_r(ng,ng),T(ng,ng),TI(ng,ng),DI(ng,ng), &
  & Lambda(ng,ng),S(ng,ng),Mm(2*ng,2*ng),Mp(2*ng,2*ng),Nm(ng,2*ng), &
  & Np(ng,2*ng),GAR1(ng,2*ng),GAR2(ng,2*ng),GAR3(ng,8*ng), &
  & GAR4(ng,8*ng),Vm(2*ng,3*ng),Vp(2*ng,3*ng),Um(ng,3*ng), &
  & Up(ng,3*ng),MAT1(ng,8*ng),MAT2(ng,8*ng))
  !
  ! quadratic leakage and boundary conditions
  xmm=xxx(1) ; xm=xxx(2) ; xp=xxx(3) ; xpp=xxx(4) ; delx=xp-xm ;
  if(xmm == -99999.) then
    ! Vacuum or zero flux node at left boundary
    xmm=2.0*xm-xp
    m0(:3,1)=1.0d0 ; m0(1,2)=(xmm+xm)/2.0d0 ; m0(1,3)=(xmm**2+xmm*xm+xm**2)/3.0d0
    m0(2,2)=(xm+xp)/2.0d0 ; m0(2,3)=(xm**2+xm*xp+xp**2)/3.0d0
    m0(3,2)=(xp+xpp)/2.0d0 ; m0(3,3)=(xp**2+xp*xpp+xpp**2)/3.0d0
    call ALINVD(3,m0,3,ier)
    if(ier /= 0) call XABORT('NSSLR2: singular matrix.(1)')
    m0(:3,1)=0.0d0
  elseif(xpp == -99999.) then
    ! Vacuum or zero flux node at right boundary
    xpp=2.0*xp-xm
    m0(:3,1)=1.0d0 ; m0(1,2)=(xmm+xm)/2.0d0 ; m0(1,3)=(xmm**2+xmm*xm+xm**2)/3.0d0
    m0(2,2)=(xm+xp)/2.0d0 ; m0(2,3)=(xm**2+xm*xp+xp**2)/3.0d0
    m0(3,2)=(xp+xpp)/2.0d0 ; m0(3,3)=(xp**2+xp*xpp+xpp**2)/3.0d0
    call ALINVD(3,m0,3,ier)
    if(ier /= 0) call XABORT('NSSLR2: singular matrix.(2)')
    m0(:3,3)=0.0d0
  else
    ! Internal node
    m0(:3,1)=1.0d0 ; m0(1,2)=(xmm+xm)/2.0d0 ; m0(1,3)=(xmm**2+xmm*xm+xm**2)/3.0d0
    m0(2,2)=(xm+xp)/2.0d0 ; m0(2,3)=(xm**2+xm*xp+xp**2)/3.0d0
    m0(3,2)=(xp+xpp)/2.0d0 ; m0(3,3)=(xp**2+xp*xpp+xpp**2)/3.0d0
    call ALINVD(3,m0,3,ier)
    if(ier /= 0) call XABORT('NSSLR2: singular matrix.(3)')
  endif
  if(bndtl == 'flat') then
    ! flat leakage approximation
    m0(:3,:3)=0.0d0 ; m0(1,2)=1.0d0
  endif
  !----
  ! compute matrices L and R
  !----
  Mm(:,:)=0.0d0
  Mp(:,:)=0.0d0
  Nm(:,:)=0.0d0
  Np(:,:)=0.0d0
  DI(:,:)=0.0d0
  Vm(:,:)=0.0d0
  Vp(:,:)=0.0d0
  Um(:,:)=0.0d0
  Up(:,:)=0.0d0
  do ig=1,ng
    do jg=1,ng
      if(ig == jg) then
        F(ig,ig)=(chi(ig)*nusigf(ig)/keff-sigr(ig))/diff(ig)
      else
        F(ig,jg)=(chi(ig)*nusigf(jg)/keff+scat(ig,jg))/diff(ig)
      endif
    enddo
    DI(ig,ig)=1./diff(ig)
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
    call XABORT('NSSLR2: complex eigenvalues.')
  endif
  T_r(:,:)=real(T(:,:),8)
  do ig=1,ng
    Lambda_r=real(Lambda(ig,ig),8)
    sqla=sqrt(abs(Lambda_r))
    m2(:3,:3)=0.0d0
    m2(1,1)=1.0d0/Lambda_r ; m2(1,3)=-2.0d0/Lambda_r**2
    m2(2,2)=1.0d0/Lambda_r ; m2(3,3)=1.0d0/Lambda_r
    m2(:3,:3)=matmul(m2(:3,:3),m0(:3,:3))
    m3(1,1)=1.0d0 ; m3(1,2)=(xm+xp)/2. ; m3(1,3)=(xm**2+xm*xp+xp**2)/3.0d0
    m3(2,1)=0.0d0 ; m3(2,2)=-1.0d0 ; m3(2,3)=-2.0d0*xm
    m3(:2,:3)=matmul(m3(:2,:3),m2(:3,:3))
    Vm(ig,ig)=m3(1,1)   ; Vm(ig,ng+ig)=m3(1,2)   ; Vm(ig,2*ng+ig)=m3(1,3)   ;
    Vm(ng+ig,ig)=m3(2,1) ; Vm(ng+ig,ng+ig)=m3(2,2) ; Vm(ng+ig,2*ng+ig)=m3(2,3) ;
    m3(1,1)=1.0d0 ; m3(1,2)=(xm+xp)/2.0d0 ; m3(1,3)=(xm**2+xm*xp+xp**2)/3.0d0
    m3(2,1)=0.0d0 ; m3(2,2)=-1.0d0 ; m3(2,3)=-2.0d0*xp
    m3(:2,:3)=matmul(m3(:2,:3),m2(:3,:3))
    Vp(ig,ig)=m3(1,1)   ; Vp(ig,ng+ig)=m3(1,2)   ; Vp(ig,2*ng+ig)=m3(1,3)   ;
    Vp(ng+ig,ig)=m3(2,1) ; Vp(ng+ig,ng+ig)=m3(2,2) ; Vp(ng+ig,2*ng+ig)=m3(2,3) ;
    m4(1,1)=1.0d0 ; m4(1,2)=xm ; m4(1,3)=xm**2
    m4(:1,:3)=matmul(m4(:1,:3),m2(:3,:3))
    Um(ig,ig)=m4(1,1) ; Um(ig,ng+ig)=m4(1,2) ; Um(ig,2*ng+ig)=m4(1,3) ;
    m4(1,1)=1.0d0 ; m4(1,2)=xp ; m4(1,3)=xp**2
    m4(:1,:3)=matmul(m4(:1,:3),m2(:3,:3))
    Up(ig,ig)=m4(1,1) ; Up(ig,ng+ig)=m4(1,2) ; Up(ig,2*ng+ig)=m4(1,3) ;
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
  if(ier /= 0) call XABORT('NSSLR2: singular matrix.(4)')
  call ALINVD(2*ng,Mp,2*ng,ier)
  if(ier /= 0) call XABORT('NSSLR2: singular matrix.(5)')
  call ALINVD(ng,TI,ng,ier)
  if(ier /= 0) call XABORT('NSSLR2: singular matrix.(6)')
  !
  GAR1=matmul(Nm,Mm)  ! ng,2*ng
  GAR2=matmul(Np,Mp)  ! ng,2*ng
  S=matmul(TI,DI)     ! ng,ng
  !
  MAT1(:ng,:2*ng)=GAR1(:ng,:2*ng)
  MAT1(:ng,2*ng+1:5*ng)=-Um(:ng,:3*ng)/dely+matmul(GAR1(:ng,:2*ng),Vm(:2*ng,:3*ng))/dely
  MAT1(:ng,5*ng+1:8*ng)=Um(:ng,:3*ng)/dely-matmul(GAR1(:ng,:2*ng),Vm(:2*ng,:3*ng))/dely
  MAT2(:ng,:2*ng)=GAR2(:ng,:2*ng)
  MAT2(:ng,2*ng+1:5*ng)=-Up(:ng,:3*ng)/dely+matmul(GAR2(:ng,:2*ng),Vp(:2*ng,:3*ng))/dely
  MAT2(:ng,5*ng+1:8*ng)=Up(:ng,:3*ng)/dely-matmul(GAR2(:ng,:2*ng),Vp(:2*ng,:3*ng))/dely
  !
  GAR3=matmul(T_r,MAT1) ! ng,8*ng
  GAR4=matmul(T_r,MAT2) ! ng,8*ng
  L(:ng,:ng)=matmul(GAR3(:ng,:ng),TI(:ng,:ng))
  R(:ng,:ng)=matmul(GAR4(:ng,:ng),TI(:ng,:ng))
  allocate(S7(7*ng,7*ng))
  S7(:,:)=0.0d0 ! 7*ng,7*ng
  do i=1,7
    S7((i-1)*ng+1:i*ng,(i-1)*ng+1:i*ng)=S(:ng,:ng)
  enddo
  L(:ng,ng+1:8*ng)=matmul(GAR3(:ng,ng+1:8*ng),S7(:7*ng,:7*ng))
  R(:ng,ng+1:8*ng)=matmul(GAR4(:ng,ng+1:8*ng),S7(:7*ng,:7*ng))
  !----
  ! scratch storage deallocation
  !----
  deallocate(S7,MAT2,MAT1,Up,Um,Vp,Vm,GAR4,GAR3,GAR2,GAR1,Np,Nm,Mp,Mm,S, &
  & Lambda,DI,TI,T,T_r,F)
end subroutine NSSLR2
