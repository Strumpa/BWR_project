subroutine NSSMXYZ(ll4f,ndim,nx,ny,nz,nmix,mat,xx,yy,zz,idl,vol,iqfr,qfr, &
& diff,drift,sigt,mux,muy,muz,imax,imay,imaz,ipy,ipz,a11x,a11y,a11z)
!
!-----------------------------------------------------------------------
!
!Purpose:
! Assembly of system matrices for coarse mesh finite differences with
! nodal correction.
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
! ll4f    total number of averaged flux unknown per energy group.
! ndim    number of dimensions (1, 2, or 3).
! nx      number of nodes in the X direction.
! ny      number of nodes in the Y direction.
! nz      number of nodes in the Z direction.
! nmix    number of mixtures.
! mat     node mixtures.
! xx      node widths in the X direction.
! yy      node widths in the Y direction.
! zz      node widths in the Z direction.
! idl     position of averaged fluxes in unknown vector.
! vol     node volumes.
! iqfr    boundary conditions.
! qfr     albedo functions.
! diff    diffusion coefficients.
! drift   drift coefficients.
! sigt    removal macroscopic cross section.
! mux     X-oriented compressed storage mode indices.
! muy     Y-oriented compressed storage mode indices.
! muz     Z-oriented compressed storage mode indices.
! imax    X-oriented position of each first non-zero column element.
! imay    Y-oriented position of each first non-zero column element.
! imaz    Z-oriented position of each first non-zero column element.
! ipy     Y-oriented permutation matrices.
! ipz     Z-oriented permutation matrices.
!
!Parameters: output
! a11x    X-directed  matrices corresponding to the divergence (i.e.
!         leakage) and removal terms. Dimensionned to imax(ll4f).
! a11y    Y-directed  matrices corresponding to the divergence (i.e.
!         leakage) and removal terms. Dimensionned to imay(ll4f).
! a11z    Z-directed  matrices corresponding to the divergence (i.e.
!         leakage) and removal terms. Dimensionned to imaz(ll4f).
!
!-----------------------------------------------------------------------
!
  !----
  !  subroutine arguments
  !----
  integer,intent(in) :: ll4f,ndim,nx,ny,nz,nmix,mat(nx,ny,nz),idl(nx,ny,nz), &
  & iqfr(6,nx,ny,nz),mux(ll4f),muy(ll4f),muz(ll4f),imax(ll4f),imay(ll4f),imaz(ll4f), &
  & ipy(ll4f),ipz(ll4f)
  real,intent(in) :: xx(nx,ny,nz),yy(nx,ny,nz),zz(nx,ny,nz),vol(nx,ny,nz), &
  & qfr(6,nx,ny,nz),diff(nmix),drift(6,nx,ny,nz),sigt(nmix)
  real,intent(out) :: a11x(*),a11y(*),a11z(*)
  !----
  !  local variables
  !----
  real :: coef(6),codr(6)
  !
  a11x(:imax(ll4f))=0.0
  if(ndim > 1) a11y(:imay(ll4f))=0.0
  if(ndim == 3) a11z(:imaz(ll4f))=0.0
  do k=1,nz
    do j=1,ny
      do i=1,nx
        ibm=mat(i,j,k)
        if(ibm <= 0) cycle
        kel=idl(i,j,k)
        if(kel == 0) cycle
        vol0=vol(i,j,k)
        call NSSCO(nx,ny,nz,nmix,i,j,k,mat,xx,yy,zz,diff,iqfr(1,i,j,k),qfr(1,i,j,k),coef)
        coef(1:2)=coef(1:2)*vol0/xx(i,j,k)
        coef(3:4)=coef(3:4)*vol0/yy(i,j,k)
        coef(5:6)=coef(5:6)*vol0/zz(i,j,k)
        codr(1:2)=drift(1:2,i,j,k)*vol0/xx(i,j,k)
        codr(3:4)=drift(3:4,i,j,k)*vol0/yy(i,j,k)
        codr(5:6)=drift(5:6,i,j,k)*vol0/zz(i,j,k)
        !
        ! x-directed couplings
        kel2=0
        kk1=iqfr(1,i,j,k)
        if(kk1 == -4) then
          kel2=idl(nx,j,k)
        else if(kk1 == 0) then
          kel2=idl(i-1,j,k)
        endif
        if(kel2 /= 0) then
          if(kel2 <= kel) then
            key=mux(kel)-kel+kel2
            a11x(key)=a11x(key)-coef(1)+codr(1)
          else
            key=mux(kel2)+kel2-kel
            a11x(key)=a11x(key)-coef(1)+codr(1)
          endif
        endif
        kel2=0
        kk2=iqfr(2,i,j,k)
        if(kk2 == -4) then
          kel2=idl(1,j,k)
        else if(kk2 == 0) then
          kel2=idl(i+1,j,k)
        endif
        if(kel2 /= 0) then
          if(kel2 <= kel) then
            key=mux(kel)-kel+kel2
            a11x(key)=a11x(key)-coef(2)-codr(2)
          else
            key=mux(kel2)+kel2-kel
            a11x(key)=a11x(key)-coef(2)-codr(2)
          endif
        endif
        key0=mux(kel)
        a11x(key0)=a11x(key0)+coef(1)+codr(1)+coef(2)-codr(2)
        a11x(key0)=a11x(key0)+coef(3)+codr(3)+coef(4)-codr(4)
        a11x(key0)=a11x(key0)+coef(5)+codr(5)+coef(6)-codr(6)
        a11x(key0)=a11x(key0)+sigt(ibm)*vol0
        !
        if(ndim > 1) then
          ! y-directed couplings
          kel2=0
          kk3=iqfr(3,i,j,k)
          if(kk3 == -4) then
            kel2=idl(i,ny,k)
          else if(kk3 == 0) then
            kel2=idl(i,j-1,k)
          endif
          ind1=ipy(kel)
          if(kel2 /= 0) then
            ind2=ipy(kel2)
            if(kel2 <= kel) then
              key=muy(ind1)-ind1+ind2
              a11y(key)=a11y(key)-coef(3)+codr(3)
            else
              key=muy(ind2)+ind2-ind1
              a11y(key)=a11y(key)-coef(3)+codr(3)
            endif
          endif
          kel2=0
          kk4=iqfr(4,i,j,k)
          if(kk4 == -4) then
            kel2=idl(i,1,k)
          else if(kk4 == 0) then
            kel2=idl(i,j+1,k)
          endif
          if(kel2 /= 0) then
            ind2=ipy(kel2)
            if(kel2 <= kel) then
              key=muy(ind1)-ind1+ind2
              a11y(key)=a11y(key)-coef(4)-codr(4)
            else
              key=muy(ind2)+ind2-ind1
              a11y(key)=a11y(key)-coef(4)-codr(4)
            endif
          endif
          key0=muy(ind1)
          a11y(key0)=a11y(key0)+coef(1)+codr(1)+coef(2)-codr(2)
          a11y(key0)=a11y(key0)+coef(3)+codr(3)+coef(4)-codr(4)
          a11y(key0)=a11y(key0)+coef(5)+codr(5)+coef(6)-codr(6)
          a11y(key0)=a11y(key0)+sigt(ibm)*vol0
        endif
        !
        if(ndim > 2) then
          ! z-directed couplings
          kel2=0
          kk5=iqfr(5,i,j,k)
          if(kk5 == -4) then
            kel2=idl(i,j,nz)
          else if(kk5 == 0) then
            kel2=idl(i,j,k-1)
          endif
          ind1=ipz(kel)
          if(kel2 /= 0) then
            ind2=ipz(kel2)
            if(kel2 <= kel) then
              key=muz(ind1)-ind1+ind2
              a11z(key)=a11z(key)-coef(5)+codr(5)
            else
              key=muz(ind2)+ind2-ind1
              a11z(key)=a11z(key)-coef(5)+codr(5)
            endif
          endif
          kel2=0
          kk6=iqfr(6,i,j,k)
          if(kk6 == -4) then
            kel2=idl(i,j,1)
          else if(kk6 == 0) then
            kel2=idl(i,j,k+1)
          endif
          if(kel2 /= 0) then
            ind2=ipz(kel2)
            if(kel2 <= kel) then
              key=muz(ind1)-ind1+ind2
              a11z(key)=a11z(key)-coef(6)-codr(6)
            else
              key=muz(ind2)+ind2-ind1
              a11z(key)=a11z(key)-coef(6)-codr(6)
            endif
          endif
          key0=muz(ind1)
          a11z(key0)=a11z(key0)+coef(1)+codr(1)+coef(2)-codr(2)
          a11z(key0)=a11z(key0)+coef(3)+codr(3)+coef(4)-codr(4)
          a11z(key0)=a11z(key0)+coef(5)+codr(5)+coef(6)-codr(6)
          a11z(key0)=a11z(key0)+sigt(ibm)*vol0
        endif
      enddo
    enddo
  enddo
  return
end subroutine NSSMXYZ
