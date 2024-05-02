!
!-----------------------------------------------------------------------
!
!Purpose:
! Fortran-2003 wrapper to calculate Molten Salt thermophysical 
! properties using data from the MSTPDB-TP Database
!
!Copyright:
! Copyright (C) 2023 Ecole Polytechnique de Montreal.
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
!Author(s): Cristian Garrido Tamm
!
!-----------------------------------------------------------------------
!
subroutine THMSGT(salt, compo, tp, impx)
    ! Get the data from MSTDB-TP csv files for the specific salt "salt" with proportions "prop"
    use t_saltdata
    character*12 :: salt, compo       ! salt formula and composition
    type(tpdata), intent(out) :: tp   ! Tupla with Thermophysical properties
    integer, intent(in) :: impx

    character(len=1000) :: line       ! line of the csv file
    character(len=255) :: filename = './MSTPDB.data'

    integer :: ios, i

    ! Initialize tp
    tp%formula=''
    tp%composition=''

    ! open file
    open(18, file=filename, status='old')

    ! Read file line by line
    i= 0
    do 
        i = i + 1
        read(18, '(a)', iostat=ios) line
        if (ios /= 0) exit
        if (impx > 5) write(*,*) "line=", line
        if (i .gt. 1) then
            read(line,*) tp
            if (impx > 5) write(*,*) "tp=", tp
        end if
        if (tp%formula .eq. salt .and. tp%composition .eq. compo) then 
            if (impx > 5) write(*,*) "Found!"
            exit
        end if
    end do
    if (ios /= 0) then
        if (impx > 5) write(*,*) "call xabort here"
        write(*,*) 'THMSGT: salt=', salt, 'compo=', compo
        call XABORT('Salt not found in MSTDB-TP')
    end if
    if (impx > 2) then
        write(*,*) "tp%formula=", tp%formula
        write(*,*) "tp%weight=", tp%weight
        write(*,*) "tp%composition=", tp%composition
        write(*,*) "tp%Tm=", tp%Tm
        write(*,*) "tp%Tb=", tp%Tb
        write(*,*) "tp%rhoA=", tp%rhoA
        write(*,*) "tp%rhoB=", tp%rhoB
        write(*,*) "tp%zmu1A=", tp%zmu1A
        write(*,*) "tp%zmu1B=", tp%zmu1B
        write(*,*) "tp%zmu2A=", tp%zmu2A
        write(*,*) "tp%zmu2B=", tp%zmu2B
        write(*,*) "tp%zmu2C=", tp%zmu2C
        write(*,*) "tp%zkA=", tp%zkA
        write(*,*) "tp%zkB=", tp%zkB
        write(*,*) "tp%cpA=", tp%cpA
        write(*,*) "tp%cpB=", tp%cpB
        write(*,*) "tp%cpC=", tp%cpC
        write(*,*) "tp%cpD=", tp%cpD
    endif
    close(18, status='keep')
end subroutine THMSGT

subroutine THMSST(salt, compo, tboil, impx)
    ! return the boiling temperature for the molten salts (If it is 0 in the MSTPDB is set to 5000 K)
    use t_saltdata
    character*12 :: salt, compo       ! salt formula and composition
    real, intent(out) :: tboil
    integer, intent(in) :: impx
    !
    type(tpdata) :: tp                ! Tupla with Thermophysical properties
    ! get the tpdata object for the specific salt
    call THMSGT(salt, compo, tp, impx) 
    if (tp%tb.eq.0.0) then
        tboil=5000
    else
        tboil=tp%tb
    endif
end subroutine THMSST

subroutine THMSPT(salt, compo, t, zrho, h, zk, zmu, zcp, impx)
    ! return the remaining thermohydraulics parameters as a function of the temperature (K) for molten salts
    use t_saltdata
    character*12 :: salt, compo       ! salt formula and composition
    real, intent(in) :: t
    real, intent(out) :: zrho, h, zk, zmu, zcp
    type(tpdata) :: tp                ! Tupla with Thermophysical properties
    integer, intent(in) :: impx
    !
    if(impx > 2) write(6,*) 'THSMPT: Molten salt thermophysical properties from MSTPDB'
    ! get the tpdata object for the specific salt
    call THMSGT(salt, compo, tp, impx) 
    zrho = rho(tp,t)
    zk = k(tp,t)
    zmu = mu(tp,t)
    zcp = cp(tp,t)
    h = zcp*t
    if (impx > 2) then
        write(*,*) 'THMSPT: ', salt, compo
        write(*,*) 'WEIGHT =', tp%weight
        write(*,*) 'TEMPERATURE =', t, '(K)'
        write(*,*) 'DENSITY =', zrho, '(kg/m3)'
        write(*,*) 'VISCOSITY =', zmu, '(kg/m2/s)'
        write(*,*) 'THERMAL CONDUCTIVITY =', zk, '(W/m/K)' 
        write(*,*) 'THERMAL CAPACITY =', zcp, '(J/K/kg)'
        write(*,*) 'SPECIFIC ENTHALPY =', h, '(J/kg)'
    endif
end subroutine THMSPT

subroutine THMSH(salt, compo, h, zrho, t, impx)
    ! return density and temperature given the entalphy
    use t_saltdata

    character*12 :: salt, compo             ! salt formula and composition
    real, intent(in) :: h
    real, intent(out) :: zrho, t
    type(tpdata) :: tp                      ! Tupla with Thermophysical properties

    integer, parameter :: rk=kind(0d0)
    integer, parameter :: degree=4
    real(rk) :: poly(degree+1), c1, c2, c3, c4, c5
    complex(rk) :: roots(degree)
    logical :: lfail
    integer, intent(in) :: impx
    !
    if (impx > 3) write(*,*) 'THMSH: Input entalpy h= ',h
    ! get the tpdata object for the specific salt
    call THMSGT(salt, compo, tp, impx) 

    ! solve polynomial h=CpT => 0 = D*T**4 + C*T**3 + B*T**2 + A*T - h
    a = tp%cpA/tp%weight*1000.0 
    b = tp%cpB/tp%weight*1000.0
    c = tp%cpC/tp%weight*1000.0
    d = tp%cpD/tp%weight*1000.0
    c1 = real(d, 8)
    c2 = real(b, 8)
    c3 = real(a, 8)
    c4 = real(-h, 8)
    c5 = real(c, 8)
    poly = [c5, c4, c3, c2, c1]
    npoly=degree
    do i=degree+1,1,-1
      if (poly(i) /= 0.) exit
      npoly=npoly-1
    enddo
    if (impx > 3) write(*,*) 'THMSH: Equation ', c1, '*T**4+', c2, '*T**3+', c3, '*T**2+', c4, '*T+', c5
    ! Note: In cmplx_roots_gen the polynomial is of the form poly(1) x^0 + poly(2) x^1 + poly(3) x^2 + ...
    call ALROOT(poly,npoly,roots,lfail)
    if(lfail) call XABORT('THMSH: foot finding failure.')
    if (impx > 3) write(*,*) 'THMSH: roots= ',roots
    do i = 1, degree
        if (aimag(roots(i)).eq.0.and.real(roots(i)).gt.0) then
            t = real(roots(i),4)
            exit
        endif
    end do
    zrho = rho(tp,t)
    if (impx > 3) write(*,*) 'THMSH: t = ', t, 'zrho = ',zrho
end subroutine THMSH
