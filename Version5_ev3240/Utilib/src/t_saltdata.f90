!
!-----------------------------------------------------------------------
!
!Purpose:
! module t_saltdata definition for the Fortran-2003 wrapper used to
! calculate Molten Salt thermophysical properties using data from the
! MSTPDB-TP Database
!
!Copyright:
! Copyright (C) 2023 Cristian Garrido Tamm
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
!Author(s): Cristian Garrido Tamm
!
!-----------------------------------------------------------------------
!
module t_saltdata
    implicit none
    type :: tpdata                         ! Salt thermophysical data structure
        character(len=32)  :: formula      ! Chemical formula of the salt (e.g 'LiF-BeF2')
        real  :: weight                    ! Molecular weight 
        character(len=32)  :: composition  ! Salt composition (e.g 'Pure Salt' or '0.66-0.34')
        real  :: Tm                        ! Melting Temperature (K)
        real  :: Tb                        ! Boiling Temperature (K)
        real  :: rhoA                      ! A coefficient in density (g/cm3):  A - BT(K)
        real  :: rhoB                      ! A coefficient in density (g/cm3):  A - BT(K)
        real  :: zmu1A                     ! A coeff in Viscosity (mN*s/m2):  A*exp(B/(R*T(K)))
        real  :: zmu1B                     ! B coeff in Viscosity (mN*s/m2):  A*exp(B/(R*T(K)))
        real  :: zmu2A                     ! A coeff in Viscosity (mN*s/m2):  10^(A + B/T + C/T**2)
        real  :: zmu2B                     ! B coeff in Viscosity (mN*s/m2):  10^(A + B/T + C/T**2)
        real  :: zmu2C                     ! C coeff in Viscosity (mN*s/m2):  10^(A + B/T + C/T**2)
        real  :: zkA                       ! A coefficient in Thermal Conductivity (W/m K):  A - B*T(K)
        real  :: zkB                       ! B coefficient in Thermal Conductivity (W/m K):  A - B*T(K)
        real :: cpA                        ! A coefficient in Heat Capacity (J/K mol) : A + B*T(K) + C*T-2(K) + D*T2(K)
        real :: cpB                        ! B coefficient in Heat Capacity (J/K mol) : A + B*T(K) + C*T-2(K) + D*T2(K)
        real :: cpC                        ! C coefficient in Heat Capacity (J/K mol) : A + B*T(K) + C*T-2(K) + D*T2(K)
        real :: cpD                        ! D coefficient in Heat Capacity (J/K mol) : A + B*T(K) + C*T-2(K) + D*T2(K)
    end type tpdata 

contains
    real function rho(self, t) result(res)
        ! Computes density (g/cm3) as A - B*T(K)
        type(tpdata), intent(in) :: self
        real, intent(in) :: t              ! Temperature (K)
        real, parameter :: f = 1000.0      ! g/cm3 to kg/m3
        res = (self%rhoA - self%rhoB * t) * f   
    end function
    real function mu(self, t) result(res)
        ! Computes Viscosity (mN*s/m2) as A*exp(B/(R*T(K))) or 10^(A + B/T + C/T**2)
        type(tpdata), intent(in) :: self
        real, intent(in) :: t              ! Temperature (K)
        real, parameter :: R = 8.314       ! J/K mol
        real, parameter :: f = 1000.0      ! mN*s/m2 to kg/m2/s
        real :: zmu1, zmu2
        zmu1 = self%zmu1A*exp(self%zmu1B/(R*t))
        zmu2 = 10**(self%zmu2A + self%zmu2B/t + self%zmu2C/T**2)
        res = max(zmu1, zmu2)/1000
    end function
    real function k(self, t) result(res)
        ! Computes Thermal Conductivity (W/m K) as A - B*T(K)
       type(tpdata), intent(in) :: self
        real, intent(in) :: t              ! Temperature (K)
        res = self%zkA - self%zkB * t
    end function
    real function cp(self, t) result(res)
        ! Computes Heat Capacity (J/K mol) as A + B*T(K) + C*T**-2(K) + D*T**2(K)
        type(tpdata), intent(in) :: self
        real, intent(in) :: t              ! Temperature (K)
        res = self%cpA + self%cpB * t + self%cpC * t**(-2) + self%cpD * t**2
        res = res / self%weight * 1000     ! J/K/kg
    end function
end module t_saltdata
