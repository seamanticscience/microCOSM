! -*- f90 -*-
MODULE MOD_BOXES
    IMPLICIT NONE
#if defined(FOURBOX)
    INTEGER, PARAMETER :: nbox = 4
#elif defined(SIXBOX)
    INTEGER, PARAMETER :: nbox = 6
#else
! default to three box model
    INTEGER, PARAMETER :: nbox = 3
#endif        
END MODULE MOD_BOXES

MODULE MOD_DIMENSIONS
USE MOD_PRECISION
USE MOD_BOXES

IMPLICIT NONE

! geometry
!REAL(KIND=wp), DIMENSION(nbox)      :: dx, dy, dz, lat
REAL(KIND=wp), DIMENSION(nbox)      :: lat
REAL(KIND=wp), DIMENSION(nbox)      :: area, vol, invol, pressure
!REAL(KIND=wp), DIMENSION(nbox,nbox) :: K, R

CONTAINS
! Routines below are specifically put here because they require hard coding 
!   using information specifically for each set up:
!
! ESTABLISH_DIMENSIONS makes dimensional info (dx, dy, dz, areas, vols, etc) 
!                      available to the main model.
! CALC_ATMOS_MOLES uses ocean area to help scale number of moles in the atmosphere
!                      for use in pCO2 calculations 
! TRANSPORT calculates advective and diffusive fluxes between boxes.
! CALC_PSTAR calculates efficiency of the biological pump associated with 
!                      fluxes going into and out of the Southern Ocean box.


!=======================================================================
SUBROUTINE ESTABLISH_DIMENSIONS(dx,dy,dz,lat,depth,area,vol,invol,     &
                                              pressure)
USE MOD_BOXES
IMPLICIT NONE
REAL(KIND=wp), DIMENSION(nbox),      intent(in)  :: dx, dy, dz, lat,   &
                                                  depth

REAL(KIND=wp), DIMENSION(nbox),      intent(out) ::  area, vol, invol, &
                                                  pressure
!REAL(KIND=wp), DIMENSION(nbox,nbox), intent(out) :: K, R
REAL(KIND=wp)                                    :: m2deg

! applied pressure in bars for carbon system coefficients
pressure = (depth/10._wp) - 1._wp
                                                    
area     = dx * dy 
vol      = area * dz 
invol    = 1._wp / vol  

RETURN
END SUBROUTINE ESTABLISH_DIMENSIONS
!=======================================================================

!!=======================================================================
FUNCTION CALC_ATMOS_MOLES(area)
! Calculate the total number of moles in the atmosphere for pCO2
USE MOD_BOXES
IMPLICIT NONE
REAL(KIND=wp) :: CALC_ATMOS_MOLES

REAL(KIND=wp), DIMENSION(nbox), intent(in) :: area

! Mass dry atmosphere = (5.1352+/-0.0003)d18 kg (Trenberth & Smith,
! Journal of Climate 2005)
! and Mean molecular mass air = 28.97 g/mol (NASA earth fact sheet)
!  but need to scale by the ratio of Earth surface area to model surface area
! Earths area = 5.10082000d8 km2 * 1.e6 m2/km2 (NOAA earth fact sheet)
!       atmos_moles = atmos_moles * area(3)/(5.10082e8 * 1.e6)
#if defined(SIXBOX)
    CALC_ATMOS_MOLES = ((5.1352e18_wp * 1000._wp) / 28.97_wp)          &
                        *((area(2)+area(4)+area(6))/                   &
!                        *((area(1)+area(3)+area(5))/                   &
                          (5.10082e8_wp * 1.e6_wp))
#elif defined(FOURBOX)
    CALC_ATMOS_MOLES = ((5.1352e18_wp * 1000._wp) / 28.97_wp)          &
                        * (area(4)/(5.10082e8_wp * 1.e6_wp))
#else
    CALC_ATMOS_MOLES = ((5.1352e18_wp * 1000._wp) / 28.97_wp)          &
                        * (area(3)/(5.10082e8_wp * 1.e6_wp))
#endif

RETURN
END FUNCTION CALC_ATMOS_MOLES
!!=======================================================================

!!=======================================================================
!FUNCTION TRANSPORT(x, kappa, psi, invol)
!!atmosphere-3-box-ocean carbon cycle model
!!evaluate rates of change due to transport
!!mick follows, march 2015/ june 2016
! USE MOD_BOXES
! IMPLICIT NONE
! REAL(KIND=wp), DIMENSION(nbox)                  :: TRANSPORT
! REAL(KIND=wp), intent(in), DIMENSION(nbox)      :: x
! REAL(KIND=wp), intent(in), DIMENSION(nbox,nbox) :: kappa
! REAL(KIND=wp), intent(in)                       :: psi
! REAL(KIND=wp), intent(in), DIMENSION(nbox)      :: invol
! !
! TRANSPORT(1) = invol(1) * (                                            &
!             psi*(x(3)-x(1))                                            &
!           + kappa(3,1)*(x(3)-x(1))                                     &
!           + kappa(2,1)*(x(2)-x(1))                                     &
!             )
! TRANSPORT(2) = invol(2) * (                                            &
!             psi*(x(1)-x(2))                                            &
!           + kappa(1,2)*(x(1)-x(2))                                     &
!           + kappa(3,2)*(x(3)-x(2))                                     &
!             )
! TRANSPORT(3) = invol(3) * (                                            &
!             psi*(x(2)-x(3))                                            &
!           + kappa(2,3)*(x(2)-x(3))                                     &
!           + kappa(1,3)*(x(1)-x(3))                                     &
!             )
! 
!        RETURN
!        END FUNCTION TRANSPORT
!!=======================================================================

!=======================================================================
FUNCTION CALC_PSTAR(nutrient)
! atmosphere-3-box-ocean carbon cycle model
! evaluate rates of change due to transport
! mick follows, march 2015/ june 2016
USE MOD_BOXES
IMPLICIT NONE
REAL(KIND=wp) :: CALC_PSTAR

REAL(KIND=wp), DIMENSION(nbox), intent(in) :: nutrient

#if defined(SIXBOX)
    CALC_PSTAR = (nutrient(4)  - nutrient(3)) / nutrient(4) 
#elif defined(FOURBOX)
    CALC_PSTAR = (nutrient(4) - nutrient(1)) / nutrient(4) 
#else
    CALC_PSTAR = (nutrient(3)  - nutrient(1)) / nutrient(3) 
#endif

RETURN
END FUNCTION CALC_PSTAR
!=======================================================================

END MODULE MOD_DIMENSIONS