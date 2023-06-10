! -*- f90 -*-
MODULE MOD_BOXES
    IMPLICIT NONE
#if defined(FOURBOX)
    INTEGER, PARAMETER :: nbox = 4
#elif defined(FOURTEENBOX)
    INTEGER, PARAMETER :: nbox = 14
#else
! default to three box model
    INTEGER, PARAMETER :: nbox = 3
#endif        
END MODULE MOD_BOXES

MODULE MOD_DIMENSIONS
#if defined(USEDUALNUMAD)
    USE DNADMOD
#endif

USE MOD_PRECISION
USE MOD_BOXES

IMPLICIT NONE

! geometry
!REAL(KIND=wp), DIMENSION(nbox)      :: dx, dy, dz, latitude
!REAL(KIND=wp), DIMENSION(nbox)      :: latitude, depth
REAL(KIND=wp), DIMENSION(nbox)      :: area, vol, invol, pressure
!REAL(KIND=wp), DIMENSION(nbox,nbox) :: K, R

CONTAINS
! Routines below are specifically put here because they require hard coding 
!   using information specifically for each set up:
!
! ESTABLISH_DIMENSIONS makes dimensional info (dx, dy, dz, areas, vols, etc) 
!                      available to the main model.
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

! applied pressure in bars for carbon system coefficients
pressure = (depth/10._wp) - 1._wp

area     = dx * dy 
vol      = area * dz 
invol    = 1._wp / vol  

RETURN
END SUBROUTINE ESTABLISH_DIMENSIONS
!=======================================================================

!!=======================================================================
!FUNCTION TRANSPORT(x, kappa, psi, invol)
!!atmosphere-3-box-ocean carbon cycle model
!!evaluate rates of change due to transport
!!mick follows, march 2015/ june 2016
! USE MOD_BOXES
! IMPLICIT NONE
! REAL(KIND=wp), DIMENSION(nbox)                  :: TRANSPORT
! REAL(KIND=wp), intent(in), DIMENSION(nbox)      :: x
! REAL(KIND=wp), intent(in), DIMENSION(nbox,nbox) :: kmask, pmask
!! REAL(KIND=wp), intent(in), DIMENSION(nbox,nbox) :: kappa
! REAL(KIND=wp), intent(in)                       :: psi, kappa
! REAL(KIND=wp), intent(in), DIMENSION(nbox)      :: invol
!
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

#if defined(FOURBOX)
    CALC_PSTAR = (nutrient(4)  - nutrient(1)) / nutrient(4) 
#elif defined(FOURTEENBOX)
    CALC_PSTAR = (nutrient(14) - nutrient(2)) / nutrient(14) 
#else
    CALC_PSTAR = (nutrient(3)  - nutrient(1)) / nutrient(3) 
#endif

RETURN
END FUNCTION CALC_PSTAR
!=======================================================================

END MODULE MOD_DIMENSIONS