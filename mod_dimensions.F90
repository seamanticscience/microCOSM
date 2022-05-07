! -*- f90 -*-
MODULE MOD_BOXES
    IMPLICIT NONE
    INTEGER, PARAMETER :: nbox = 3
END MODULE MOD_BOXES

MODULE MOD_DIMENSIONS
#if defined(USEDUALNUMAD)
    USE DNADMOD
#endif

USE MOD_PRECISION
USE MOD_BOXES

IMPLICIT NONE

! geometry
!REAL(KIND=wp), DIMENSION(nbox)      :: dx, dy, dz, lat
REAL(KIND=wp), DIMENSION(nbox)      :: lat
REAL(KIND=wp), DIMENSION(nbox)      :: area, vol, invol, depth, pressure
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
SUBROUTINE ESTABLISH_DIMENSIONS(dx,dy,dz,lat,area,vol,invol,           &
                                              depth,pressure)
USE MOD_BOXES
IMPLICIT NONE
REAL(KIND=wp), DIMENSION(nbox),      intent(in)  :: dx, dy, dz

REAL(KIND=wp), DIMENSION(nbox),      intent(out) ::                    &
                                                 lat, area, vol, invol,&
                                                 depth, pressure
!REAL(KIND=wp), DIMENSION(nbox,nbox), intent(out) :: K, R
REAL(KIND=wp)                                    :: m2deg

!dx   = [ 17.0e6_wp, 17.0e6_wp, 17.0e6_wp ]
!dy   = [  4.0e6_wp, 12.0e6_wp, 16.0e6_wp ]  
!dz   = [ 50.0_wp,   50.0_wp, 5050.0_wp   ]

! depth in m or decibars
depth = [ 25.0_wp,   25.0_wp, 2575.0_wp ]
! applied pressure in bars for carbon system coefficients
pressure = (depth/10._wp) - 1._wp

m2deg = 180._wp/(dy(1)+dy(2))  
lat = [  -90._wp+(dy(1)       /2._wp) *m2deg,                          &        
         -90._wp+(dy(1)+(dy(2)/2._wp))*m2deg,                          &
         -90._wp+(dy(3)       /2._wp) *m2deg                           &
       ]
                                                    
area     = dx * dy 
vol      = area * dz 
invol    = 1._wp / vol  

RETURN
END SUBROUTINE ESTABLISH_DIMENSIONS
!=======================================================================

!=======================================================================
FUNCTION TRANSPORT(x, kappa, psi, invol)
! atmosphere-3-box-ocean carbon cycle model
! evaluate rates of change due to transport
! mick follows, march 2015/ june 2016
USE MOD_BOXES
IMPLICIT NONE
REAL(KIND=wp), DIMENSION(nbox)                  :: TRANSPORT
REAL(KIND=wp), intent(in), DIMENSION(nbox)      :: x
REAL(KIND=wp), intent(in), DIMENSION(nbox,nbox) :: kappa
REAL(KIND=wp), intent(in)                       :: psi
REAL(KIND=wp), intent(in), DIMENSION(nbox)      :: invol
!
TRANSPORT(1) = invol(1) * (                                            &
            psi*(x(3)-x(1))                                            &
          + kappa(3,1)*(x(3)-x(1))                                     &
          + kappa(2,1)*(x(2)-x(1))                                     &
            )
TRANSPORT(2) = invol(2) * (                                            &
            psi*(x(1)-x(2))                                            &
          + kappa(1,2)*(x(1)-x(2))                                     &
          + kappa(3,2)*(x(3)-x(2))                                     &
            )
TRANSPORT(3) = invol(3) * (                                            &
            psi*(x(2)-x(3))                                            &
          + kappa(2,3)*(x(2)-x(3))                                     &
          + kappa(1,3)*(x(1)-x(3))                                     &
            )

       RETURN
       END FUNCTION TRANSPORT
!=======================================================================

! !=======================================================================
! FUNCTION TRANSPORT(x, kappa, psi, invol)
! !atmosphere-5-box-ocean carbon cycle model
! ! a la Follows, Ito and Marotzke, 2002 GBC 
! !evaluate rates of change due to transport
! 
! USE MOD_BOXES
! IMPLICIT NONE
! REAL(KIND=wp), DIMENSION(nbox)                  :: TRANSPORT
! REAL(KIND=wp), intent(in), DIMENSION(nbox)      :: x
! REAL(KIND=wp), intent(in), DIMENSION(nbox,nbox) :: kappa
! REAL(KIND=wp), intent(in)                       :: psi
! REAL(KIND=wp), intent(in), DIMENSION(nbox)      :: invol
! !
! TRANSPORT(1) = invol(1) * (                                            &
!             psi(1)*(x(5)-x(1))                                         &
!           + kappa(5,1)*(x(5)-x(1))                                     &
!           + kappa(2,1)*(x(2)-x(1))                                     &
!           + kappa(3,1)*(x(3)-x(1)))                                    &
!             )
! TRANSPORT(2) = invol(2) * (                                            &
!             psi(1)*(x(1)-x(2))                                         &
!           + psi(2)*(x(4)-x(2))  
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
! !=======================================================================

!=======================================================================
FUNCTION CALC_PSTAR(nutrient)
! atmosphere-3-box-ocean carbon cycle model
! evaluate rates of change due to transport
! mick follows, march 2015/ june 2016
USE MOD_BOXES
IMPLICIT NONE
REAL(KIND=wp) :: CALC_PSTAR

REAL(KIND=wp), DIMENSION(nbox), intent(in) :: nutrient

CALC_PSTAR = (nutrient(3) - nutrient(1)) / nutrient(3) 

RETURN
END FUNCTION CALC_PSTAR
!=======================================================================

END MODULE MOD_DIMENSIONS