MODULE MOD_COMMON
!variable declarations etc for microCOSM model

USE MOD_PRECISION
USE MOD_BOXES
IMPLICIT NONE

! timestepping variables
REAL(KIND=wp), PARAMETER :: dt         = 86400._wp 
REAL(KIND=wp), PARAMETER :: s_per_d    = 86400._wp
REAL(KIND=wp), PARAMETER :: d_per_yr   = 365._wp
REAL(KIND=wp), PARAMETER :: s_per_yr   = 31536000._wp
INTEGER                  :: nstepmax

! biogeochemical tracers internal
REAL(KIND=wp), DIMENSION(nbox) :: theta   
REAL(KIND=wp), DIMENSION(nbox) :: salt   
REAL(KIND=wp), DIMENSION(nbox) :: dic   
REAL(KIND=wp), DIMENSION(nbox) :: alk   
REAL(KIND=wp), DIMENSION(nbox) :: po4   
REAL(KIND=wp), DIMENSION(nbox) :: no3   
REAL(KIND=wp), DIMENSION(nbox) :: fet   
REAL(KIND=wp), DIMENSION(nbox) :: lt    

! extra biogeochem...
REAL(KIND=wp), DIMENSION(nbox) :: ph
REAL(KIND=wp), DIMENSION(nbox) :: sit
REAL(KIND=wp), DIMENSION(nbox) :: carb
REAL(KIND=wp), DIMENSION(nbox) :: feprime
REAL(KIND=wp), DIMENSION(nbox) :: bioP
REAL(KIND=wp)                  :: pstar 

! some arrays for conversions...
REAL(KIND=wp), DIMENSION(nbox) :: dicM  
REAL(KIND=wp), DIMENSION(nbox) :: alkM  
REAL(KIND=wp), DIMENSION(nbox) :: po4M  
REAL(KIND=wp), DIMENSION(nbox) :: no3M  
REAL(KIND=wp), DIMENSION(nbox) :: fetM  
REAL(KIND=wp), DIMENSION(nbox) :: ltM   
REAL(KIND=wp), DIMENSION(nbox) :: sitM   
REAL(KIND=wp), DIMENSION(nbox) :: pco2M
REAL(KIND=wp) :: pco2A
! time derivatives 
REAL(KIND=wp), DIMENSION(nbox) :: dthetadt  
REAL(KIND=wp), DIMENSION(nbox) :: dsaltdt  
REAL(KIND=wp), DIMENSION(nbox) :: ddicdt  
REAL(KIND=wp), DIMENSION(nbox) :: dalkdt  
REAL(KIND=wp), DIMENSION(nbox) :: dpo4dt  
REAL(KIND=wp), DIMENSION(nbox) :: dno3dt  
REAL(KIND=wp), DIMENSION(nbox) :: dfetdt  
REAL(KIND=wp), DIMENSION(nbox) :: dltdt  
REAL(KIND=wp), DIMENSION(nbox) :: dsitdt  

! Redfield ratios      
REAL(KIND=wp), PARAMETER       :: rCP  = 106._wp
REAL(KIND=wp), PARAMETER       :: rNP  = 16._wp
REAL(KIND=wp), PARAMETER       :: rPO2 = 170._wp     
REAL(KIND=wp), PARAMETER       :: rCN  = rCP/rNP 
REAL(KIND=wp), PARAMETER       :: rCO2 = rCP/rPO2 
REAL(KIND=wp), PARAMETER       :: rCFe = rCP/1.e-3_wp 
REAL(KIND=wp), PARAMETER       :: rSIP = 15._wp
REAL(KIND=wp), PARAMETER       :: rCACO3 = 10.e-2_wp 

! conversion factor, average density (kg m-3)
REAL(KIND=wp), PARAMETER       :: conv = 1024.5_wp 
REAL(KIND=wp), PARAMETER       :: umolkg2molm3 = conv * 1.e-6_wp 
REAL(KIND=wp), PARAMETER       :: nmolkg2molm3 = conv * 1.e-9_wp 
REAL(KIND=wp), PARAMETER       :: uatm2atm     = 1.e-6_wp 
REAL(KIND=wp), PARAMETER       :: molps2gtcyr  =                       &
                                    rCP * 12._wp * s_per_yr * 1.e-15_wp 

REAL(KIND=wp), DIMENSION(nbox) :: pco2ocean  
REAL(KIND=wp), DIMENSION(nbox) :: fluxCO2  
REAL(KIND=wp)                  :: netco2flux
REAL(KIND=wp)                  :: atmos_moles
REAL(KIND=wp)                  :: atmos_carbon
REAL(KIND=wp)                  :: pco2atmos
! Iron cycle parameters 
! atomic weight of iron = 56
REAL(KIND=wp), PARAMETER       :: weight_fe = 56._wp  
!solubility of iron:
REAL(KIND=wp), PARAMETER       :: fe_sol = 0.0025_wp 
! conditional stability FeL: (mol kg-1)-1
REAL(KIND=wp), PARAMETER       :: beta   = 1.0e9_wp 
! Free Fe scavenging rate: (s-1) 
REAL(KIND=wp), PARAMETER       :: Kscav  = 1.0e-7_wp    
! relaxfe (s) 
REAL(KIND=wp), PARAMETER       :: relaxfe = 0.01_wp * s_per_yr  
! multiplier to test sensitivity to dust deposition
!REAL(KIND=wp), intent(in)      :: depfactor  
!REAL(KIND=wp), intent(in)      :: ventfactor
! iron input rate
REAL(KIND=wp), DIMENSION(nbox) :: fe_depo
! Dynamic Ligand variables
REAL(KIND=wp), DIMENSION(nbox) :: lambda 

! export related
! half saturation for iron limitation (mol m-3)
REAL(KIND=wp), PARAMETER       :: kfe    =   0.1e-9_wp*conv
! half saturation for phosphate limitation (mol m-3)
REAL(KIND=wp), PARAMETER       :: kpo4   =  0.1e-6_wp*conv
REAL(KIND=wp), PARAMETER       :: kno3   =  0.1e-6_wp*conv*rNP
! half saturation for light limitation (W m-2)
REAL(KIND=wp), PARAMETER       :: klight = 30._wp
REAL(KIND=wp), DIMENSION(nbox) :: export
! max export prodution rate: phosphorus units! (mol P m-3 s-1)
!REAL(KIND=wp), intent(in)      :: alpha_yr 
REAL(KIND=wp)                  :: alpha 
REAL(KIND=wp)                  :: itim
REAL(KIND=wp), DIMENSION(nbox) ::  light  
REAL(KIND=wp), DIMENSION(nbox) :: ilimit
REAL(KIND=wp), DIMENSION(nbox) :: plimit
REAL(KIND=wp), DIMENSION(nbox) :: nlimit
REAL(KIND=wp), DIMENSION(nbox) :: flimit
! nutrient limitation codes
INTEGER                        :: lim

END MODULE MOD_COMMON
