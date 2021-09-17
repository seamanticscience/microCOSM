!variable declarations etc for microCOSM

! Declare all parameters first
! timestepping variables
        REAL(KIND=wp), PARAMETER :: dt = 86400.0 
        REAL(KIND=wp), PARAMETER :: s_per_d = 86400.0
        REAL(KIND=wp), PARAMETER :: d_per_yr = 365.0
        REAL(KIND=wp), PARAMETER :: s_per_yr = 31536000.0
!        REAL(KIND=wp), PARAMETER :: 
!       &         s_per_yr = s_per_d*d_per_yr
        INTEGER, PARAMETER  :: outstepmax = 1000
        REAL(KIND=wp), intent(in)  :: maxyears, outputyears
        INTEGER :: nstepmax
        
! Which iteration of model parameters is this
        INTEGER, intent(in) :: id
        
! geometry
       INTEGER, PARAMETER :: nbox = 3 

       REAL(KIND=wp) :: m2deg
       REAL(KIND=wp), DIMENSION(nbox) :: dx
       REAL(KIND=wp), DIMENSION(nbox) :: dy
       REAL(KIND=wp), DIMENSION(nbox) :: dz
       REAL(KIND=wp), DIMENSION(nbox) :: area 
       REAL(KIND=wp), DIMENSION(nbox) :: vol
       REAL(KIND=wp), DIMENSION(nbox) :: invol
       REAL(KIND=wp), DIMENSION(nbox) :: lat
       
! circulation and physics
! overturning circulation (m3 s-1) provided as an input
       REAL(KIND=wp), intent(in) :: psi
! define array of mixing rates
       REAL(KIND=wp), DIMENSION(nbox,nbox) :: K

! biogeochemical tracers input   
       REAL(KIND=wp), intent(in) :: t1in 
       REAL(KIND=wp), intent(in) :: t2in 
       REAL(KIND=wp), intent(in) :: t3in 
       REAL(KIND=wp), intent(in) :: s1in 
       REAL(KIND=wp), intent(in) :: s2in 
       REAL(KIND=wp), intent(in) :: s3in 
       REAL(KIND=wp), intent(in) :: c1in 
       REAL(KIND=wp), intent(in) :: c2in 
       REAL(KIND=wp), intent(in) :: c3in 
       REAL(KIND=wp), intent(in) :: a1in 
       REAL(KIND=wp), intent(in) :: a2in 
       REAL(KIND=wp), intent(in) :: a3in 
       REAL(KIND=wp), intent(in) :: p1in 
       REAL(KIND=wp), intent(in) :: p2in 
       REAL(KIND=wp), intent(in) :: p3in 
       REAL(KIND=wp), intent(in) :: n1in 
       REAL(KIND=wp), intent(in) :: n2in 
       REAL(KIND=wp), intent(in) :: n3in 
       REAL(KIND=wp), intent(in) :: f1in 
       REAL(KIND=wp), intent(in) :: f2in 
       REAL(KIND=wp), intent(in) :: f3in 
       REAL(KIND=wp), intent(in) :: l1in 
       REAL(KIND=wp), intent(in) :: l2in 
       REAL(KIND=wp), intent(in) :: l3in 
       REAL(KIND=wp), intent(in) :: atpco2in 

! biogeochemical tracers internal
       REAL(KIND=wp), DIMENSION(nbox) :: theta   
       REAL(KIND=wp), DIMENSION(nbox) :: salt   
       REAL(KIND=wp), DIMENSION(nbox) :: dic   
       REAL(KIND=wp), DIMENSION(nbox) :: alk   
       REAL(KIND=wp), DIMENSION(nbox) :: po4   
       REAL(KIND=wp), DIMENSION(nbox) :: no3   
       REAL(KIND=wp), DIMENSION(nbox) :: fet   
       REAL(KIND=wp), DIMENSION(nbox) :: lt    

! biogeochemical tracers output 
! Cannot use allocatable arrays with f2py :(
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: tout  
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: t1out 
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: t2out 
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: t3out 
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: s1out 
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: s2out 
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: s3out 
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: c1out 
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: c2out 
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: c3out 
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: a1out 
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: a2out 
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: a3out 
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: p1out 
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: p2out 
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: p3out 
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: n1out 
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: n2out 
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: n3out 
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: f1out 
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: f2out 
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: f3out 
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: l1out 
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: l2out 
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: l3out 
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: ep1out
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: ep2out
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: nlout 
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: psout 
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: o1pco2out
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: o2pco2out  
       REAL(KIND=wp), intent(out), DIMENSION(outstepmax+1) :: atpco2out

! extra biogeochem...
       REAL(KIND=wp), DIMENSION(nbox) :: ph
       REAL(KIND=wp), DIMENSION(nbox) :: sit
       REAL(KIND=wp), DIMENSION(nbox) :: carb
       REAL(KIND=wp), DIMENSION(nbox) :: feprime
       REAL(KIND=wp) :: remin
       REAL(KIND=wp) :: pstar 

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
       REAL(KIND=wp), PARAMETER :: rCP  = 106.0
       REAL(KIND=wp), PARAMETER :: rNP  = 16.0
       REAL(KIND=wp), PARAMETER :: rPO2 = 170.0     
       REAL(KIND=wp), PARAMETER :: rCN  = rCP/rNP 
       REAL(KIND=wp), PARAMETER :: rCO2 = rCP/rPO2 
       REAL(KIND=wp), PARAMETER :: rCFe = rCP/1.e-3
       REAL(KIND=wp), PARAMETER :: rSIP = 15.0
       REAL(KIND=wp), PARAMETER :: rCACO3 = 10.e-2

! conversion factor, average density (kg m-3)
       REAL(KIND=wp), PARAMETER :: conv = 1024.5
! gas transfer piston velocity       
       REAL(KIND=wp) :: Kg
       REAL(KIND=wp), DIMENSION(nbox) :: wind  
       REAL(KIND=wp), DIMENSION(nbox) :: fice  
       REAL(KIND=wp), DIMENSION(nbox) :: fluxCO2  
       REAL(KIND=wp), DIMENSION(nbox) :: pco2oce  
       REAL(KIND=wp) :: atmos_moles
       REAL(KIND=wp) :: atmos_carbon
       REAL(KIND=wp) :: pco2atm
       
! Iron cycle parameters 
! atomic weight of iron = 56
       REAL(KIND=wp), PARAMETER  :: weight_fe = 56.0  
!solubility of iron:
       REAL(KIND=wp), PARAMETER  :: fe_sol = 0.0025
! conditional stability FeL: (mol kg-1)-1       
!       REAL(KIND=wp), PARAMETER  :: beta   = 1.0e9
       REAL(KIND=wp), PARAMETER  :: beta   = 1.0e9
! Free Fe scavenging rate: (s-1) 
       REAL(KIND=wp), PARAMETER  :: Kscav  = 1.0e-7      
! relaxfe (s) 
       REAL(KIND=wp), PARAMETER  :: relaxfe = 0.01*31536000.0 
! multiplier to test sensitivity to dust deposition
!       REAL(KIND=wp), PARAMETER  :: depfactor = 1.0  
       REAL(KIND=wp), intent(in) :: depfactor  
       REAL(KIND=wp), intent(in) :: ventfactor
! iron input rate
       REAL(KIND=wp), DIMENSION(nbox) :: fe_depo
! Dynamic Ligand variables
! gamma_Fe is fraction of "export"/"remin" as ligand. Must be < 1!
       REAL(KIND=wp), intent(in) :: gamma_Fe
!      lt_lifetime is the degradation rate of ligand (s) 
!          2.0*3.0e7 (ie 2 yrs default)
       REAL(KIND=wp), intent(in) :: lt_lifetime
! dlambdadz is the gradient in timescale with depth with a default
!       of 0.01 (ie 100 longer in the deep ocean)
       REAL(KIND=wp), intent(in) :: dlambdadz
       REAL(KIND=wp), DIMENSION(nbox) :: lambda 
       
! export related
! half saturation for iron limitation (mol m-3)
       REAL(KIND=wp), PARAMETER :: kfe    =   0.1e-9*conv
! half saturation for phosphate limitation (mol m-3)
       REAL(KIND=wp), PARAMETER :: kpo4   =  0.1e-6*conv
       REAL(KIND=wp), PARAMETER :: kno3   =  0.1e-6*conv*rNP
! half saturation for light limitation (W m-2)       
       REAL(KIND=wp), PARAMETER :: klight = 30.0
       REAL(KIND=wp), DIMENSION(nbox) :: export
! max export prodution rate: phosphorus units! (mol P m-3 s-1)
       REAL(KIND=wp), intent(in) ::  alpha_yr 
       REAL(KIND=wp) :: alpha 
       REAL(KIND=wp) :: itim
       REAL(KIND=wp) :: ilat
       REAL(KIND=wp) :: ilight  
       REAL(KIND=wp) :: ilimit
       REAL(KIND=wp) :: plimit
       REAL(KIND=wp) :: nlimit
       REAL(KIND=wp) :: flimit



       

       

       
