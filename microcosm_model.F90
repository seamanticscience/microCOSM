! atmospherd-3-box-ocean carbon cycle model
! mick follows, march 2015
! convert matlab to f90 - march/june 2016
! significant work by jonathan lauderale june 2016-oct 2019
! refresh, parallelization, and expansion by jonathan lauderdale 2020-2021
!
!         --------------------------------
!         |          ATMOSPHERE          |
!         |                              |
!         --------------------------------
!         | SUR- |                       |
!         | FACE ->   SURFACE  2         |
!         |  1   |                   |   |
!         |---^----------------------v---|
!         |   |                          | 
!         |          DEEP OCEAN  3       |
!         |                              | 
!         |                              | 
!         --------------------------------
!=======================================================================
       PROGRAM MICROCOSM_MODEL
!=======================================================================
       USE MOD_MODELMAIN
       USE MOD_PRECISION
       
       IMPLICIT NONE
       REAL(KIND=wp) :: t1in 
       REAL(KIND=wp) :: t2in 
       REAL(KIND=wp) :: t3in 
       REAL(KIND=wp) :: s1in 
       REAL(KIND=wp) :: s2in 
       REAL(KIND=wp) :: s3in 
       REAL(KIND=wp) :: c1in 
       REAL(KIND=wp) :: c2in 
       REAL(KIND=wp) :: c3in 
       REAL(KIND=wp) :: a1in 
       REAL(KIND=wp) :: a2in 
       REAL(KIND=wp) :: a3in 
       REAL(KIND=wp) :: p1in 
       REAL(KIND=wp) :: p2in 
       REAL(KIND=wp) :: p3in 
       REAL(KIND=wp) :: n1in 
       REAL(KIND=wp) :: n2in 
       REAL(KIND=wp) :: n3in 
       REAL(KIND=wp) :: f1in 
       REAL(KIND=wp) :: f2in 
       REAL(KIND=wp) :: f3in 
       REAL(KIND=wp) :: l1in 
       REAL(KIND=wp) :: l2in 
       REAL(KIND=wp) :: l3in 
       REAL(KIND=wp) :: atpco2in 

       INTEGER :: id
       REAL(KIND=wp) :: maxyears
       REAL(KIND=wp) :: outputyears
       REAL(KIND=wp) :: gaovla_opt
       REAL(KIND=wp) :: gamma_Fe
       REAL(KIND=wp) :: lt_lifetime
       REAL(KIND=wp) :: depfactor
       REAL(KIND=wp) :: ventfactor
       REAL(KIND=wp) :: alpha_yr
       REAL(KIND=wp) :: dlambdadz
       REAL(KIND=wp) :: psi
       
       ! Cannot use allocatable arrays with f2py :(
       INTEGER, PARAMETER  :: outstepmax = 1000
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: tout  
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: t1out 
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: t2out 
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: t3out 
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: s1out 
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: s2out 
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: s3out 
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: c1out 
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: c2out 
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: c3out 
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: a1out 
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: a2out 
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: a3out 
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: p1out 
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: p2out 
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: p3out 
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: n1out 
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: n2out 
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: n3out 
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: f1out 
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: f2out 
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: f3out 
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: l1out 
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: l2out 
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: l3out 
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: ep1out
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: ep2out
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: nlout 
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: psout        
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: o1pco2out
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: o2pco2out
       REAL(KIND=wp), DIMENSION(outstepmax+1) :: atpco2out

       ! Input some initial parameters
       maxyears   = 1.e4
       outputyears= 1.e1
       t1in   =    2.0
       t2in   =   20.0
       t3in   =    4.0
       s1in   =   34.00
       s2in   =   35.50
       s3in   =   34.75
       ! Initial concentrations in mol/kg
       ! Make sure to compile with -DFIXATMPCO2 first
       !c1in   = 2100.e-6
       !c2in   = 2100.e-6
       !c3in   = 2350.e-6
       !a1in   = 2300.e-6
       !a2in   = 2300.e-6
       !a3in   = 2400.e-6
       !p1in   =    2.0e-6
       !p2in   =    0.0e-6
       !p3in   =    2.5e-6
       !n1in   =   32.0e-6
       !n2in   =    0.0e-6
       !n3in   =   36.0e-6
       !f1in   =    0.e-9
       !f2in   =    0.e-9
       !f3in   =    0.e-9
       !l1in   =    0.e-9
       !l2in   =    0.e-9
       !l3in   =    0.e-9
       !atpco2in =  280.

       ! Initial concentrations in mol/kg
       ! Here are some equilibrated values run for 100,000 yrs (round-off error notwithstanding)
       ! Make sure to compile without -DFIXATMPCO2
       c1in   = 2264.67564e-6 
       c2in   = 2103.48757e-6
       c3in   = 2364.66971e-6
       a1in   = 2395.54471e-6
       a2in   = 2387.42965e-6
       a3in   = 2399.11941e-6
       p1in   =    1.81089e-6
       p2in   =    0.25031e-6
       p3in   =    2.49834e-6
       n1in   =   25.01353e-6
       n2in   =    0.04412e-6
       n3in   =   36.01262e-6
       f1in   =    0.00377e-9
       f2in   =    0.49776e-9
       f3in   =    0.58847e-9
       l1in   =    2.08548e-9
       l2in   =    1.56387e-9
       l3in   =    1.62029e-9
       atpco2in =  280.00000

       ! Input some example parameters
       ! Gamma over lambda for ligands "optimum" value (Lauderdale et al 2020)
       gaovla_opt = 4398.
       ! Gamma ligand production rate (in phosphate, not carbon, units)
       gamma_Fe    = 5.e-5*106.
       ! Lambda ligand lifetime (s)
       lt_lifetime =  1/((gamma_Fe/106.)/gaovla_opt)
       ! Dust deposition in g Fe m-2 year-1
       depfactor   = 0.15
       ! Hydrothermal vent input of 1 Gmol/yr (Tagliabue et al., 2010)
       ! mol Fe/yr * g/mol * 1/area  == g Fe m-2 year-1....
       !divide by 2.5e-3 because fe_sol is multiplied again within model.
       ventfactor  = (1.e9*56.)/(17.e6*16.e6*2.5e-3)
       ! Biological production maximum rate (mol P/yr)
       alpha_yr    = 6.e-6
       ! Deep ocean box lifetime modifier to capture the gradient due to
       ! photodegradation near the surface and slower loss in the deep
       dlambdadz   = 1.e-2
       ! Overturning rate (m3/s)
       psi         = 20.0e6
       ! File number identifier
       id          = 1
            
       call model(maxyears, outputyears,                               &
            t1in,t2in,t3in,                                            &
            s1in,s2in,s3in,                                            &
            c1in,c2in,c3in,                                            &
            a1in,a2in,a3in,                                            &
            p1in,p2in,p3in,                                            &
            n1in,n2in,n3in,                                            &
            f1in,f2in,f3in,                                            &
            l1in,l2in,l3in,                                            &
            atpco2in,                                                  &
            tout,                                                      &
            t1out,t2out,t3out,                                         &
            s1out,s2out,s3out,                                         &
            c1out,c2out,c3out,                                         &
            a1out,a2out,a3out,                                         &
            p1out,p2out,p3out,                                         &
            n1out,n2out,n3out,                                         &
            f1out,f2out,f3out,                                         &
            l1out,l2out,l3out,                                         &
            ep1out,ep2out,nlout,psout,                                 &
            o1pco2out,o2pco2out,atpco2out,                             &
            gamma_Fe,lt_lifetime,depfactor,ventfactor,                 &
            alpha_yr,dlambdadz,psi,id)

!=======================================================================
       END PROGRAM MICROCOSM_MODEL
!=======================================================================
