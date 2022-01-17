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
       USE MOD_PRECISION
       USE MOD_BOXES
       USE MOD_MODELMAIN

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

       INTEGER :: outstepmax, id
       REAL(KIND=wp) :: maxyears
       REAL(KIND=wp) :: outputyears
       REAL(KIND=wp) :: gaovla_opt
       REAL(KIND=wp) :: gamma_Fe
       REAL(KIND=wp) :: lt_lifetime
       REAL(KIND=wp) :: fe_input1
       REAL(KIND=wp) :: fe_input2
       REAL(KIND=wp) :: fe_input3
       REAL(KIND=wp) :: alpha_yr
       REAL(KIND=wp) :: dlambdadz1
       REAL(KIND=wp) :: dlambdadz2
       REAL(KIND=wp) :: dlambdadz3
       REAL(KIND=wp) :: psi
       REAL(KIND=wp) :: wind1
       REAL(KIND=wp) :: wind2
       REAL(KIND=wp) :: wind3
       REAL(KIND=wp) :: fopen1
       REAL(KIND=wp) :: fopen2
       REAL(KIND=wp) :: fopen3
       
       REAL(KIND=wp), dimension (:,:), allocatable ::                  &
            thout,                                                     &
            sout,                                                      &
            cout,                                                      &
            aout,                                                      &
            pout,                                                      &
            nout,                                                      &
            fout,                                                      &
            lout,                                                      &
            expout,                                                    &
            pco2out
            
       REAL(KIND=wp), dimension (:), allocatable   ::                  &
            tout,                                                      &
            psout,                                                     &
            atpco2out

       INTEGER, dimension (:), allocatable   ::                        &
            nlout
       
       ! Input some initial parameters
       maxyears   = 1.e5_wp
       outputyears= 1.e2_wp
       outstepmax = int((maxyears/outputyears)+1)
       
       ! allocate memory
       allocate ( tout      (outstepmax) )
       allocate ( nlout     (outstepmax) )
       allocate ( psout     (outstepmax) )
       allocate ( atpco2out (outstepmax) )
       allocate ( thout     (outstepmax,nbox) )
       allocate ( sout      (outstepmax,nbox) )
       allocate ( cout      (outstepmax,nbox) )
       allocate ( aout      (outstepmax,nbox) )
       allocate ( pout      (outstepmax,nbox) )
       allocate ( nout      (outstepmax,nbox) )
       allocate ( fout      (outstepmax,nbox) )
       allocate ( lout      (outstepmax,nbox) )
       allocate ( expout    (outstepmax,nbox) )
       allocate ( pco2out   (outstepmax,nbox) )
       
       ! Input arguements
       t1in   =    2._wp
       t2in   =   20._wp
       t3in   =    4._wp
       s1in   =   34._wp
       s2in   =   35.50_wp
       s3in   =   34.75_wp
       ! Initial concentrations in (u/n)mol/kg
       ! Make sure to compile with -DFIXATMPCO2 first
       c1in   = 2100._wp
       c2in   = 2100._wp
       c3in   = 2350._wp
       a1in   = 2300._wp
       a2in   = 2300._wp
       a3in   = 2400._wp
       p1in   =    2._wp
       p2in   =    0._wp
       p3in   =    2.5_wp
       n1in   =   32._wp
       n2in   =    0._wp
       n3in   =   36._wp
       f1in   =    0._wp
       f2in   =    0._wp
       f3in   =    0._wp
       l1in   =    0._wp
       l2in   =    0._wp
       l3in   =    0._wp
       atpco2in =  280._wp

       ! Initial concentrations in (u/n)mol/kg
       ! Here are some equilibrated values run for 100,000 yrs (round-off error notwithstanding)
       ! Make sure to compile without -DFIXATMPCO2
       !c1in   = 2262.29678_wp 
       !c2in   = 2102.93797_wp
       !c3in   = 2363.90418_wp
       !a1in   = 2395.63558_wp
       !a2in   = 2387.42536_wp
       !a3in   = 2399.11408_wp
       !p1in   =    1.82962_wp
       !p2in   =    0.25073_wp
       !p3in   =    2.49856_wp
       !n1in   =   25.31309_wp
       !n2in   =    0.05089_wp
       !n3in   =   36.01617_wp
       !f1in   =    0.00273_wp
       !f2in   =    0.32988_wp
       !f3in   =    0.57430_wp
       !l1in   =    1.70537_wp
       !l2in   =    1.62990_wp
       !l3in   =    1.62589_wp
       !atpco2in= 280.00000_wp

       ! Input some example parameters
       ! Wind speed (m/s)for CO2 gas fluxes
       wind1 = 10._wp
       wind2 = 5._wp
       wind3 = 0._wp
       
       ! Open surface fraction in contact with atmoshpere 
       !  1 => fully open, <1 => flux impeded (e.g. by sea ice)
       fopen1 = 1._wp
       fopen2 = 1._wp
       fopen3 = 0._wp
       
       ! Gamma over lambda for ligands "optimum" value (Lauderdale et al 2020)
       gaovla_opt = 4398._wp
       ! Gamma ligand production rate (in phosphate, not carbon, units)
       gamma_Fe    = 5.e-5_wp*106._wp
       ! Lambda ligand lifetime (s)
       lt_lifetime =  1._wp/((gamma_Fe/106._wp)/gaovla_opt)
       ! Dust deposition in g Fe m-2 year-1
       fe_input1    = 1.5e-3_wp
       fe_input2    = 1.5e-1_wp
       ! Hydrothermal vent input of 1 Gmol/yr (Tagliabue et al., 2010)
       ! mol Fe/yr * g/mol * 1/area  == g Fe m-2 year-1....
       !divide by 2.5e-3 because fe_sol is multiplied again within model.
       fe_input3    = (1.e9_wp*56._wp)/(17.e6_wp*16.e6_wp*2.5e-3_wp)

       ! Biological production maximum rate (mol P/yr)
       alpha_yr     = 6.e-6_wp
       ! Deep ocean box lifetime modifier to capture the gradient due to
       ! photodegradation near the surface and slower loss in the deep
       dlambdadz1   = 1._wp
       dlambdadz2   = 1._wp
       dlambdadz3   = 1.e-2_wp
       ! Overturning rate (m3/s)
       psi         = 20.0e6_wp
       ! File number identifier
       id          = 1
            
       call model(id, maxyears, outputyears, outstepmax,               &
            psi, alpha_yr,                                             &
            gamma_Fe, lt_lifetime,                                     &
            (/dlambdadz1,dlambdadz2,dlambdadz3/),                      &
            (/fe_input1 ,fe_input2 ,fe_input3 /),                      &
            (/wind1,wind2,wind3/),                                     &
            (/fopen1,fopen2,fopen3/),                                  &
            (/t1in,t2in,t3in/),                                        &
            (/s1in,s2in,s3in/),                                        &
            (/c1in,c2in,c3in/),                                        &
            (/a1in,a2in,a3in/),                                        &
            (/p1in,p2in,p3in/),                                        &
            (/n1in,n2in,n3in/),                                        &
            (/f1in,f2in,f3in/),                                        &
            (/l1in,l2in,l3in/),                                        &
            atpco2in,                                                  &
            tout,                                                      &
            thout,                                                     &
            sout,                                                      &
            cout,                                                      &
            aout,                                                      &
            pout,                                                      &
            nout,                                                      &
            fout,                                                      &
            lout,                                                      &
            expout,                                                    &
            nlout,                                                     &
            psout,                                                     &
            pco2out,                                                   &
            atpco2out)

!=======================================================================
       END PROGRAM MICROCOSM_MODEL
!=======================================================================
