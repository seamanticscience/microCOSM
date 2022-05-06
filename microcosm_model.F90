! -*- f90 -*-
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

       INTEGER :: outstepmax, id
       
       REAL(KIND=wp) ::                                                &
            maxyears,                                                  &
            outputyears,                                               &
            gaovla_opt,                                                &
            gamma_Fe,                                                  &
            lt_lifetime,                                               &
            alpha_yr,                                                  &
            psi,                                                       &
            atpco2in 

! Input arrays (nbox dimensions)
       REAL(KIND=wp), dimension (nbox) ::                              & 
            dx,                                                        &
            dy,                                                        &
            dz,                                                        &
            thin,                                                       & 
            sain,                                                       &
            cain,                                                       & 
            alin,                                                       & 
            phin,                                                       & 
            niin,                                                       & 
            fein,                                                       & 
            liin,                                                       & 
            fe_input,                                                  &
            dlambdadz,                                                 &
            wind,                                                      &
            fopen

       REAL(KIND=wp), dimension (nbox, nbox) ::                        & 
            K,                                                         &
            R
            
! Output arrays (nbox, by timestep dimensions)
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
       maxyears   = 1.e4_wp
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

       ! Geometry and array inputs
       dx   = [ 17.0e6_wp, 17.0e6_wp,   17.0e6_wp ]
       dy   = [  4.0e6_wp, 12.0e6_wp,   16.0e6_wp ]  
       dz   = [ 50.0_wp  , 50.0_wp  , 5050.0_wp   ]

       ! define array of mixing rates
       K = RESHAPE( [  0.e6_wp, 1.e6_wp, 1.e6_wp,                       &
                       1.e6_wp, 0.e6_wp, 1.e6_wp,                       &
                       1.e6_wp, 1.e6_wp, 0.e6_wp ],                     &
                       [ nbox, nbox ] )
       
       ! define array of remineralization coefficients (Columnwise)
       ! -1 indicates all of export is lost from cell, while 
       ! +1 indicates all of export is remineralized (gained) by cell
       ! Thus [-1.0, 0.0, 1.0,                                   &
       !        0.0,-1.0, 1.0,                                   &
       !        0.0, 0.0, 0.0 ],  
       ! indicates the first box (column one) loses export from box 1,
       !           the second box (col two) loses export from box 2,
       !       and the third box (col three) gains export from boxes 1 and 2 
       R = RESHAPE([ -1._wp, 0._wp, 1._wp,                              &
                      0._wp,-1._wp, 1._wp,                              &
                      0._wp, 0._wp, 0._wp ],                            &
                      [ nbox, nbox ] )
                    
       ! Initialize input arguements, and set some reasonable values
       thin     =      0._wp
       sain     =     34._wp
       cain     =   2150._wp
       alin     =   2350._wp
       phin     =      1._wp
       niin     =     16._wp
       fein     =      0._wp
       liin     =      0._wp
       
       thin(1:3)= [    2._wp,   20._wp  ,    4._wp   ]
       sain(1:3)= [   34._wp,   35.50_wp,   34.75_wp ]
       ! Initial concentrations in (u/n)mol/kg
       ! Make sure to compile with -DFIXATMPCO2 first
       cain(1:3)= [ 2100._wp, 2100._wp  , 2350._wp   ]
       alin(1:3)= [ 2300._wp, 2300._wp  , 2400._wp   ]
       phin(1:3)= [    2._wp,    0._wp  ,    2.5_wp  ]
       niin(1:3)= [   32._wp,    0._wp  ,   36._wp   ]

       atpco2in =    280._wp

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
       wind         =    0._wp
       wind(1:3)    = [ 10._wp, 5._wp, 0._wp ]
       
       ! Open surface fraction in contact with atmoshpere 
       !  1 => fully open, <1 => flux impeded (e.g. by sea ice)
       fopen        =   0._wp 
       fopen(1:3)   = [ 1._wp, 1._wp, 0._wp ]
       
       ! Gamma over lambda for ligands "optimum" value (Lauderdale et al 2020)
       gaovla_opt   = 4398._wp
       ! Gamma ligand production rate (in phosphate, not carbon, units)
       gamma_Fe     =   5.e-5_wp*106._wp
       ! Lambda ligand lifetime (s)
       lt_lifetime  =   1._wp/((gamma_Fe/106._wp)/gaovla_opt)
       ! Dust deposition in g Fe m-2 year-1
       fe_input     =   0._wp
       fe_input(1:3)= [ 1.5e-3_wp, 1.5e-1_wp,                          &
       ! Hydrothermal vent input of 1 Gmol/yr (Tagliabue et al., 2010)
       ! mol Fe/yr * g/mol * 1/area  == g Fe m-2 year-1....
       !divide by 2.5e-3 because fe_sol is multiplied again within model.
        (1.e9_wp*56._wp)/(17.e6_wp*16.e6_wp*2.5e-3_wp) ]

       ! Biological production maximum rate (mol P/yr)
       alpha_yr      = 6.e-6_wp
       ! Deep ocean box lifetime modifier to capture the gradient due to
       ! photodegradation near the surface and slower loss in the deep
       dlambdadz     =   1._wp
       dlambdadz(1:3)= [ 1._wp, 1._wp, 1.e-2_wp ]
       ! Overturning rate (m3/s)
       psi           = 20.0e6_wp
       ! File number identifier
       id            = 1
            
       call model(id, maxyears, outputyears, outstepmax,               &
            dx, dy, dz,                                                &
            K, R,                                                      &
            psi, alpha_yr,                                             &
            gamma_Fe, lt_lifetime,                                     &
            dlambdadz,                                                 &
            fe_input,                                                  &
            wind,                                                      &
            fopen,                                                     &
            thin,                                                      &
            sain,                                                      &
            cain,                                                      &
            alin,                                                      &
            phin,                                                      &
            niin,                                                      &
            fein,                                                      &
            liin,                                                      &
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
            atpco2out                                                  &           
            )

!=======================================================================
       END PROGRAM MICROCOSM_MODEL
!=======================================================================
