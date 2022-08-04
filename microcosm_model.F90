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
#if defined(USEDUALNUMAD)
       USE DNADMOD
#endif
       USE MOD_PRECISION
       USE MOD_BOXES
       USE MOD_MODELMAIN

       IMPLICIT NONE

       INTEGER :: outstepmax, id
       
       REAL(KIND=wp) ::                                                &
            maxyears,                                                  &
            outputyears,                                               &
            gaovla_opt,                                                &
            gamma_in,                                                  &
            lt_lifein,                                                 &
            alpha_yr,                                                  &
            atpco2in,                                                  &
            psi_in,                                                    &
            dif_in

! Input arrays (nbox dimensions)
       REAL(KIND=wp), dimension (nbox) ::                              & 
            dx,                                                        &
            dy,                                                        &
            dz,                                                        &
            thin,                                                      & 
            sain,                                                      &
            cain,                                                      & 
            alin,                                                      & 
            phin,                                                      & 
            niin,                                                      & 
            fein,                                                      & 
            liin,                                                      & 
            fe_input,                                                  &
            dldz_in,                                                   &
            wind_in,                                                   &
            fopen_in

       REAL(KIND=wp), dimension (nbox, nbox) ::                        & 
            Kin,                                                       &
            Rin,                                                       &
            Pin
            
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

#if defined(USEDUALNUMAD)
       REAL(KIND=wp), dimension (:,:,:), allocatable ::                &
            thdxout,                                                   &
            sdxout,                                                    &
            cdxout,                                                    &
            adxout,                                                    &
            pdxout,                                                    &
            ndxout,                                                    &
            fdxout,                                                    &
            ldxout,                                                    &
            expdxout,                                                  &
            pco2dxout

       REAL(KIND=wp), dimension (:,:), allocatable   ::                &
            tdxout,                                                    &
            psdxout,                                                   &
            atpco2dxout
#endif
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

#if defined(USEDUALNUMAD)
       ! allocate memory
       allocate ( tdxout      (outstepmax,ndv) )
       allocate ( psdxout     (outstepmax,ndv) )
       allocate ( atpco2dxout (outstepmax,ndv) )
       allocate ( thdxout     (outstepmax,nbox,ndv) )
       allocate ( sdxout      (outstepmax,nbox,ndv) )
       allocate ( cdxout      (outstepmax,nbox,ndv) )
       allocate ( adxout      (outstepmax,nbox,ndv) )
       allocate ( pdxout      (outstepmax,nbox,ndv) )
       allocate ( ndxout      (outstepmax,nbox,ndv) )
       allocate ( fdxout      (outstepmax,nbox,ndv) )
       allocate ( ldxout      (outstepmax,nbox,ndv) )
       allocate ( expdxout    (outstepmax,nbox,ndv) )
       allocate ( pco2dxout   (outstepmax,nbox,ndv) )
#endif

       ! Geometry and array inputs
       dx   = [ 17.0e6_wp, 17.0e6_wp,   17.0e6_wp ]
       dy   = [  4.0e6_wp, 12.0e6_wp,   16.0e6_wp ]  
       dz   = [ 50.0_wp  , 50.0_wp  , 5050.0_wp   ]
       
       ! Overturning and mixing rates (m3/s)
       psi_in = 20.e6_wp
       dif_in =  1.e6_wp

       ! mixing mask array
       Kin = RESHAPE( [ 0._wp, 1._wp, 1._wp,                           &
                        1._wp, 0._wp, 1._wp,                           &
                        1._wp, 1._wp, 0._wp ],                         &
                      [ nbox, nbox ] )

       ! Overturning mask array
       Pin = RESHAPE([ 0._wp, 1._wp, 0._wp,                            &
                       0._wp, 0._wp, 1._wp,                            &
                       1._wp, 0._wp, 0._wp ],                          &
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
       Rin = RESHAPE([ -1._wp, 0._wp, 1._wp,                           &
                        0._wp,-1._wp, 1._wp,                           &
                        0._wp, 0._wp, 0._wp ],                         &
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
       cain(1:3)= [ 2200._wp, 2100._wp  , 2400._wp   ]
       alin(1:3)= [ 2350._wp, 2350._wp  , 2400._wp   ]
       phin(1:3)= [    2._wp,    0._wp  ,    2.5_wp  ]
       niin(1:3)= [   25._wp,    0._wp  ,   35._wp   ]

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

       cain(1:3)= [ 2263.27105_wp, 2104.02729_wp, 2358.11830_wp ]
       alin(1:3)= [ 2396.24755_wp, 2388.24068_wp, 2399.60156_wp ]
       phin(1:3)= [ 1.85304_wp   , 0.31325_wp   , 2.49804_wp    ]
       niin(1:3)= [ 24.68043_wp  , 0.04392_wp   , 35.00046_wp   ]
       fein(1:3)= [ 0.01007_wp   , 0.37382_wp   , 0.55782_wp    ]
       liin(1:3)= [ 1.62217_wp   , 1.58451_wp   , 1.57992_wp    ]
       atpco2in = 280.00000_wp

       ! Input some example parameters
       ! Wind speed (m/s)for CO2 gas fluxes
       wind_in      =    0._wp
       wind_in(1:3) = [ 10._wp, 5._wp, 0._wp ]
       
       ! Open surface fraction in contact with atmoshpere 
       !  1 => fully open, <1 => flux impeded (e.g. by sea ice)
       fopen_in     =   0._wp
       fopen_in(1:3)= [ 1._wp, 1._wp, 0._wp ]
       
       ! Gamma over lambda for ligands "optimum" value (Lauderdale et al 2020)
       gaovla_opt   = 4398._wp
       ! Gamma ligand production rate (in phosphate, not carbon, units)
       gamma_in     = 5.e-5_wp*106._wp
       ! Lambda ligand lifetime (s)
       lt_lifein    = 1._wp/((gamma_in/106._wp)/gaovla_opt)
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
       dldz_in       =   1._wp
       dldz_in(1:3)  = [ 1._wp, 1._wp, 1.e-2_wp ]

       ! File number identifier
       id            = 1
            
       call model(id, maxyears, outputyears, outstepmax,               &
            dx, dy, dz,                                                &
            Kin, Rin, Pin,                                             &
            psi_in, dif_in,                                            &
            alpha_yr, gamma_in, lt_lifein,                             &
            dldz_in,                                                   &
            fe_input,                                                  &
            wind_in,                                                   &
            fopen_in,                                                  &
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
#if defined(USEDUALNUMAD)
            ,tdxout,                                                   &
            thdxout,                                                   &
            sdxout,                                                    &
            cdxout,                                                    &
            adxout,                                                    &
            pdxout,                                                    &
            ndxout,                                                    &
            fdxout,                                                    &
            ldxout,                                                    &
            expdxout,                                                  &
            psdxout,                                                   &
            pco2dxout,                                                 &
            atpco2dxout                                                &
#endif
            )

!=======================================================================
       END PROGRAM MICROCOSM_MODEL
!=======================================================================
