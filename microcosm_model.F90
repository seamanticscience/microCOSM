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
#if defined(USEJSONIN)
       USE MOD_MODELIO, ONLY: MODELIO_JSON_INPUT
#endif
       IMPLICIT NONE

       CHARACTER*64 :: filename_init            

       INTEGER :: outstepmax, id
       
       REAL(KIND=wp) ::                                                &
            maxyears,                                                  &
            outputyears,                                               &
            m2deg,                                                     &
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
            depth,                                                     &
            latitude,                                                  &
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

       INTEGER(KIND=ip), dimension (:), allocatable   ::               &
            nlout

! Local variable definitions

#if defined(USEJSONIN)
       write (filename_init , '(a,a)') 'microcosm_pickup.json'

       call modelio_json_input(                                       &
            filename_init, id,                                         &
            maxyears, outputyears, outstepmax,                         &
            dx, dy, dz, depth, latitude,                               &
            Kin, Rin, Pin,                                             &
            psi_in, dif_in,                                            &
            alpha_yr, gamma_in, lt_lifein,                             &
            dldz_in, fe_input, wind_in, fopen_in,                      &
            thin, sain, cain, alin, phin, niin, fein, liin,            &
            atpco2in                                                   &
            )          
#else
       ! Input some initial parameters
       maxyears   = 1.e4_wp
       outputyears= 1.e2_wp
       outstepmax = int((maxyears/outputyears)+1)

! Initialize input arguements manually...depends on how many boxes there are
       thin     =      0._wp
       sain     =     34._wp
       cain     =   2150._wp
       alin     =   2350._wp
!       phin     =      1._wp
       phin     =      2._wp
!       niin     =     16._wp
       niin     =     36._wp
       fein     =      0._wp
       liin     =      2._wp
       atpco2in =    280._wp
       
! Overturning and mixing rates (m3/s)
!       psi_in = 20.e6_wp
       dif_in =  1.e6_wp
       
! Wind speed (m/s)for CO2 gas fluxes
       wind_in      =   0._wp
       
! Open surface fraction in contact with atmoshpere 
!  1 => fully open, <1 => flux impeded (e.g. by sea ice)
       fopen_in     =  0._wp
       
! Gamma over lambda for ligands "optimum" value (Lauderdale et al 2020)
!       gaovla_opt   = 0._wp !4398._wp
       gaovla_opt   = 4398._wp

! Gamma ligand production rate (in phosphate, not carbon, units)
!       gamma_in     = 0._wp !5.e-5_wp*106._wp
       gamma_in     = 5.e-5_wp*106._wp
       
! Lambda ligand lifetime (s)
!       lt_lifein    = 0._wp ! 1._wp/((gamma_in/106._wp)/gaovla_opt)
       lt_lifein    = 1._wp/((gamma_in/106._wp)/gaovla_opt)
       
! Dust deposition in g Fe m-2 year-1
       fe_input     =  0._wp

! Biological production maximum rate (mol P/yr)
       alpha_yr      = 6.e-6_wp

! Deep ocean box lifetime modifier to capture the gradient due to
! photodegradation near the surface and slower loss in the deep
       dldz_in       =   0._wp

! File number identifier
       id            = 1

! Array inputs
#if defined(FOURBOX)
! For a 20SV AMOC, psi_in (i.e. Southern Ocean upwelling) needs to be 2x
       psi_in= 40.e6_wp

       dx    = [17.0e6_wp, 17.0e6_wp, 17.0e6_wp, 17.0e6_wp]
       dy    = [ 1.0e6_wp,  3.0e6_wp, 12.0e6_wp, 16.0e6_wp]
       dz    = [50._wp,    50._wp,    50._wp,     5050._wp]
       depth = [       dz(1)/2._wp,                                    &
                       dz(2)/2._wp,                                    &
                       dz(3)/2._wp,                                    &
                 dz(1)+dz(4)/2._wp]
                 
       m2deg    = 180._wp/dy(4)  

       latitude = [(           +(dy(1)/2._wp)),                        &
                   (dy(1)      +(dy(2)/2._wp)),                        &
                   (dy(1)+dy(2)+(dy(3)/2._wp)),                        &
                   (           +(dy(4)/2._wp))                         &
                  ]
       latitude = -90._wp+(latitude*m2deg)

! define arrays (nbox*nbox long) of box connectivity for mixing and overturning (by rows)
! Box 1 mixes with box 2 and 4; 
! Box 2 mixes with box 1, 3 and 4; 
! Box 3 mixes with box 2 and 4; 
! Box 4 mixes with box 1, 2, and 3.
!                       Box1    Box2    Box3    Box4 
       Kin = RESHAPE([ 0.0_wp, 1.0_wp, 0.0_wp, 1.0_wp,                 & ! Box1
                       1.0_wp, 0.0_wp, 1.0_wp, 1.0_wp,                 & ! Box2
                       0.0_wp, 1.0_wp, 0.0_wp, 1.0_wp,                 & ! Box3
                       1.0_wp, 1.0_wp, 1.0_wp, 0.0_wp                  & ! Box4
                     ], [ nbox, nbox ] )
       
! Box 1 is upstream of box 4 (ie AABW downwelling); 
! Box 2 is upstream of box 1 and 3 (ie Antarctic divergence); 
! Box 3 is upstream of box 4 (ie NADW downwelling); 
! Box 4 is upstream of box 2 (ie SO upwelling).
!                       Box1    Box2    Box3    Box4 
       Pin = RESHAPE([ 0.0_wp, 0.0_wp, 0.0_wp, 0.5_wp,                 & ! Box1
                       0.5_wp, 0.0_wp, 0.5_wp, 0.0_wp,                 & ! Box2
                       0.0_wp, 0.0_wp, 0.0_wp, 0.5_wp,                 & ! Box3
                       0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp                  & ! Box4
                     ], [ nbox, nbox ] )

! define array of remineralization coefficients (by rows)
! -1 indicates all of export is lost from cell, while 
! +1 indicates all of export is remineralized (gained) by cell
! Box 1 loses export from Box 1, which is completely remineralized in Box 3 (Box 2 is adjacent)
! Box 2 loses export from Box 2, which is also completely remineralized in Box 3 (Box 1 is adjacent)
!                       Box1    Box2    Box3    Box4 
       Rin = RESHAPE([-1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp,                 & ! Box1
                       0.0_wp,-1.0_wp, 0.0_wp, 1.0_wp,                 & ! Box2
                       0.0_wp, 0.0_wp,-1.0_wp, 1.0_wp,                 & ! Box3
                       0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp                  & ! Box4
                     ], [ nbox, nbox ] )
                     
! Initial conditions
       thin(1:4)= [   -1._wp,    2._wp,   20._wp  ,    4._wp   ]
       sain(1:4)= [   35._wp,   34._wp,   35.50_wp,   34.75_wp ]
       cain(1:4)= [ 2100._wp, 2100._wp, 2100._wp  , 2400._wp   ]
       alin(1:4)= [ 2350._wp, 2300._wp, 2300._wp  , 2400._wp   ]
       phin(1:4)= [    2._wp,    2._wp,    0.0_wp ,    2.5_wp  ]
       niin(1:4)= [   32._wp,   32._wp,    0._wp  ,   36._wp   ]

! Wind speed (m/s)for CO2 gas fluxes
       wind_in(1:4) = [ 10._wp, 10._wp, 5._wp, 0._wp ]
       
! Open surface fraction in contact with atmoshpere 
!  1 => fully open, <1 => flux impeded (e.g. by sea ice)
       fopen_in(1:4)= [ 0.5_wp, 1._wp, 1._wp, 0._wp ]       
       
! Dust deposition in g Fe m-2 year-1
       fe_input(1:4)= [ 1.5e-3_wp, 1.5e-3_wp, 1.5e-1_wp,               &
! Hydrothermal vent input of 1 Gmol/yr (Tagliabue et al., 2010)
! mol Fe/yr * g/mol * 1/area  == g Fe m-2 year-1....
!divide by 2.5e-3 because fe_sol is multiplied again within model.
        (1.e9_wp*56._wp)/(2.5e-3_wp*dx(4)*dy(4)) ]

! Deep ocean box lifetime modifier to capture the gradient due to
! photodegradation near the surface and slower loss in the deep
       dldz_in(1:4)  = [ 1._wp, 1._wp, 1._wp, 1.e-2_wp ]
       
#else

! Default to the three box model
       psi_in = 20.e6_wp
       
       dx    = [ 17.0e6_wp, 17.0e6_wp,   17.0e6_wp ]
       dy    = [  4.0e6_wp, 12.0e6_wp,   16.0e6_wp ]  
       dz    = [ 50.0_wp  , 50.0_wp  , 5050.0_wp   ]

       depth = [ dz(1)/2._wp, dz(2)/2._wp, dz(1)+dz(3)/2._wp]

       m2deg    = 180._wp/dy(3)  
       latitude = [(dy(1)       /2._wp) ,                              &
                   (dy(1)+(dy(2)/2._wp)),                              &
                   (dy(3)       /2._wp)                                &
                   ]
       latitude = -90._wp+(latitude*m2deg)

! Mixing mask array
       Kin = RESHAPE( [ 0._wp, 1._wp, 1._wp,                           &
                        1._wp, 0._wp, 1._wp,                           &
                        1._wp, 1._wp, 0._wp ],                         &
                      [ nbox, nbox ] )
       
! Overturning mask array
       Pin = RESHAPE([ 0._wp, 1._wp, 0._wp,                            &
                       0._wp, 0._wp, 1._wp,                            &
                       1._wp, 0._wp, 0._wp ],                          &
                        [ nbox, nbox ] )
       
! Remineralization coefficients 
! -1 indicates all of export is lost from cell, while 
! +1 indicates all of export is remineralized (gained) by cell
!  the first box (column one) loses export from box 1,
!  the second box (col two) loses export from box 2,
!  and the third box (col three) gains export from boxes 1 and 2 
       Rin = RESHAPE([ -1._wp, 0._wp, 1._wp,                           &
                        0._wp,-1._wp, 1._wp,                           &
                        0._wp, 0._wp, 0._wp ],                         &
                      [ nbox, nbox ] )

! Initial conditions
       thin(1:3)= [    2._wp,   20._wp  ,    4._wp   ]
       sain(1:3)= [   34._wp,   35.50_wp,   34.75_wp ]
       cain(1:3)= [ 2200._wp, 2100._wp  , 2400._wp   ]
       alin(1:3)= [ 2350._wp, 2350._wp  , 2400._wp   ]
       phin(1:3)= [    2._wp,    0._wp  ,    2.5_wp  ]
       niin(1:3)= [   25._wp,    0._wp  ,   35._wp   ]

! Initial concentrationesix in (u/n)mol/kg
! Here are some equilibrated values run for 100,000 yrs (round-off error notwithstanding)
! Make sure to compile without -DFIXATMPCO2
!       cain(1:3)= [ 2263.27105_wp, 2104.02729_wp, 2358.11830_wp ]
!       alin(1:3)= [ 2396.24755_wp, 2388.24068_wp, 2399.60156_wp ]
!       phin(1:3)= [ 1.85304_wp   , 0.31325_wp   , 2.49804_wp    ]
!       niin(1:3)= [ 24.68043_wp  , 0.04392_wp   , 35.00046_wp   ]
!       fein(1:3)= [ 0.01007_wp   , 0.37382_wp   , 0.55782_wp    ]
!       liin(1:3)= [ 1.62217_wp   , 1.58451_wp   , 1.57992_wp    ]

! Wind speed (m/s)for CO2 gas fluxes
       wind_in(1:3) = [ 10._wp, 5._wp, 0._wp ]
       
! Open surface fraction in contact with atmoshpere 
!  1 => fully open, <1 => flux impeded (e.g. by sea ice)
       fopen_in(1:3)= [ 1._wp, 1._wp, 0._wp ]

! Dust deposition in g Fe m-2 year-1
       fe_input(1:3)= [ 1.5e-3_wp, 1.5e-1_wp,                          &
! Hydrothermal vent input of 1 Gmol/yr (Tagliabue et al., 2010)
! mol Fe/yr * g/mol * 1/area  == g Fe m-2 year-1....
!divide by 2.5e-3 because fe_sol is multiplied again within model.
        (1.e9_wp*56._wp)/(2.5e-3_wp*dx(3)*dy(3)) ]
! Deep ocean box lifetime modifier to capture the gradient due to
! photodegradation near the surface and slower loss in the deep
       dldz_in(1:3)  = [ 1._wp, 1._wp, 1.e-2_wp ]
#endif        
#endif

       ! allocate memory for output
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

       call model(id, maxyears, outputyears, outstepmax,               &
            dx, dy, dz, depth, latitude,                               &
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
            )

!=======================================================================
       END PROGRAM MICROCOSM_MODEL
!=======================================================================
