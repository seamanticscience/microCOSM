! -*- f90 -*-
       MODULE MOD_MODELMAIN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION 
       USE MOD_PRECISION
       USE MOD_BOXES
       USE MOD_DIMENSIONS
       use MOD_COMMON
       USE MOD_CARBONCHEM
IMPLICIT NONE
! --------------------------------------------------------
! List of (PRIVATE) routines/functions
! --------------------------------------------------------

!       PRIVATE INSOL
!       PRIVATE FE_EQUIL
!       PRIVATE TRANSPORT
            
       CONTAINS

!=======================================================================
       SUBROUTINE MODEL(                                               &
            id,                                                        &
            maxyears,                                                  &
            outputyears,                                               &
            outstepmax,                                                &
            dx,                                                        &
            dy,                                                        &
            dz,                                                        &
            K,                                                         &
            R,                                                         &
            psi,                                                       &
            alpha_yr,                                                  &
            gamma_Fe,                                                  &
            lt_lifetime,                                               &
            dlambdadz,                                                 &
            fe_input,                                                  &
            wind,                                                      &
            fopen,                                                     &
            thin,                                                      &
            sain,                                                      &
            cin,                                                       &
            ain,                                                       &
            pin,                                                       &
            nin,                                                       &
            fin,                                                       &
            lin,                                                       &
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
            ocpco2out,                                                 &
            atpco2out                                                  &        
            )

!-----------------------------------------------------------------------         
! input arguments   
       INTEGER, intent(in)       :: outstepmax, id

       REAL(KIND=wp), intent(in) ::                                    &
            maxyears,                                                  &
            outputyears,                                               &
            psi,                                                       &
            alpha_yr,                                                  &
            gamma_Fe,                                                  & 
            lt_lifetime,                                               &
            atpco2in

       REAL(KIND=wp), intent(in), dimension (nbox) ::                  &
            dx, dy, dz,                                                &
            thin,                                                      &
            sain,                                                      &
            cin,                                                       &
            ain,                                                       &
            pin,                                                       &
            nin,                                                       &
            fin,                                                       &
            lin,                                                       &
            fe_input,                                                  &
            dlambdadz,                                                 &
            wind,                                                      &
            fopen

       REAL(KIND=wp), intent(in), dimension (nbox, nbox) ::            &
            K, R
     
       REAL(KIND=wp), intent(out), dimension (outstepmax,nbox) ::      &
            thout,                                                     &
            sout,                                                      &
            cout,                                                      &
            aout,                                                      &
            pout,                                                      &
            nout,                                                      &
            fout,                                                      &
            lout,                                                      &
            expout,                                                    &
            ocpco2out
            
       REAL(KIND=wp), intent(out), dimension (outstepmax) ::           &
            tout,                                                      &
            psout,                                                     &
            atpco2out

       INTEGER, intent(out), dimension (outstepmax) ::                 &
            nlout                                                     

! local variables
!       include "comdeck.h"
       INTEGER       :: nstep, outstep
       REAL(KIND=wp) :: time
       CHARACTER*64  :: fmt, inifmt, varfmt, frep, filename
!-----------------------------------------------------------------------         
        
       CALL common_assignments()
       
! set some parameters
       nstepmax   = int(maxyears*(speryr/dt))
! initialize outstep
       outstep = 1
! initial time
       time = 0._wp

       CALL establish_dimensions(dx,dy,dz,lat,area,vol,invol,          &
                                          depth,pressure)

       theta = thin
       salt  = sain
! convert from u/nmol kg-1 to moles m-3
       dic = cin * umolkg2molm3
       alk = ain * umolkg2molm3
       po4 = pin * umolkg2molm3
       no3 = nin * umolkg2molm3
       fet = fin * nmolkg2molm3
       lt  = lin * nmolkg2molm3
       sit = pin * umolkg2molm3 * rSIP
       
       ph  = eight
! initialize tracer rates of change
! temp, salt, and si are passive for now, just for co2 solubility
       dthetadt = zero 
       dsaltdt  = zero 
       ddicdt   = zero 
       dalkdt   = zero 
       dpo4dt   = zero 
       dno3dt   = zero 
       dfetdt   = zero 
       dltdt    = zero 
       dsitdt   = zero 

!! Iron cycle parameters ......... 
! Iron external source
!  convert to mol Fe m-2 s-1
       fe_depo = fe_input / (weight_fe*speryr)

! ligand parameters 
! longer lifetime in deep ocean (Ye et al, 2009; Yamaguchi et al, 2002)
       if (lt_lifetime.LE.zero) then
          lambda = zero
       else
          lambda = dlambdadz/lt_lifetime
       endif
       
! Export production parameters (Parekh et al, 2005):
! max export prodution rate: (again, phosphorus units, mol P m-3 s-1)
!      alpha = 0.5d-6 * conv / (30.0*86400.0) ! Recover with alpha_yr=6e-6
       alpha = alpha_yr * conv / (speryr) 
       
! Initial export production and nutrient limitation code
       light  = zero
       ilimit = zero
       plimit = zero
       nlimit = zero
       flimit = zero
       export = zero
       lim    = 0

! evaluate pstar, consistent with Harvardton Bears SO sensitivity
       pstar = calc_pstar(po4) 
       
! Initialize atmospheric carbon content
! Mass dry atmosphere = (5.1352+/-0.0003)d18 kg (Trenberth & Smith,
! Journal of Climate 2005)
! and Mean molecular mass air = 28.97 g/mol (NASA earth fact sheet)
!  but need to scale by the ratio of Earth surface area to model surface area
! Earths area = 5.10082000d8 km2 * 1.e6 m2/km2 (NOAA earth fact sheet)
!       atmos_moles = atmos_moles * area(3)/(5.10082e8 * 1.e6)
       atmos_moles = ((5.1352e18_wp * 1000._wp) / 28.97_wp)            &
                        * (area(3)/(5.10082e8_wp * 1.e6_wp))
       pco2atmos    = atpco2in  * uatm2atm
       atmos_carbon = pco2atmos * atmos_moles
       sitM         = (po4 * rSIP) / umolkg2molm3 

!! Find out initial conditions of the carbon system for given input values
       call carbon_fluxes(theta,                                       &
                           salt,                                       &
                            dic,                                       &
                            alk,                                       &
                            po4,                                       &
                            sit,                                       & 
                             ph,                                       &
                      pco2atmos,                                       & 
                           wind,                                       &
                          fopen,                                       &
                       pressure,                                       &
                      pco2ocean,                                       &
                        fluxCO2 )

! write initial values to the average accumulators....
       timeM  = time
       thetaM = theta
       saltM  = salt
! convert to nmol kg-1 for iron, umol kg-1 for PO4
       dicM   = dic / umolkg2molm3  
       alkM   = alk / umolkg2molm3
       po4M   = po4 / umolkg2molm3
       no3M   = no3 / umolkg2molm3
       fetM   = fet / nmolkg2molm3
       ltM    = lt  / nmolkg2molm3
       exportM= export * vol * molps2gtcyr
       pstarM = pstar
       pco2M  = pco2ocean / uatm2atm
       pco2A  = pco2atmos / uatm2atm
       
#if defined(WRITEOUTFILE)    
! open an output file and write initial values to file
          write (filename, '(a,I0.6,a)') 'microCOSM',id,'.dat'
          open(14,file=filename,status='unknown')

! write column header output 
           write(14,*)'     t(yr)     Limits  ',                       &
               repeat('  THETA    ',nbox),                             &
               repeat(' SALT      ',nbox),                             &
               repeat('DIC        ',nbox),                             &
               repeat('ALK        ',nbox),                             &
               repeat('   PO4     ',nbox),                             &
               repeat('  NO3      ',nbox),                             &
               repeat('   FET     ',nbox),                             &
               repeat('   LIG     ',nbox),                             &
               repeat('   EXPORT  ',nbox),                             &
               '   P*      ',                                          &
               repeat(' OCPCO2    ',nbox),                             &
               ' ATPCO2    '

! Construct fortran format string
! Output the time and nutrient limitation code
           inifmt='1x, i10.1, 1x, i10.0,'
! Each variable then is a space and a 10 position float with 5 decimal places
           varfmt='1x, f10.5'
! This is the number of repeats (10 variables of nbox dimensions plus pstar and atmpco2)
           write(frep ,'(I4)') 10*nbox+2
! Combine everything together
!          fmt='('//trim(fmt)//trim(frep)//'('//trim(varfmt)//'))'
           write(fmt,'(6A)') '(',trim(inifmt),trim(frep),'(',trim(varfmt),'))'


! Write initial conditions to file
           write(14,fmt) int(timeM),                                   & 
                               lim,                                    & 
                            thetaM,                                    & 
                             saltM,                                    &
                              dicM,                                    & 
                              alkM,                                    &
                              po4M,                                    &
                              no3M,                                    &
                              fetM,                                    &
                               ltM,                                    &
                           exportM,                                    &
                            pstarM,                                    &
                             pco2M,                                    &
                             pco2A
#endif

! output to array
       thout     (outstep,1:nbox) = thetaM
       sout      (outstep,1:nbox) = saltM
       cout      (outstep,1:nbox) = dicM
       aout      (outstep,1:nbox) = alkM
       pout      (outstep,1:nbox) = po4M
       nout      (outstep,1:nbox) = no3M
       fout      (outstep,1:nbox) = fetM
       lout      (outstep,1:nbox) = ltM
       expout    (outstep,1:nbox) = exportM
       ocpco2out (outstep,1:nbox) = pco2M
       tout      (outstep) = timeM
       nlout     (outstep) = lim
       psout     (outstep) = pstarM
       atpco2out (outstep) = pco2A
! Increment outstep
       outstep=outstep+1
         
! timestepping .........................................
       do 200 nstep = 1,nstepmax
! evaluate rates of change due to transport
!         dthetadt = transport(nbox, theta, K, psi, invol) 
!         dsaltdt  = transport(nbox, salt,  K, psi, invol) 
         ddicdt = TRANSPORT(dic, K, psi, invol) 
         dalkdt = TRANSPORT(alk, K, psi, invol) 
         dpo4dt = TRANSPORT(po4, K, psi, invol) 
         dno3dt = TRANSPORT(no3, K, psi, invol) 
         dfetdt = TRANSPORT(fet, K, psi, invol) 
         dltdt  = TRANSPORT(lt,  K, psi, invol) 

         time=nstep*dt / (speryr) 

! evaluate biogeochemical rates of change
! Surface boxes...
         netco2flux = zero
                  
! Calculate surface air-sea gas exchange of CO2      
         call carbon_fluxes(theta,                                     &
                             salt,                                     &
                              dic,                                     &
                              alk,                                     &
                              po4,                                     &
                              sit,                                     &
                               ph,                                     &
                        pco2atmos,                                     &
                             wind,                                     &
                            fopen,                                     &
                         pressure,                                     &
                        pco2ocean,                                     &
                          fluxCO2) 
                
         netco2flux = netco2flux + sum(fluxCO2 * area)
         
! Make sure subsurface boxes are masked by fopen = 0         
         ddicdt     = ddicdt     + fluxCO2 / dz
         
! biological terms
         light  = INSOL(time * speryr, lat)
         
         ilimit = light / (light + klight)
         plimit = po4   / (po4   +  kpo4 ) 
         nlimit = no3   / (no3   +  kno3 ) 
         flimit = fet   / (fet   +   kfe ) 

! -ve export is uptake by phytoplankton, +ve export is net remineralization
       bioP = CALC_PRODUCTION(nlimit, plimit, flimit, ilimit, alpha)

       lim  = NUTRIENT_LIMIT_CODE(plimit, nlimit, flimit, ilimit)

! scale rate of nutrient export with rate of phosphorus export
! R matrix determines export flux and remineralization locations
! Spread broadcasts export and volume arrays to matrices
! Each is volume weighted (technically for a single box, this is not necessary,
!    but it works for accumulation of several boxes too.)
        export = CALC_EXPORT(R, bioP, vol, invol)

! carbonate flux depends on rain ratio
       carb = export * rCP * rCACO3
       
       dpo4dt = dpo4dt + export  
       dno3dt = dno3dt + export * (rCP/rCN)
       dfetdt = dfetdt + export * (rCP/rCFe) 
  
! For DIC carbonate is the export of 1 mol C (_C_O32-)   
! -ve bioP is uptake by phytoplankton, +ve bioP is net remineralization
       ddicdt = ddicdt + one * carb + export * rCP 

! Whereas for ALK carbonate is the export of 2 mol ions (CO3_2-_) 
!     there is also change in ions due to consumption of nitrate
       dalkdt = dalkdt + two * carb - export * rNP

! Dynamic ligand production is based on exudation in the surface layers depending on 
!   production and release during remineralization in the ocean interior
       dltdt  = dltdt + (abs(export) * gamma_Fe)  - (lambda * lt) 

! end of surface boxes loop

! input of iron (can include (vent source)/fe_sol)
       dfetdt = dfetdt + fe_sol * fe_depo / dz 

! scavenging and complexation of iron
! evaluate local feprime from fet and lt
! determine scavenging rate and add to total tendency
       feprime=FE_EQUIL(fet, lt, beta)

       dfetdt = dfetdt - Kscav*feprime 

! if FeT > LT, then all excess iron is Fe-prime and precipitates out quickly
! Liu and Millero indicate very small "soluble" free iron
       WHERE (fet > lt ) dfetdt = dfetdt - (one/relaxfe)*(fet - lt) 

#ifndef FIXATMPCO2
! Update atmospheric CO2 (but only if you want to)

        netco2flux=netco2flux*dt
        call calc_atmos_pco2(atmos_moles,                              &
                             atmos_carbon,                             &
                             netco2flux,                               &
                             pco2atmos)
#endif

! Euler forward step concentrations
       theta = theta + dthetadt * dt 
       salt  = salt  + dsaltdt  * dt 
       dic   = dic   + ddicdt   * dt 
       alk   = alk   + dalkdt   * dt 
       po4   = po4   + dpo4dt   * dt
       no3   = no3   + dno3dt   * dt         
       fet   = fet   + dfetdt   * dt 
       lt    = lt    + dltdt    * dt 

! evaluate pstar
       pstar  = calc_pstar(po4) 
       time   = nstep*dt / speryr 

! Increment the average accumulators
       timeM  = (timeM  + time              )
       thetaM = (thetaM + theta             )
       saltM  = (saltM  + salt              )

       dicM   = (dicM + dic / umolkg2molm3  )
       alkM   = (alkM + alk / umolkg2molm3  )
       po4M   = (po4M + po4 / umolkg2molm3  )
       no3M   = (no3M + no3 / umolkg2molm3  )
       fetM   = (fetM + fet / nmolkg2molm3  )
       ltM    = (ltM  + lt  / nmolkg2molm3  )
       pstarM = (pstarM + pstar             )
       pco2M  = (pco2M+ pco2ocean / uatm2atm)
       pco2A  = (pco2A+ pco2atmos / uatm2atm)
         
       exportM= (exportM+export*vol*molps2gtcyr)


! if an output time, write some output to screen and file
       if (mod(time,outputyears) .eq. 0)then
! For output, work out what the average is
         timeM  = timeM  /(outputyears*speryr/dt)
         thetaM = thetaM /(outputyears*speryr/dt)
         saltM  = saltM  /(outputyears*speryr/dt)

         dicM   = dicM   /(outputyears*speryr/dt)
         alkM   = alkM   /(outputyears*speryr/dt)
         po4M   = po4M   /(outputyears*speryr/dt)
         no3M   = no3M   /(outputyears*speryr/dt)
         fetM   = fetM   /(outputyears*speryr/dt)
         ltM    = ltM    /(outputyears*speryr/dt)
         pstarM = pstarM /(outputyears*speryr/dt)
         pco2M  = pco2M  /(outputyears*speryr/dt)
         pco2A  = pco2A  /(outputyears*speryr/dt)
         
         exportM= exportM/(outputyears*speryr/dt)

#if defined(WRITEOUTFILE)           
! Write model state to file                             
           write(14,fmt) int(timeM),                                   & 
                               lim,                                    & 
                            thetaM,                                    & 
                             saltM,                                    &
                              dicM,                                    & 
                              alkM,                                    &
                              po4M,                                    &
                              no3M,                                    &
                              fetM,                                    &
                               ltM,                                    &
                          -exportM,                                    &
                            pstarM,                                    &
                             pco2M,                                    &
                             pco2A
#endif
       
! output to array
       thout     (outstep,1:nbox) = theta
       sout      (outstep,1:nbox) = salt
       cout      (outstep,1:nbox) = dicM
       aout      (outstep,1:nbox) = alkM
       pout      (outstep,1:nbox) = po4M
       nout      (outstep,1:nbox) = no3M
       fout      (outstep,1:nbox) = fetM
       lout      (outstep,1:nbox) = ltM
       expout    (outstep,1:nbox) =-exportM
       ocpco2out (outstep,1:nbox) = pco2M
       tout      (outstep) = timeM
       nlout     (outstep) = lim
       psout     (outstep) = pstarM
       atpco2out (outstep) = pco2A

! Reset the average accumulators
       timeM  = 0._wp
       thetaM = 0._wp
       saltM  = 0._wp

       dicM   = 0._wp
       alkM   = 0._wp
       po4M   = 0._wp
       no3M   = 0._wp
       fetM   = 0._wp
       ltM    = 0._wp
       pstarM = 0._wp
       pco2M  = 0._wp
       pco2A  = 0._wp
         
       exportM= 0._wp

! Increment outstep
         outstep=outstep+1
       endif
! end timestepping loop
 200   enddo

#if defined(WRITEOUTFILE)    
! close the output file
       close(14)
#endif
       RETURN
       END SUBROUTINE MODEL
!=======================================================================

!=======================================================================
! find light as function of date and latitude
! based on paltridge and parson
       FUNCTION INSOL(boxtime,boxlat)
       USE MOD_BOXES

       IMPLICIT NONE
       REAL(KIND=wp), DIMENSION(nbox) :: INSOL
       REAL(KIND=wp), intent(in), DIMENSION(nbox) :: boxlat
       REAL(KIND=wp), intent(in) :: boxtime
! Local variables       
       REAL(KIND=wp), DIMENSION(nbox) :: dayfrac
       REAL(KIND=wp), DIMENSION(nbox) :: yday
       REAL(KIND=wp), DIMENSION(nbox) :: delta
       REAL(KIND=wp), DIMENSION(nbox) :: dayhrs
       REAL(KIND=wp), DIMENSION(nbox) :: frac
       REAL(KIND=wp), DIMENSION(nbox) :: fluxi
       REAL(KIND=wp), DIMENSION(nbox) :: latrad
       REAL(KIND=wp), DIMENSION(nbox) :: sun
       REAL(KIND=wp), DIMENSION(nbox) :: cosz
       REAL(KIND=wp), PARAMETER :: pi = 3.14159265358979323844_wp
       REAL(KIND=wp), PARAMETER :: deg2rad = pi/180._wp
!fraction of sunlight that is photosynthetically active
       REAL(KIND=wp), PARAMETER :: parfrac = 0.4_wp
!solar constant
       REAL(KIND=wp), PARAMETER :: solar = 1360._wp 
!planetary albedo
       REAL(KIND=wp), PARAMETER :: albedo = 0.60_wp   
       REAL(KIND=wp), PARAMETER :: minsun =-0.999_wp   
       REAL(KIND=wp), PARAMETER :: mincosz= 0.005_wp   
       REAL(KIND=wp), PARAMETER :: mininso= 0.00001_wp   

! find day (****NOTE for year starting in winter*****)
       dayfrac=mod(boxtime  ,speryr)/(speryr) !fraction of year
       yday = two*pi*dayfrac                 !convert to radians
       delta = (0.006918_wp                                           &
           -(0.399912_wp*cos(yday))                                   &
           +(0.070257_wp*sin(yday))                                   &
           -(0.006758_wp*cos(two*yday))                               &
           +(0.000907_wp*sin(two*yday))                               &
           -(0.002697_wp*cos(three*yday))                             &
           +(0.001480_wp*sin(three*yday)))                                   

! latitude in radians
       latrad = boxlat*deg2rad
       
       sun    = -sin(delta)/cos(delta) * sin(latrad)/cos(latrad)

       where ( sun .LT. minsun )
          sun = minsun
       elsewhere ( sun .GE. abs(minsun) ) 
          sun =  abs(minsun)
       end where

       dayhrs = abs(acos(sun))
!      average zenith angle
       cosz = ( sin(delta)*sin(latrad)+                                &
              ( cos(delta)*cos(latrad)*sin(dayhrs)/dayhrs) )
       
       where ( cosz .LT. mincosz )
           cosz = mincosz
       end where
       
       frac = dayhrs/pi                       !fraction of daylight in day

! daily average photosynthetically active solar radiation just below surface
       INSOL = solar*(one-albedo)*cosz*frac*parfrac

       where ( INSOL .LT. mininso )
           INSOL = mininso
       end where
       
       RETURN
       END FUNCTION INSOL
!=======================================================================

!=======================================================================
! Calculate surface primary production given macro/micronutrient/light limitation
       FUNCTION CALC_PRODUCTION(nlimit, plimit, flimit, ilimit, alpha)
       USE MOD_BOXES

       IMPLICIT NONE
       REAL(KIND=wp), DIMENSION(nbox):: CALC_PRODUCTION

       REAL(KIND=wp), intent(in) , DIMENSION(nbox):: nlimit
       REAL(KIND=wp), intent(in) , DIMENSION(nbox):: plimit
       REAL(KIND=wp), intent(in) , DIMENSION(nbox):: flimit
       REAL(KIND=wp), intent(in) , DIMENSION(nbox):: ilimit 
       REAL(KIND=wp), intent(in)                  :: alpha

! Non-linear model can use array operations
! minval accepts an array of values and then finds the minimum along dim arguement
!   need to reshape the concatenated nutrient arrays here to stack them by box
! -ve export is uptake by phytoplankton, +ve export is net remineralization
       CALC_PRODUCTION = alpha * ilimit * minval(                      &
                      RESHAPE([ plimit, nlimit, flimit ],[ nbox, 3 ])  &
                                                 ,2)
       RETURN
       END FUNCTION CALC_PRODUCTION
!=======================================================================

!=======================================================================
! calculate export flux of primary production
       FUNCTION CALC_EXPORT(R, bioP, vol, invol)
       USE MOD_BOXES

       IMPLICIT NONE
       REAL(KIND=wp), DIMENSION(nbox):: CALC_EXPORT

       REAL(KIND=wp), intent(in) , DIMENSION(nbox,nbox) :: R
       REAL(KIND=wp), intent(in) , DIMENSION(nbox)      :: bioP
       REAL(KIND=wp), intent(in) , DIMENSION(nbox)      :: vol 
       REAL(KIND=wp), intent(in) , DIMENSION(nbox)      :: invol 

! scale rate of nutrient export with rate of phosphorus export
! R matrix determines export flux and remineralization locations
! Spread broadcasts export and volume arrays to matrices
! Each is volume weighted (technically for a single box, this is not necessary,
!    but it works for accumulation of several boxes too.)
        CALC_EXPORT = SUM(R                                            &
                            * SPREAD(bioP,1,nbox)                      &
                            * SPREAD(vol ,1,nbox)                      &
                          ,2)                                          &
                      * invol
       RETURN
       END FUNCTION CALC_EXPORT
!=======================================================================

!=======================================================================
! solve quadratic for iron speciation
! mick follows, March 2015
       FUNCTION FE_EQUIL(iron, ligand, lig_beta)
       USE MOD_BOXES

       IMPLICIT NONE
       REAL(KIND=wp), DIMENSION(nbox):: FE_EQUIL

       REAL(KIND=wp), intent(in) , DIMENSION(nbox):: iron
       REAL(KIND=wp), intent(in) , DIMENSION(nbox):: ligand 
!!       REAL(KIND=wp), DIMENSION(nbox_dum):: beta_dum
       REAL(KIND=wp), intent(in)  :: lig_beta
! Local variables       
       REAL(KIND=wp), DIMENSION(nbox):: invbeta
       REAL(KIND=wp), DIMENSION(nbox):: a,b,c,discriminant,x1,x2
!
       invbeta = one/lig_beta
       a  = one 
       b  = (ligand + invbeta - iron) 
       c = -one * iron * invbeta 
! standard quadratic solution for roots
       discriminant = ( b*b - four*a*c )**0.5
       x1 = (-b + discriminant) / (two*a) 
       x2 = (-b - discriminant) / (two*a) 
! which root?
       FE_EQUIL = x1 
! 
       RETURN
       END FUNCTION FE_EQUIL
!=======================================================================

!=======================================================================
! solve quadratic for iron speciation
! mick follows, March 2015
       FUNCTION NUTRIENT_LIMIT_CODE(plimit, nlimit, flimit, ilimit)
       USE MOD_BOXES

       IMPLICIT NONE
       INTEGER :: NUTRIENT_LIMIT_CODE

       REAL(KIND=wp), intent(in) , DIMENSION(nbox) :: plimit
       REAL(KIND=wp), intent(in) , DIMENSION(nbox) :: nlimit
       REAL(KIND=wp), intent(in) , DIMENSION(nbox) :: flimit
       REAL(KIND=wp), intent(in) , DIMENSION(nbox) :: ilimit
       REAL(KIND=wp), DIMENSION(nbox, 4)           :: leibig

       INTEGER,                    DIMENSION(nbox) :: lim
       
       INTEGER           :: i, limout
       CHARACTER(nbox*2) :: clim
       CHARACTER(2)      :: tmp
       
       lim=0
! Nutrient Limitation codes:
!  0 = Initial condition
!  1 = phosphate
!  2 = nitrate
!  3 = iron
!  4 = light
       leibig = RESHAPE([ plimit,nlimit,flimit,ilimit ],[ nbox, 4 ])

       lim=minloc(leibig,2)

!      write out array integers and concatenate as a string       
       write(clim,'(I0)') lim(1)
       do i = 2,nbox
          write(tmp ,'(I0)') lim(i)
          write(clim,'(2A)') trim(clim),trim(tmp)
       end do

!      Read the string back in to an integer
       read(clim,*) limout

       NUTRIENT_LIMIT_CODE = limout
       RETURN
       END FUNCTION NUTRIENT_LIMIT_CODE
!=======================================================================

      END MODULE MOD_MODELMAIN
