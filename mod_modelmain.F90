       MODULE MOD_MODELMAIN
#define NPY_NO_DEPRECATED_API      
       USE MOD_PRECISION
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
       SUBROUTINE MODEL(maxyears, outputyears,                         &
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

! local variables
       include "comdeck.h"
       INTEGER :: i, nstep, outstep
       REAL(KIND=wp) :: time, exp1, exp2, lim, press0, netco2flux
       CHARACTER*64 :: filename
            
! set some parameters
       nstepmax   = int(maxyears*d_per_yr)        
! initialize outstep
       outstep = 1
! initial time
       time = 0.0

! box dimensions (m)
       dx=(/ 17.0e6, 17.0e6, 17.0e6 /) 
       dy=(/ 3.e6, 13.0e6, 16.e6 /) 
       dz=(/  50.0,  50.0, 5050.0 /) 

! Calculate average latitude for each surface box, depending on dy
       lat=(/ 0.0, 0.0, 0.0 /)
       m2deg=180.0/(dy(1)+dy(2))
       lat(1)=-90.0+(dy(1)       /2.0)*m2deg
       lat(2)=-90.0+(dy(1)+(dy(2)/2.0))*m2deg
       lat(3)=-90.0+(dy(3)       /2.0)*m2deg

! box volumes
       area  = dx * dy 
       vol   = area * dz 
       invol = 1.0 / vol   
       
! define array of mixing rates
       K = RESHAPE( (/ 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0 /), &
                   (/ nbox, nbox /) )
       K = K * 1.0e6 
       
       theta = (/t1in,t2in,t3in/)
       salt  = (/s1in,s2in,s3in/) 
! convert from mol kg-1 to moles m-3
       dic = (/c1in,c2in,c3in/)*conv 
       alk = (/a1in,a2in,a3in/)*conv 
       po4 = (/p1in,p2in,p3in/)*conv 
       no3 = (/n1in,n2in,n3in/)*conv 
       fet = (/f1in,f2in,f3in/)*conv
       lt  = (/l1in,l2in,l3in/)*conv 
       sit = (/p1in,p2in,p3in/)*conv*rSIP
       ph  = (/  8 ,  8 ,  8 /)
           
! initialize tracer rates of change
! temp, salt, and si are passive for now, just for co2 solubility
       dthetadt= (/0.0,0.0,0.0/)
       dsaltdt = (/0.0,0.0,0.0/)
       ddicdt= (/0.0,0.0,0.0/)
       dalkdt= (/0.0,0.0,0.0/)
       dpo4dt= (/0.0,0.0,0.0/)
       dno3dt= (/0.0,0.0,0.0/)
       dfetdt= (/0.0,0.0,0.0/)
       dltdt = (/0.0,0.0,0.0/)
       dsitdt = (/0.0,0.0,0.0/)

!! Iron cycle parameters ......... 
! Iron external source
! (gradient based on Parekh et al, 2004; also Mahowald et al 2006)
!  g Fe m-2 year-1 (1 - “Southern Ocean”, 2 - "N Atlantic")
       fe_depo(1)= depfactor * 0.01   
       fe_depo(2)= depfactor * 1.0 
       fe_depo(3)= ventfactor* 1.0 
       
!  convert to mol Fe m-2 s-1
       fe_depo = fe_depo / (weight_fe*s_per_yr) 

! ligand parameters 
! longer lifetime in deep ocean (Ye et al, 2009; Yamaguchi et al, 2002)
       if (lt_lifetime.LE.0.0) then
          lambda = (/ 0.0, 0.0, 0.0 /)
       else
          lambda = (/ (1.0/lt_lifetime),(1.0/lt_lifetime),             &
                   (dlambdadz/lt_lifetime) /) 
       endif
       
! Export production parameters (Parekh et al, 2005):
! max export prodution rate: (again, phosphorus units, mol P m-3 s-1)
!      alpha = 0.5d-6 * conv / (30.0*86400.0) ! Recover with alpha_yr=6e-6
       alpha = alpha_yr * conv / (s_per_yr) 
       
! Initial export production and nutrient limitation code
       export = (/ 0.0, 0.0, 0.0 /)
       lim    = 0.0

! evaluate pstar, consistent with Harvardton Bears SO sensitivity
       pstar = (po4(3) - po4(1)) / po4(3) 
       
! Initialize atmospheric carbon content
! Mass dry atmosphere = (5.1352+/-0.0003)d18 kg (Trenberth & Smith,
! Journal of Climate 2005)
! and Mean molecular mass air = 28.97 g/mol (NASA earth fact sheet)
       atmos_moles = (5.1352e18*1.e3)/28.97
! but need to scale by the ratio of Earth surface area to model surface area
! Earths area = 5.10082000d8 km2 * 1.e6 m2/km2 (NOAA earth fact sheet)
!       atmos_moles = atmos_moles * area(3)/(5.10082e8 * 1.e6)

       pco2atm      = atpco2in * 1.e-6
       atmos_carbon = pco2atm*atmos_moles
       Kg           = 5.0e-5 
       press0       = 0.
       sitM = (po4/conv) * 1.0e6 * rSIP  
       wind = (/5.0,5.0,0.0/)
       fice = (/0.0,0.0,0.0/)

! Find out initial conditions of the carbon system for given input values
       do 100 i = 1,2 
! calculate surface coefficients SETUP_API4PHSWS(temp_in_k, salt, Applied Pressure)
         CALL SETUP_API4PHSWS( 273.15+theta(i), salt(i), press0)
         CALL SETUP_FLUXCOEFFS(273.15+theta(i), salt(i))

! Calculate surface air-sea gas exchange of CO2      
         call carbon_flux(theta(i),salt(i),dic(i),alk(i),po4(i),       &
                         sit(i), ph(i), pco2atm, wind(i),fice(i),      &
                         pco2oce(i), fluxCO2(i))
 100   end do     
 
! write initial values to the output....
! convert to nmol kg-1 for iron, micromol kg-1 for PO4
       dicM = (dic/conv) * 1.0e6  
       alkM = (alk/conv) * 1.0e6  
       po4M = (po4/conv) * 1.0e6  
       no3M = (no3/conv) * 1.0e6  
       fetM = (fet/conv) * 1.0e9
       ltM  = (lt /conv) * 1.0e9
       pco2M= pco2oce    * 1.0e6
       pco2A= pco2atm    * 1.0e6
       
#ifdef WRITEOUTFILE    
! open an output file and write initial values to file
          write (filename, '(a,I0.6,a)') 'microCOSM',id,'.dat'
          open(14,file=filename,status='unknown')
           
! write column header output 
           write(14,*)'t(yr) ',                                        &
                      'THETA(1) THETA(2) THETA(3) ',                   &
                      'SALT(1)  SALT(2)  SALT(3)  ',                   &
                      'DIC(1)   DIC(2)   DIC(3)   ',                   &
                      'ALK(1)   ALK(2)   ALK(3)   ',                   &
                      'PO4(1)   PO4(2)   PO4(3)   ',                   &
                      'NO3(1)   NO3(2)   NO3(3)   ',                   &
                      'FeT(1)   FeT(2)   FeT(3)   ',                   &
                      'LT(1)    LT(2)    LT(3)    ',                   &
                      'Exp(1)   Exp(2)   Limit    ',                   &
                      'P*  OCPCO2(1) OCPCO2(1) ATPCO2'
     
 300       format(1x, i7.1, 1x,                                        &
                  f8.5, 1x, f8.5, 1x, f8.5, 1x,                        &
                  f8.5, 1x, f8.5, 1x, f8.5, 1x,                        &
                  f10.5, 1x, f10.5, 1x, f10.5, 1x,                     &
                  f10.5, 1x, f10.5, 1x, f10.5, 1x,                     &
                  f8.5, 1x, f8.5, 1x, f8.5, 1x,                        &
                  f8.5, 1x, f8.5, 1x, f8.5, 1x,                        &
                  f8.5, 1x, f8.5, 1x, f8.5, 1x,                        &
                  f8.5, 1x, f8.5, 1x, f8.5, 1x,                        &
                  f10.3, 1x, f10.3, 1x, i3.0, 1x,                      &
                  f8.5, 1x, f9.5, 1x, f9.5, 1x, f9.5)

! Write initial conditions to file
           write(14,300) int(time), theta, salt,                       & 
                   dicM, alkM, po4M, no3M, fetM, ltM,                  &
                   export(1), export(2), int(lim), pstar,              &
                   pco2M(1), pco2M(2), pco2A          
#endif

! output to array
       tout  (outstep) = time
       t1out (outstep) = theta(1)
       t2out (outstep) = theta(2)
       t3out (outstep) = theta(3)
       s1out (outstep) = salt(1)
       s2out (outstep) = salt(2)
       s3out (outstep) = salt(3)
       c1out (outstep) = dicM(1)
       c2out (outstep) = dicM(2)
       c3out (outstep) = dicM(3)
       a1out (outstep) = alkM(1)
       a2out (outstep) = alkM(2)
       a3out (outstep) = alkM(3)
       p1out (outstep) = po4M(1)
       p2out (outstep) = po4M(2)
       p3out (outstep) = po4M(3)
       n1out (outstep) = no3M(1)
       n2out (outstep) = no3M(2)
       n3out (outstep) = no3M(3)
       f1out (outstep) = fetM(1)
       f2out (outstep) = fetM(2)
       f3out (outstep) = fetM(3)
       l1out (outstep) = ltM (1)
       l2out (outstep) = ltM (2)
       l3out (outstep) = ltM (3)
       ep1out(outstep) = export(1)*vol(1) 
       ep2out(outstep) = export(2)*vol(2) 
       nlout (outstep) = lim
       psout (outstep) = pstar
       o1pco2out(outstep)= pco2M(1)
       o2pco2out(outstep)= pco2M(2)
       atpco2out(outstep)= pco2A

! Increment outstep
         outstep=outstep+1
         
! timestepping .........................................
       do 200 nstep = 1,INT(nstepmax) 
! evaluate rates of change due to transport
!         dthetadt = transport(nbox, theta, K, psi, invol) 
!         dsaltdt  = transport(nbox, salt,  K, psi, invol) 
         ddicdt = TRANSPORT(nbox, dic, K, psi, invol) 
         dalkdt = TRANSPORT(nbox, alk, K, psi, invol) 
         dpo4dt = TRANSPORT(nbox, po4, K, psi, invol) 
         dno3dt = TRANSPORT(nbox, no3, K, psi, invol) 
         dfetdt = TRANSPORT(nbox, fet, K, psi, invol) 
         dltdt  = TRANSPORT(nbox, lt,  K, psi, invol) 

         time=nstep*dt / (s_per_yr) 

! evaluate biogeochemical rates of change
! Surface boxes...
       netco2flux=0.

       do 210 i = 1,2 
! calculate surface coefficients SETUP_API4PHSWS(temp_in_k, salt, Applied Pressure)
         CALL SETUP_API4PHSWS( 273.15+theta(i), salt(i), press0)
         CALL SETUP_FLUXCOEFFS(273.15+theta(i), salt(i))

! Calculate surface air-sea gas exchange of CO2      
         call carbon_flux(theta(i),salt(i),dic(i),alk(i),po4(i),       &
                         sit(i), ph(i), pco2atm, wind(i),fice(i),      &
                         pco2oce(i), fluxCO2(i))
                
         netco2flux=netco2flux+fluxCO2(i)*area(i)
         
         ddicdt(i) = ddicdt(i) + fluxCO2(i)/dz(i)
         
! biological terms
         ilat=lat(i)
         itim=time*s_per_yr
         ilight=INSOL(itim, s_per_yr, ilat)
         
         ilimit = ilight/(ilight + klight)
         plimit = po4(i)/(po4(i) + kpo4) 
         nlimit = no3(i)/(no3(i) + kno3) 
         flimit = fet(i)/(fet(i) + kfe) 

         export(i) = alpha * ilimit * min(plimit,nlimit,flimit)

! Nutrient Limitation codes:
!  0 = Initial condition
!  1 = macronutrient
!  2 = iron
!  3 = light
!  4 = colimited
! So for box 1 "tens" and box 2 "units" values (-9.0*i+19 allows this)
!   the possible codes are 00, 11-14, 21-24, 31-34 and 41-44, 
!   with the most common being:
!   00 - initial, 
!   11 - global macronutrient limitation, 
!   22 - global iron limitation, and
!   21 - box 1 iron limited and box 2 macronutrient limited.
         if (plimit.LT.flimit.OR.nlimit.LT.flimit                      &
                            .AND.plimit.LT.ilimit) then
!           Phosphate limits production (code 1)       
            lim=lim+(1*(-9.0*i+19))
         elseif (flimit.LT.plimit.OR.flimit.LT.nlimit                  &
                                .AND.flimit.LT.ilimit) then
!           Iron limits production (code 2))
            lim=lim+(2*(-9.0*i+19))
         elseif (ilimit.LT.plimit.AND.ilimit.LT.nlimit                 &
                                 .AND.ilimit.LT.flimit) then
!           Light limits production (but also there will be nutrient limitation too)
!           (code 3)
            lim=lim+(3*(-9.0*i+19))
         else
!           colimitation (code 4)         
            lim=lim+(4*(-9.0*i+19))
         endif

! scale rate of export with rate of phosphorus export
         dpo4dt(i) = dpo4dt(i) - export(i) 
         dno3dt(i) = dno3dt(i) - export(i) * (rCP/rCN) 
         dfetdt(i) = dfetdt(i) - export(i) * (rCP/rCFe) 

! carbonate flux depends on rain ratio
         carb(i)   = export(i) * rCP * rCACO3
         
! For DIC carbonate is the export of 1 mol C (_C_O32-)         
         ddicdt(i) = ddicdt(i) -1.0*carb(i) - export(i)*rCP 
! Whereas for ALK carbonate is the export of 2 mol ions (CO3_2-_) 
!     there is also change in ions due to consumption of nitrate
         dalkdt(i) = dalkdt(i) -2.0*carb(i) + export(i)*rNP

! Dynamic ligand .........
         dltdt(i) = dltdt(i) + (export(i)*gamma_Fe  - lambda(i)*lt(i)) 

 210   end do 
! end of surface boxes loop

! input of iron (can include (vent source)/fe_sol)
       dfetdt = dfetdt + fe_sol * fe_depo / dz 
       
! Deep box...assume whats lost from surface boxes is gained by deep box
       remin = (export(1)*vol(1) + export(2)*vol(2))/vol(3)
       dpo4dt(3) = dpo4dt(3) + remin 
       dno3dt(3) = dno3dt(3) + (rCP/rCN) *remin 
       dfetdt(3) = dfetdt(3) + (rCP/rCFe)*remin 

! carbonate flux depends on rain ratio
       carb(3)   = remin * rCP * rCACO3

       ddicdt(3) = ddicdt(3) + rCP * remin  + 1.0 * carb(3)
       dalkdt(3) = dalkdt(3) - rNP * remin  + 2.0 * carb(3)

! scavenging and complexation of iron
! evaluate local feprime from fet and lt
! determine scavenging rate and add to total tendency
       feprime=FE_EQUIL(nbox, fet, lt, beta)

       dfetdt = dfetdt - Kscav*feprime 

! if FeT > LT, then all excess iron is Fe' and precipitates out quickly
! Liu and Millero indicate very small "soluble" free iron
      do i=1,nbox
        if(fet(i) .gt. lt(i))then
            dfetdt(i) =  dfetdt(i) - (1.0/relaxfe)*(fet(i) - lt(i)) 
        else
            dfetdt(i) =  dfetdt(i) 
        endif
      end do

! Dynamic ligand ..............
      dltdt(3) = dltdt(3) + gamma_Fe*remin - lambda(3)*lt(3)
!..............................

#ifndef FIXATMPCO2
! Update atmospheric CO2 (but only if you want to)

        netco2flux=netco2flux*dt
        call calc_atmos_pco2(atmos_moles,atmos_carbon,netco2flux,      &
                                pco2atm)
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
       pstar = (po4M(3) - po4M(1)) / po4M(3) 
       exp1 = export(1)*vol(1) 
       exp2 = export(2)*vol(2)
       time = nstep*dt / s_per_yr 

! if an output time, write some output to screen and file
       if (mod(time,outputyears) .eq. 0)then

! For output ! convert to nmol kg-1 (Fe, LT) micromol kg-1 (other)
         dicM = (dic/conv) * 1.0e6  
         alkM = (alk/conv) * 1.0e6  
         po4M = (po4/conv) * 1.0e6  
         no3M = (no3/conv) * 1.0e6  
         fetM = (fet/conv) * 1.0e9
         ltM  = (lt /conv) * 1.0e9
         pco2M= pco2oce    * 1.0e6
         pco2A= pco2atm    * 1.0e6

#ifdef WRITEOUTFILE    
! write to file
           write(14,300) int(time), theta, salt,                       & 
                   dicM, alkM, po4M, no3M, fetM, ltM,                  &
                   exp1, exp2, int(lim), pstar,                        &
                   pco2M(1), pco2M(2), pco2A         
#endif
       
! output to array
         tout (outstep)  = time
         t1out (outstep) = theta(1)
         t2out (outstep) = theta(2)
         t3out (outstep) = theta(3)
         s1out (outstep) = salt(1)
         s2out (outstep) = salt(2)
         s3out (outstep) = salt(3)
         c1out (outstep) = dicM(1)
         c2out (outstep) = dicM(2)
         c3out (outstep) = dicM(3)
         a1out (outstep) = alkM(1)
         a2out (outstep) = alkM(2)
         a3out (outstep) = alkM(3)
         p1out (outstep) = po4M(1)
         p2out (outstep) = po4M(2)
         p3out (outstep) = po4M(3)
         n1out (outstep) = no3M(1)
         n2out (outstep) = no3M(2)
         n3out (outstep) = no3M(3)       
         f1out (outstep) = fetM(1)
         f2out (outstep) = fetM(2)
         f3out (outstep) = fetM(3)
         l1out (outstep) = ltM (1)
         l2out (outstep) = ltM (2)
         l3out (outstep) = ltM (3)
         ep1out(outstep) = exp1 
         ep2out(outstep) = exp2
         nlout(outstep)  = lim
         psout(outstep)  = pstar
         o1pco2out(outstep)= pco2M(1)
         o2pco2out(outstep)= pco2M(2)
         atpco2out(outstep)= pco2A
! Python Pandas really doesnt like this output
!         do i=1,nbox
!           pout (outstep,i) = po4M(i)
!           fout (outstep,i) = fetM(i)
!           lout (outstep,i) = ltM (i)
!           epout(outstep,i)= export(i)*vol(i) 
!         enddo

! Increment outstep
         outstep=outstep+1
       endif
       
!      reset lim
       lim = 0.0
       
! end timestepping loop
 200   enddo

#ifdef WRITEOUTFILE    
! close the output file
       close(14)
#endif
       RETURN
       END SUBROUTINE MODEL
!=======================================================================

!=======================================================================
       FUNCTION TRANSPORT(nbox_dum, x, K_dum, psi_dum, invol_dum)
! atmosphere-3-box-ocean carbon cycle model
! evaluate rates of change due to transport
! mick follows, march 2015/ june 2016
       IMPLICIT NONE
       REAL(KIND=wp), DIMENSION(nbox_dum):: TRANSPORT
       
       INTEGER :: nbox_dum
       REAL(KIND=wp), intent(in), DIMENSION(nbox_dum):: x
       REAL(KIND=wp), intent(in), DIMENSION(nbox_dum,nbox_dum):: K_dum
       REAL(KIND=wp), intent(in) :: psi_dum
       REAL(KIND=wp), intent(in), DIMENSION(nbox_dum):: invol_dum
!
       TRANSPORT(1) = invol_dum(1) * (                                 &
            psi_dum*(x(3)-x(1))                                        &
          + K_dum(3,1)*(x(3)-x(1))                                     &
          + K_dum(2,1)*(x(2)-x(1))                                     &
            )
       TRANSPORT(2) = invol_dum(2) * (                                 &
            psi_dum*(x(1)-x(2))                                        &
          + K_dum(1,2)*(x(1)-x(2))                                     &
          + K_dum(3,2)*(x(3)-x(2))                                     &
            )
       TRANSPORT(3) = invol_dum(3) * (                                 &
            psi_dum*(x(2)-x(3))                                        &
          + K_dum(2,3)*(x(2)-x(3))                                     &
          + K_dum(1,3)*(x(1)-x(3))                                     &
            )

       RETURN
       END FUNCTION TRANSPORT
!=======================================================================

!=======================================================================
       FUNCTION INSOL(boxtime,sinayr,boxlat)
! find light as function of date and latitude
! based on paltridge and parson
       IMPLICIT NONE
       REAL(KIND=wp) :: INSOL

       REAL(KIND=wp), intent(in) :: boxtime
       REAL(KIND=wp), intent(in) :: sinayr
       REAL(KIND=wp), intent(in) :: boxlat
! Local variables       
       REAL(KIND=wp) :: dayfrac
       REAL(KIND=wp) :: yday
       REAL(KIND=wp) :: delta
       REAL(KIND=wp) :: sun
       REAL(KIND=wp) :: dayhrs
       REAL(KIND=wp) :: cosz
       REAL(KIND=wp) :: frac
       REAL(KIND=wp) :: fluxi
       REAL(KIND=wp) :: latrad
       REAL(KIND=wp), PARAMETER :: pi = 3.14159265358979323844e0
!fraction of sunlight that is photosynthetically active
       REAL(KIND=wp), PARAMETER :: parfrac = 0.4
!solar constant
       REAL(KIND=wp), PARAMETER :: solar=1360.00 
!planetary albedo
       REAL(KIND=wp), PARAMETER :: albedo=0.60   
             
! find day (****NOTE for year starting in winter*****)
       dayfrac=mod(boxtime,sinayr)/(sinayr) !fraction of year
       yday = 2.0*pi*dayfrac                 !convert to radians
       delta = (0.006918                                               &
           -(0.399912*cos(yday))                                       &
           +(0.070257*sin(yday))                                       &
           -(0.006758*cos(2.0*yday))                                   &
           +(0.000907*sin(2.0*yday))                                   &
           -(0.002697*cos(3.0*yday))                                   &
           +(0.001480*sin(3.0*yday)))                                   

! latitude in radians
       latrad=boxlat*(pi/180.d0)
       sun  = max(-0.999,min(                                          &
            -sin(delta)/cos(delta) * sin(latrad)/cos(latrad)           &
            ,0.999))
       
!       IF (sun1.LE.-0.999) sun1=-0.999
!       IF (sun1.GE. 0.999) sun1= 0.999

       dayhrs = abs(acos(sun))
!      average zenith angle
       cosz = ( sin(delta)*sin(latrad)+                                 &
             ( cos(delta)*cos(latrad)*sin(dayhrs)/dayhrs) )
!       IF (cosz.LE.5.d-3) cosz= 5.d-3
       cosz=max(cosz,5.d-3)
       frac = dayhrs/pi                    !fraction of daylight in day

! daily average photosynthetically active solar radiation just below surface
       fluxi = solar*(1.0-albedo)*cosz*frac*parfrac

! convert to sfac
       INSOL = max(1.d-5,fluxi)
          
       RETURN
       END FUNCTION INSOL
!=======================================================================

!=======================================================================
! solve quadratic for iron speciation
! mick follows, March 2015
       FUNCTION FE_EQUIL(nbox_dum, fet_dum, lt_dum, beta_dum)
       IMPLICIT NONE
       REAL(KIND=wp), DIMENSION(nbox_dum):: FE_EQUIL

       INTEGER,intent(in)  :: nbox_dum
       REAL(KIND=wp), intent(in) , DIMENSION(nbox_dum):: fet_dum
       REAL(KIND=wp), intent(in) , DIMENSION(nbox_dum):: lt_dum 
!!       REAL(KIND=wp), DIMENSION(nbox_dum):: beta_dum
       REAL(KIND=wp), intent(in)  :: beta_dum
! Local variables       
       REAL(KIND=wp), DIMENSION(nbox_dum):: betam
       REAL(KIND=wp), DIMENSION(nbox_dum):: a,b,c,dummy,dummy1,x1,x2
!
       betam = 1.0/beta_dum
       a  = 1.0 
       b  = (lt_dum + betam - fet_dum) 
       c = -1.0 * fet_dum * betam 
! standard quadratic solution for roots
       dummy = b*b - 4.0*a*c 
       dummy1 = dummy**0.5 
       x1 = (-b + dummy1) / (2.0*a) 
       x2 = (-b - dummy1) / (2.0*a) 
! which root?
       FE_EQUIL = x1 
! 
       RETURN
       END FUNCTION FE_EQUIL
!=======================================================================

      END MODULE MOD_MODELMAIN
