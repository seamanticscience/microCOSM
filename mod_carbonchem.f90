       MODULE MOD_CARBONCHEM
! Import modules from the SOLVESAPHE package of Munhoven (2013)
!
!    Copyright 2013 Guy Munhoven
!
!    SolveSAPHE is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SolveSAPHE is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with SolveSAPHE.  If not, see <http://www.gnu.org/licenses/>.
!
       USE MOD_PRECISION
       USE MOD_CHEMCONST
       USE MOD_PHSOLVERS
       USE MOD_CHEMSPECIATION
       USE MOD_BOXES
       IMPLICIT NONE
      
! --------------------------------------------------------
! List of (PRIVATE) routines/functions
! --------------------------------------------------------

!       PRIVATE AK_FFATM_WEIS80
!       PRIVATE AK_FFOCE_WEIS74
!       PRIVATE AK_CARB_0_WEIS74

!       REAL(KIND=wp), PARAMETER :: gasconst_bar_cm3_o_mol_k = 83.14510_wp ! libthdyct
!REAL(KIND=wp), PARAMETER, PRIVATE :: gasconst_bar_cm3_o_mol_k = 83.14472_wp ! Handbook (2007)

! 0 degrees centigrade in Kelvin
!       REAL(KIND=wp), PARAMETER :: t_k_zerodegc = 273.15_wp ! Handbook (2007)

! Pressure at one atmosphere (bar)
       REAL(KIND=wp), PARAMETER :: p_bar_oneatmosphere = 1.01325_wp ! Handbook (2007)

! convert from mol/kg to mol/m3
       REAL(KIND=wp), PARAMETER :: permil = 1.0_wp/1024.5_wp
       
! --------------------------------------------------------------
! variable for usage by users of the module
! --------------------------------------------------------------
       REAL(KIND=wp) :: apiff_atm, apiff_oce
       REAL(KIND=wp) :: api0_dic

       CONTAINS

!=======================================================================
       SUBROUTINE SETUP_FLUXCOEFFS(t, s)
! ------------------
! Argument variables
! ------------------
!     t_k    : temperature in Kelvin
!     s      : salinity
       REAL(KIND=wp), INTENT(IN) :: t
       REAL(KIND=wp), INTENT(IN) :: s

! ---------------
! Local variables
! ---------------
       REAL(KIND=wp) :: t_k

       t_k=t_k_zerodegc+t

       apiff_atm =  AK_FFATM_WEIS80 (t_k, s)
       apiff_oce =  AK_FFOCE_WEIS74 (t_k, s)
       api0_dic  =  AK_CARB_0_WEIS74(t_k, s)
       RETURN
       END SUBROUTINE SETUP_FLUXCOEFFS
!=======================================================================
      
!=======================================================================
       FUNCTION AK_FFATM_WEIS80(t_k, s)
!      Calculate f = k0(1-pH2O)*correction term for non-ideality
!                in (mol/kg-SW)/atmosphere
!      References: Weiss & Price (1980, Mar. Chem., 8, 347-359
!                  Eq 13 with table 6 values)
! Note      : currently no pressure correction
       IMPLICIT NONE

       REAL(KIND=wp) :: AK_FFATM_WEIS80
! ------------------
! Argument variables
! ------------------
!     s      : salinity
!     t_k    : temperature in K
       REAL(KIND=wp), INTENT(IN) :: t_k
       REAL(KIND=wp), INTENT(IN) :: s

! ---------------
! Local variables
! ---------------
!     zt_k_o_100   : zt_k/100
       REAL(KIND=wp) :: t_k_o_100
       REAL(KIND=wp) :: t_k_o_100_2

       t_k_o_100   = t_k/100._wp
       t_k_o_100_2 = t_k_o_100*t_k_o_100

       AK_FFATM_WEIS80                                                 &
          = exp(-162.8301_wp + 218.2968_wp/t_k_o_100                   &
               + 90.9241_wp*log(t_k_o_100) - 1.47696_wp*t_k_o_100_2    &
               + s * (.025695_wp - .025225_wp*t_k_o_100                &
               + 0.0049867_wp*t_k_o_100_2))

       RETURN
       END FUNCTION AK_FFATM_WEIS80
!=======================================================================

!=======================================================================
       FUNCTION AK_FFOCE_WEIS74(t_k, s)
!      Calculate Fugacity Factor needed for non-ideality in ocean
!                 in [(mol/kg-SW)/atm]
!      References: Weiss (1974) Marine Chemistry
!      pH scale  : N/A
! Note      : currently no pressure correction
       IMPLICIT NONE
       REAL(KIND=wp) :: AK_FFOCE_WEIS74
! ------------------
! Argument variables
! ------------------
!     s      : salinity
!     t_k    : temperature in K
       REAL(KIND=wp), INTENT(IN) :: t_k
       REAL(KIND=wp), INTENT(IN) :: s

! ---------------
! Local variables
! ---------------
!     zt_k_o_100   : zt_k/100
       REAL(KIND=wp) :: zt_k_o_100
       REAL(KIND=wp) :: delta
       REAL(KIND=wp) :: B1
       REAL(KIND=wp) :: B

       zt_k_o_100 = t_k/100.0_wp
              
       delta = (57.7_wp - 0.118_wp*t_k)
       B1 = -1636.75_wp + 12.0408_wp*t_k - 0.0327957_wp*t_k*t_k
       B  = B1 + 3.16528_wp*t_k*t_k*t_k*1.e-5_wp
 
!  "x2" term often neglected (assumed=1) in applications of Weiss's (1974) eq.9
!   x2 = 1 - x1 = 1 - xCO2 (it is very close to 1, but not quite)
       AK_FFOCE_WEIS74                                                 &
            = exp( (B+2.0_wp*delta) * p_bar_oneatmosphere /            &
                         (gasconst_bar_cm3_o_mol_k*t_k))

       RETURN
       END FUNCTION AK_FFOCE_WEIS74
!=======================================================================

!!=======================================================================
!       FUNCTION AK_CARB_0_WEIS74(t_k, s)
!
!! Function calculates K0 in (mol/kg-SW)/atmosphere
!
!! References: Weiss (1979) [(mol/kg-SW)/atm]
!! pH scale  : N/A
!! Note      : currently no pressure correction
!       IMPLICIT NONE
!       REAL(KIND=wp) :: AK_CARB_0_WEIS74
!
!! ------------------
!! Argument variables
!! ------------------
!!     s      : salinity
!!     t_k    : temperature in K
!       REAL(KIND=wp), INTENT(IN) :: t_k
!       REAL(KIND=wp), INTENT(IN) :: s
!
!! ---------------
!! Local variables
!! ---------------
!!     zt_k_o_100   : zt_k/100
!       REAL(KIND=wp) :: zt_k_o_100
!
!       zt_k_o_100 = t_k/100._wp
!
!       AK_CARB_0_WEIS74                                                &
!          = EXP( -60.2409_wp + 93.4517_wp/zt_k_o_100                   &
!                + 23.3585_wp*LOG(zt_k_o_100)                           &
!                + (   0.023517_wp - 0.023656_wp*zt_k_o_100             &
!                   + 0.0047036_wp*zt_k_o_100*zt_k_o_100)*s )           
!
!       RETURN
!       END FUNCTION AK_CARB_0_WEIS74
!!=======================================================================

!=======================================================================
      SUBROUTINE CALC_PCO2(theta,salt,dic,alk,po4,sit,                 &
                         ph,pco2oce,co3,hco3,co2aq)

!General parameters
!      dictot = total inorgani!carbon (mol/m^3)
!            where 1 T = 1 metri!ton = 1000 kg
!      alktot = total alkalinity (eq/m^3)
!      po4tot = inorganic phosphate (mol/^3)
!      siltot = inorganic silicate (mol/^3)
!      t   = temperature (degrees C)
!      s   = salinity (PSU)
      REAL(KIND=wp), INTENT(IN)    :: theta, salt
      REAL(KIND=wp), INTENT(IN)    :: dic, alk, po4, sit
      REAL(KIND=wp), INTENT(INOUT) :: ph
      REAL(KIND=wp), INTENT(OUT)   :: pco2oce, co3, hco3, co2aq

!    Local variables
      REAL(KIND=wp) hini, z_val, hnew
      REAL(KIND=wp) bor, nh4, h2s, so4, flu

      nh4 =  0.E-3_wp
      h2s =  0.E-3_wp
      bor =  A_BTOT_SALIN  (salt)
      so4 =  A_SO4TOT_SALIN(salt)
      flu =  A_FTOT_SALIN  (salt)

      hini = 10._wp**(-1.0_wp * ph)

      hnew = SOLVE_AT_GENERAL(alk, dic, bor, po4, sit, nh4, h2s,       &   
                            so4, flu, p_hini=hini, p_val=z_val)
     
!Return update pH to main routine
      ph = -log10(hnew)

!now determine [CO2*], HCO3- and CO32- , carbonate ion concentration
      CALL SPECIATION_DIC(dic,hnew,co2aq,hco3,co3)

! get K0 and fugacity            
      CALL SETUP_FLUXCOEFFS(theta, salt)

!      z_fco2 = z_co2aq/z_k0
!     co2aq will be in mol/m3, convert to mol/kg to use with K0 and ff
      pco2oce = co2aq*permil/(api0_dic*apiff_oce)
      
      RETURN
      END SUBROUTINE CALC_PCO2
!=======================================================================

!=======================================================================
      SUBROUTINE CARBON_FLUXES(theta,salt,dic,alk,po4,si,              &
                               ph,pco2atmos,wind,fopen,pressure,       &
                               pco2ocean,fluxCO2)
!Calculate air-sea CO2 flux

!    Argument variables
      REAL(KIND=wp), INTENT(IN),    DIMENSION(nbox) :: theta,          &
                                                        salt,          &
                                                         dic,          &
                                                         alk,          &
                                                         po4,          &
                                                         si
      REAL(KIND=wp), INTENT(IN),    DIMENSION(nbox) :: wind,           &
                                                       fopen,          &
                                                       pressure
      REAL(KIND=wp), INTENT(INOUT), DIMENSION(nbox) :: ph
      REAL(KIND=wp), INTENT(OUT),   DIMENSION(nbox) :: pco2ocean,      &
                                                       fluxCO2
      REAL(KIND=wp), INTENT(IN)                     :: pco2atmos
          
!    Local variables       
      REAL(KIND=wp),   DIMENSION(nbox) :: scadic, Kwexch
      REAL(KIND=wp),   DIMENSION(nbox) :: co3, hco3, co2aq
      
      INTEGER                          :: i
      
! Initialize
      pco2ocean = 0._wp
      fluxCO2   = 0._wp
      scadic    = 0._wp
      Kwexch    = 0._wp
      co3       = 0._wp
      hco3      = 0._wp
      co2aq     = 0._wp
           
!     calculate SCHMIDT NO. for CO2 (4th order, Wanninkhof 2014)
      scadic =  2116.8_wp                                              &
              -  136.25_wp     * theta                                 &
              +    4.7353_wp   * theta * theta                         &
              -    0.092307_wp * theta * theta * theta                 &
              +    7.555e-4_wp * theta * theta * theta * theta

      Kwexch =  (0.337_wp * wind**2._wp/3.6e5_wp) * fopen              &
                / sqrt(scadic/660._wp)    
               
      do i = 1,nbox
          ! calculate surface coefficients 
          CALL SETUP_API4PHSWS( t_k_zerodegc + theta(i),               &
                                                salt(i),               &
                                            pressure(i))
          CALL SETUP_FLUXCOEFFS(t_k_zerodegc + theta(i),               & 
                                                salt(i))

          CALL CALC_PCO2(theta(i),                                     &
                          salt(i),                                     &
                           dic(i),                                     &
                           alk(i),                                     &
                           po4(i),                                     &
                            si(i),                                     &
                            ph(i),                                     &
                     pco2ocean(i),                                     &
                           co3(i),                                     &
                          hco3(i),                                     &
                         co2aq(i))
           
!Flux = kw*rho*(ff*pCO2atm-k0*FugFac*pCO2ocean)
          fluxCO2(i) = Kwexch(i)*(                                     &
                               pco2atmos    * apiff_atm -              &
                               pco2ocean(i) * apiff_oce                &
                             * api0_dic )
                  
      end do

!      convert flux (mol kg-1 m s-1) to (mol m-2 s-1)
      fluxCO2 = fluxCO2/permil

      RETURN
      END SUBROUTINE CARBON_FLUXES
!=======================================================================

!=======================================================================
      SUBROUTINE CALC_ATMOS_PCO2(atmos_moles,atmos_carbon,netCO2flux,  &
                                pco2atm)
!Integrate atmospheri!carbon content       

!    Argument variables
      REAL(KIND=wp), INTENT(IN)   :: atmos_moles, netCO2flux
      REAL(KIND=wp), INTENT(INOUT):: pco2atm, atmos_carbon
      
!    How much carbon (moles) is currently in the atmosphere
!     atmos_carbon=atmos_moles*pco2atm

!    What is the change resulting from air-sea flux in surface boxes
      atmos_carbon=atmos_carbon-netCO2flux
      
!    Update atmospheric CO2 level
      pco2atm = atmos_carbon/atmos_moles

      RETURN
      END SUBROUTINE CALC_ATMOS_PCO2    
!=======================================================================

      END MODULE MOD_CARBONCHEM
