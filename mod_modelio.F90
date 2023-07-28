! -*- f90 -*-
       MODULE MOD_MODELIO
#if defined(USEJSONOUT)
       USE JSON_MODULE
#endif
       USE MOD_PRECISION
       USE MOD_BOXES
       USE MOD_COMMON
IMPLICIT NONE

! --------------------------------------------------------
! List of (PRIVATE) routines/functions
! --------------------------------------------------------

       CONTAINS

!=======================================================================

       SUBROUTINE MODELIO_FMT(fmt)
       
       CHARACTER*64, INTENT(OUT)  :: fmt
       CHARACTER*64               :: inifmt, varfmt, frep
       
       ! Construct fortran format string
       
! Output the time and nutrient limitation code
           inifmt='1x, i10.1, 1x, i14.0,'
! Each variable then is a space and a 10 position float with 5 decimal places
           varfmt='1x, f10.5'
! This is the number of repeats (10 variables of nbox dimensions plus pstar and atmpco2)
           write(frep ,'(I4)') 10*nbox+2
! Combine everything together
!          fmt='('//trim(fmt)//trim(frep)//'('//trim(varfmt)//'))'
           write(fmt,'(6A)') '(',trim(inifmt),trim(frep),'(',trim(varfmt),'))'

       RETURN
       END SUBROUTINE MODELIO_FMT     
!=======================================================================
       
       SUBROUTINE MODELIO_HEADER(filename_avg                          &
                                )
       
       CHARACTER*64, INTENT(IN)  :: filename_avg
       INTEGER                   :: iounit

       iounit = 14
       
#if defined(WRITESTOUTFILE)      
! open an output file and write initial values to file
       open(iounit, file=filename_avg, status="unknown")          

! write column header output 
       write(iounit,*)'     t(yr)     Limits  ',                       &
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

! close the output file
       close(iounit)
#endif
       RETURN
       END SUBROUTINE MODELIO_HEADER
!=======================================================================       

       SUBROUTINE MODELIO_OUTPUT(                                      &
                      filename_avg,                                    &
                           outstep,                                    &
                        outstepmax,                                    &
                             timeM,                                    & 
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
                             pco2A,                                    &
                             thout,                                    &
                              sout,                                    &
                              cout,                                    &
                              aout,                                    &
                              pout,                                    &
                              nout,                                    &
                              fout,                                    &
                              lout,                                    &
                            expout,                                    &
                         ocpco2out,                                    &
                              tout,                                    &
                             nlout,                                    &
                             psout,                                    &
                         atpco2out                                     &                
                             )

       REAL(KIND=wp), INTENT(IN), DIMENSION(nbox) ::                   &
                               thetaM, saltM, exportM, pco2M,          &
                               dicM, alkM, po4M, no3M, fetM, ltM                                   
       REAL(KIND=wp), INTENT(IN)                  ::                   &
                               timeM, pco2A, pstarM

       CHARACTER*64, INTENT(IN)          :: filename_avg
       INTEGER,      INTENT(IN)          :: outstep, outstepmax
       INTEGER(kind=ip),      INTENT(IN) :: lim

       REAL(KIND=wp), intent(inout), dimension (outstepmax,nbox) ::    &
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

       REAL(KIND=wp), intent(inout), dimension (outstepmax) ::         &
            tout,                                                      &
            psout,                                                     &
            atpco2out

       INTEGER(KIND=ip), intent(inout), dimension (outstepmax) ::      &
            nlout     
            
! Always fill the arrays
       CALL MODELIO_ARRAY_OUTPUT(                                      &
                           outstep,                                    &
                        outstepmax,                                    &
                             timeM,                                    & 
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
                             pco2A,                                    &
                             thout,                                    &
                              sout,                                    &
                              cout,                                    &
                              aout,                                    &
                              pout,                                    &
                              nout,                                    &
                              fout,                                    &
                              lout,                                    &
                            expout,                                    &
                         ocpco2out,                                    &
                              tout,                                    &
                             nlout,                                    &
                             psout,                                    &
                         atpco2out                                     &                 
                             )

! Only write out to textfile if needed
#if defined(WRITESTOUTFILE)  
       CALL MODELIO_STOUT_OUTPUT(                                      &
                           filename_avg,                               &
                           outstep,                                    &
                        outstepmax,                                    &
                             timeM,                                    & 
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
                             pco2A                                     &
                             )
#endif

#if defined(USEJSONOUT)
! Only write to JSON file at the end of the run, if requested
       if (mod(outstep,outstepmax).eq.zero) then
        call MODELIO_FJSON_OUTPUT(                                     &
                      filename_avg,                                    &
                             thout,                                    &
                              sout,                                    &
                              cout,                                    &
                              aout,                                    &
                              pout,                                    &
                              nout,                                    &
                              fout,                                    &
                              lout,                                    &
                            expout,                                    &
                         ocpco2out,                                    &
                              tout,                                    &
                             nlout,                                    &
                             psout,                                    &
                         atpco2out                                     &             
                             )
       endif ! output time check
#endif
       RETURN
       END SUBROUTINE MODELIO_OUTPUT
!=======================================================================

       SUBROUTINE MODELIO_STOUT_OUTPUT(                                &
                           filename_avg,                               &
                           outstep,                                    &
                        outstepmax,                                    &
                             timeM,                                    & 
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
                             pco2A                                     &
                             )

       REAL(KIND=wp), INTENT(IN), DIMENSION(nbox) ::                   &
                               thetaM, saltM, exportM, pco2M,          &
                               dicM, alkM, po4M, no3M, fetM, ltM                                   
       REAL(KIND=wp), INTENT(IN)                  ::                   &
                               timeM, pco2A, pstarM

       CHARACTER*64, INTENT(IN)          :: filename_avg
       INTEGER,      INTENT(IN)          :: outstep, outstepmax
       INTEGER(kind=ip),      INTENT(IN) :: lim
       
       INTEGER      :: iounit
       CHARACTER*64 :: fmt, avg_fname             

       iounit = 14

       write (avg_fname , '(a,a)') trim(filename_avg) ,'.dat'

! Set up the I/O numerical format for files
       call modelio_fmt(fmt)

#if defined(WRITESTOUTFILE)           
! Write out the file header if necessary
       if ( timeM .eq.zero) then
           call modelio_header(filename_avg                            &
                          )
       endif
       
! Append model state to file  
       open(iounit, file=filename_avg, status="old",                   &
                position="append")

       write(iounit,fmt) int(timeM),                                   & 
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
       close(iounit)                     
#endif
       RETURN
       END SUBROUTINE MODELIO_STOUT_OUTPUT
!=======================================================================

       SUBROUTINE MODELIO_ARRAY_OUTPUT(                                &
                           outstep,                                    &
                        outstepmax,                                    &
                             timeM,                                    & 
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
                             pco2A,                                    &
                             thout,                                    &
                              sout,                                    &
                              cout,                                    &
                              aout,                                    &
                              pout,                                    &
                              nout,                                    &
                              fout,                                    &
                              lout,                                    &
                            expout,                                    &
                         ocpco2out,                                    &
                              tout,                                    &
                             nlout,                                    &
                             psout,                                    &
                         atpco2out                                     &              
                             )

       REAL(KIND=wp), INTENT(IN), DIMENSION(nbox) ::                   &
                               thetaM, saltM, exportM, pco2M,          &
                               dicM, alkM, po4M, no3M, fetM, ltM                                   
       REAL(KIND=wp), INTENT(IN)                  ::                   &
                               timeM, pco2A, pstarM

       INTEGER,               INTENT(IN) :: outstep, outstepmax
       INTEGER(kind=ip),      INTENT(IN) :: lim

       REAL(KIND=wp), intent(inout), dimension (outstepmax,nbox) ::    &
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

       REAL(KIND=wp), intent(inout), dimension (outstepmax) ::         &
            tout,                                                      &
            psout,                                                     &
            atpco2out

       INTEGER(KIND=ip), intent(inout), dimension (outstepmax) ::      &
            nlout     
            
! output to array
       thout     (outstep,1:nbox) = thetaM
       sout      (outstep,1:nbox) = saltM
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
       RETURN
       END SUBROUTINE MODELIO_ARRAY_OUTPUT
!=======================================================================

       SUBROUTINE MODELIO_FJSON_OUTPUT(                                &
                      filename_avg,                                    &
                           outstep,                                    &
                        outstepmax,                                    &
                             thout,                                    &
                              sout,                                    &
                              cout,                                    &
                              aout,                                    &
                              pout,                                    &
                              nout,                                    &
                              fout,                                    &
                              lout,                                    &
                            expout,                                    &
                         ocpco2out,                                    &
                              tout,                                    &
                             nlout,                                    &
                             psout,                                    &
                         atpco2out                                     &            
                             )
! Write output using the json-fortran library 
! 
! JSON-Fortran: A Modern Fortran JSON API
! <https://github.com/jacobwilliams/json-fortran>
! 
! Copyright (c) 2014-2021, Jacob Williams
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without modification,
! are permitted provided that the following conditions are met:
! 
! * Redistributions of source code must retain the above copyright notice, this
!   list of conditions and the following disclaimer.
! 
! * Redistributions in binary form must reproduce the above copyright notice, this
!   list of conditions and the following disclaimer in the documentation and/or
!   other materials provided with the distribution.
! 
! * The names of its contributors may not be used to endorse or promote products
!   derived from this software without specific prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
! ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
! ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! 
! > -----------------------------------------------------------------------------------------
! >
! >	Original FSON License:
! >
! > Copyright (c) 2012 Joseph A. Levin
! >
! > Permission is hereby granted, free of charge, to any person obtaining a copy of this
! > software and associated documentation files (the "Software"), to deal in the Software
! > without restriction, including without limitation the rights to use, copy, modify, merge,
! > publish, distribute, sublicense, and/or sell copies of the Software, and to permit
! > persons to whom the Software is furnished to do so, subject to the following conditions:
! >
! > The above copyright notice and this permission notice shall be included in all copies or
! > substantial portions of the Software.
! >
! > THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
! > INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
! > PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
! > LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
! > OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
! > DEALINGS IN THE SOFTWARE.
! >
! > -----------------------------------------------------------------------------------------
       CHARACTER*64, INTENT(IN)          :: filename_avg

       INTEGER,               intent(in) :: outstep, outstepmax

       REAL(KIND=wp), intent(in), dimension (outstepmax,nbox) ::       &
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

       REAL(KIND=wp), intent(in), dimension (outstepmax) ::            &
            tout,                                                      &
            psout,                                                     &
            atpco2out

       INTEGER(KIND=ip), intent(in), dimension (outstepmax) ::         &
            nlout     
            

#if defined(USEJSONOUT)
       CHARACTER*64 :: avg_fname            

       write (avg_fname , '(a,a)') filename_avg ,'.json'
! Write it out here






#endif ! defined(USEJSONOUT))
       RETURN
       END SUBROUTINE MODELIO_FJSON_OUTPUT
!=======================================================================

       SUBROUTINE MODELIO_FJSON_INPUT(                                 &
            filename, id,                                              &
            maxyears, outputyears, outstepmax,                         &
            dx, dy, dz, depth, latitude,                               &
            K, R, P,                                                   &
            psi, dif,                                                  &
            alpha_yr, gamma, lt_life,                                  &
            dldz, fe_input, wind, fopen,                               &
            thin, sain, cain, alin, phin, niin, fein, liin,            &
            atpco2                                                     &
            )
! Read input using the json-fortran library
! 
! JSON-Fortran: A Modern Fortran JSON API
! <https://github.com/jacobwilliams/json-fortran>
! 
! Copyright (c) 2014-2021, Jacob Williams
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without modification,
! are permitted provided that the following conditions are met:
! 
! * Redistributions of source code must retain the above copyright notice, this
!   list of conditions and the following disclaimer.
! 
! * Redistributions in binary form must reproduce the above copyright notice, this
!   list of conditions and the following disclaimer in the documentation and/or
!   other materials provided with the distribution.
! 
! * The names of its contributors may not be used to endorse or promote products
!   derived from this software without specific prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
! ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
! ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! 
! > -----------------------------------------------------------------------------------------
! >
! >	Original FSON License:
! >
! > Copyright (c) 2012 Joseph A. Levin
! >
! > Permission is hereby granted, free of charge, to any person obtaining a copy of this
! > software and associated documentation files (the "Software"), to deal in the Software
! > without restriction, including without limitation the rights to use, copy, modify, merge,
! > publish, distribute, sublicense, and/or sell copies of the Software, and to permit
! > persons to whom the Software is furnished to do so, subject to the following conditions:
! >
! > The above copyright notice and this permission notice shall be included in all copies or
! > substantial portions of the Software.
! >
! > THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
! > INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
! > PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
! > LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
! > OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
! > DEALINGS IN THE SOFTWARE.
! >
! > -----------------------------------------------------------------------------------------
       CHARACTER*64, INTENT(IN)          :: filename

       INTEGER, INTENT(OUT) :: outstepmax, id, maxyears, outputyears   
       
       REAL(KIND=wp), INTENT(OUT) ::                                   &
            gamma,                                                     &
            lt_life,                                                   &
            alpha_yr,                                                  &
            atpco2,                                                    &
            psi,                                                       &
            dif

! Input arrays (nbox dimensionless)
!       REAL(KIND=wp), DIMENSION(:), ALLOCATABLE, INTENT(OUT) ::        & 
       REAL(KIND=wp), DIMENSION(nbox), INTENT(OUT) ::                  & 
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
            dldz,                                                      &
            wind,                                                      &
            fopen

!       REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) ::      & 
       REAL(KIND=wp), DIMENSION(nbox,nbox), INTENT(OUT) ::             & 
            K,                                                         &
            R,                                                         &
            P

#if defined(USEJSONOUT)
! Local variable definitions
       TYPE(json_file)           :: json
       TYPE(json_value), POINTER :: p_matrix, p_child
       TYPE(json_core)           :: core
       LOGICAL                   :: found
       INTEGER                   :: n_cols, n_rows, var_type
       REAL(KIND=wp)                              :: isca      
       REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: ivec
       REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: imat


! initialize the class
       call json%initialize()
        
! read the file
!json_data => fson_parse("run_microCOSM_4boxvariablelt_pickup.json")
       call json%load_file(filename); if (json%failed()) stop
        
! extract scalar and vector data from the file
! [found can be used to check if the data was really there]
! Iteration ID
       call json%get('niter', isca, found)
       if ( .not. found ) then
          stop 1
       else
          id=isca
       endif

! Length of simulation in years
       call json%get('nyrs', isca, found)
       if ( .not. found ) then
          stop 1
       else
          nyrs=isca
       endif

! Output frequency
       call json%get('tout', isca, found)
       if ( .not. found ) then
          stop 1
       else
          tout=isca
       endif

! Number of outputs, including initial conditions
       call json%get('nout', isca, found)
       if ( .not. found ) then 
          stop 1
       else
          nout=isca
       endif

! Geometry and geography
! Box size in x-direction
       call json%get('dx', ivec, found)
       if ( .not. found ) then 
          stop 1
       else
          dx=ivec
       endif

! Box size in y-direction
       call json%get('dy', ivec, found)
       if ( .not. found ) then 
          stop 1
       else
          dy=ivec
       endif

! Box size in z-direction
       call json%get('dz', ivec, found)
       if ( .not. found ) then 
          stop 1
       else
          dz=ivec
       endif

! Box depth
       call json%get('depth', ivec, found)
       if ( .not. found ) then 
          stop 1
       else
          depth=ivec
       endif

! Box latitude
       call json%get('latitude', ivec, found)
       if ( .not. found ) then 
          stop 1
       else
          latitude=ivec
       endif

! Parameter values
! Overturning circulations strength
       call json%get('psi', isca, found)
       if ( .not. found ) then 
          stop 1
       else
          psi=isca
       endif

! Diffusive mixing strength
       call json%get('dif', isca, found)
       if ( .not. found ) then 
          stop 1
       else
          dif=isca
       endif

! Biological productivity max rate
       call json%get('alphabio', isca, found)
       if ( .not. found ) then 
          stop 1
       else
          alphabio=isca
       endif

! Fraction of biological matter that can be ligand
       call json%get('gamma', isca, found)
       if ( .not. found ) then 
          stop 1
       else
          gamma=isca
       endif

! Ligand degradation timescale
       call json%get('lambda', isca, found)
       if ( .not. found ) then 
          stop 1
       else
          lambda=isca
       endif

! Ligand degradation rate vertical gradient
       call json%get('dlambdadz', ivec, found)
       if ( .not. found ) then 
          stop 1
       else
          dlambdadz=ivec
       endif

! Iron input rate
       call json%get('sourceFe', ivec, found)
       if ( .not. found ) then 
          stop 1
       else
          sourceFe=ivec
       endif

! Wind speed
       call json%get('wind', ivec, found)
       if ( .not. found ) then 
          stop 1
       else
          wind=ivec
       endif

! Open ocean fraction (ice or subsurface)
       call json%get('fopen', ivec, found)
       if ( .not. found ) then 
          stop 1
       else
          fopen=ivec
       endif

! Initial concentrations
! Initial temperatures
       call json%get('theta', ivec, found)
       if ( .not. found ) then 
          stop 1
       else
          theta=ivec
       endif

! Initial salinities
       call json%get('salt', ivec, found)
       if ( .not. found ) then 
          stop 1
       else
          salt=ivec
       endif

! Initial carbon concentrations
       call json%get('carbon', ivec, found)
       if ( .not. found ) then 
          stop 1
       else
          carbon=ivec
       endif

! Initial alkalinities
       call json%get('alkalinity', ivec, found)
       if ( .not. found ) then 
          stop 1
       else
          alkalinity=ivec
       endif

! Initial phosphate concentrations
       call json%get('phosphate', ivec, found)
       if ( .not. found ) then 
          stop 1
       else
          phosphate=ivec
       endif

! Initial nitrate concentrations
       call json%get('nitrate', ivec, found)
       if ( .not. found ) then 
          stop 1
       else
          nitrate=ivec
       endif

! Initial iron concentrations
       call json%get('iron', ivec, found)
       if ( .not. found ) then 
          stop 1
       else
          iron=ivec
       endif

! Initial ligand concentrations
       call json%get('ligand', ivec, found)
       if ( .not. found ) then 
          stop 1
       else
          ligand=ivec
       endif

! Initial atmospheric pCO2
       call json%get('atmpco2', isca, found)
       if ( .not. found ) then 
          stop 1
       else
          atmpco2=isca
       endif
      
! Read circulation matrix (Wish i could subroutine this, but fails compilation 
!   with "Error: Found no matching specific binding for the call to the GENERIC
!   'info'")
        call json%info('Pcir',found,var_type,n_cols)
        if (.not. found) error stop 'error: Pcir not found'
        
        call json%info('Pcir(1)',found,var_type,n_rows)
        if (.not. found) error stop 'error: Pcir(1) not found'
        
!get a pointer to the matrix:
        call json%get('Pcir',p_matrix)
        if (.not. associated(p_matrix)) error stop 'error: Pcir not found'
        
!size the array:
        allocate(imat(n_rows,n_cols))
        
!grab each column of the p_matrix:
! [we need a json_core for this so we can use json_value_get_by_index]
        do i=1,n_cols
            call core%get_child(p_matrix,i,p_child)
            if (.not. associated(p_child)) error stop 'error: column not found'
            call core%get(p_child,ivec) !get the vector (column of the matrix)
            if (.not. allocated(ivec)) error stop 'error: could not get integer column'
            if (size(ivec)/=n_rows) error stop 'error: column is wrong size'
            imat(:,i) = ivec
            deallocate(ivec)
            nullify(p_child)
        end do
        nullify(p_matrix)

        P=imat
        deallocate(imat)

!Read Mixing matrix (Wish i could subroutine this, but fails compilation with
! "Error: Found no matching specific binding for the call to the GENERIC 
! 'info'")
        call json%info('Kmix',found,var_type,n_cols)
        if (.not. found) error stop 'error: Kmix not found'
        
        call json%info('Kmix(1)',found,var_type,n_rows)
        if (.not. found) error stop 'error: Kmix(1) not found'
        
!get a pointer to the matrix:
        call json%get('Kmix',p_matrix)
        if (.not. associated(p_matrix)) error stop 'error: Kmix not found'
        
!size the array:
        allocate(imat(n_rows,n_cols))
        
!grab each column of the p_matrix:
! [we need a json_core for this so we can use json_value_get_by_index]
        do i=1,n_cols
            call core%get_child(p_matrix,i,p_child)
            if (.not. associated(p_child)) error stop 'error: column not found'
            call core%get(p_child,ivec) !get the vector (column of the matrix)
            if (.not. allocated(ivec)) error stop 'error: could not get integer column'
            if (size(ivec)/=n_rows) error stop 'error: column is wrong size'
            imat(:,i) = ivec
            deallocate(ivec)
            nullify(p_child)
        end do
        nullify(p_matrix)

        K=imat
        deallocate(imat)

!Read Rrmemin matrix (Wish i could subroutine this, but fails compilation with
! "Error: Found no matching specific binding for the call to the GENERIC 
! 'info'")
        call json%info('Rremin',found,var_type,n_cols)
        if (.not. found) error stop 'error: Kmix not found'
        
        call json%info('Rremin(1)',found,var_type,n_rows)
        if (.not. found) error stop 'error: Kmix(1) not found'
        
!get a pointer to the matrix:
        call json%get('Rremin',p_matrix)
        if (.not. associated(p_matrix)) error stop 'error: Kmix not found'
        
!size the array:
        allocate(imat(n_rows,n_cols))
        
!grab each column of the p_matrix:
! [we need a json_core for this so we can use json_value_get_by_index]
        do i=1,n_cols
            call core%get_child(p_matrix,i,p_child)
            if (.not. associated(p_child)) error stop 'error: column not found'
            call core%get(p_child,ivec) !get the vector (column of the matrix)
            if (.not. allocated(ivec)) error stop 'error: could not get integer column'
            if (size(ivec)/=n_rows) error stop 'error: column is wrong size'
            Rin(:,i) = ivec
            deallocate(ivec)
            nullify(p_child)
        end do
        nullify(p_matrix)
 
        R=imat
        deallocate(imat)

!! Write out values to check
!         write(6,*)"niter"
!         write(6,'(I0)') niter
!         write(6,*)"nyrs"
!         write(6,'(I0)') nyrs
!         write(6,*)"tout"
!         write(6,'(I0)') tout
!         write(6,*)"nout"
!         write(6,'(I0)') nout
!         
!         write(6,*)"dx"
!         do i=1,nbox
!            write(6,'(f17.4)') dx(i)
!         end do
!         write(6,*)"dy"
!         do i=1,nbox
!            write(6,'(f17.4)') dy(i)
!         end do
!         write(6,*)"dz"
!         do i=1,nbox
!            write(6,'(f17.4)') dz(i)
!         end do
!         write(6,*)"depth"
!         do i=1,nbox
!            write(6,'(f17.4)') depth(i)
!         end do
!         write(6,*)"latitude"
!         do i=1,nbox
!            write(6,'(f17.4)') latitude(i)
!         end do
!         
!         write(6,*)"psi"
!         write(6,'(f17.4)') psi
!         write(6,*)"dif"
!         write(6,'(f17.4)') dif
!         write(6,*)"alphabio"
!         write(6,'(f15.8)') alphabio
!         write(6,*)"gamma"
!         write(6,'(f15.8)') gamma
!         write(6,*)"lambda"
!         write(6,'(f17.4)') lambda
!         
!         write(6,*)"dlambdadz"
!         do i=1,nbox
!            write(6,'(f15.8)') dlambdadz(i)
!         end do
!         write(6,*)"sourceFe"
!         do i=1,nbox
!            write(6,'(f15.8)') sourceFe(i)
!         end do
!         write(6,*)"wind"
!         do i=1,nbox
!            write(6,'(f15.8)') wind(i)
!         end do
!         write(6,*)"fopen"
!         do i=1,nbox
!            write(6,'(f15.8)') fopen(i)
!         end do
!         
!         write(6,*)"theta"
!         do i=1,nbox
!            write(6,'(f15.8)') theta(i)
!         end do
!         write(6,*)"salt"
!         do i=1,nbox
!            write(6,'(f15.8)') salt(i)
!         end do
!         write(6,*)"carbon"
!         do i=1,nbox
!            write(6,'(f15.8)') carbon(i)
!         end do
!         write(6,*)"nitrate"
!         do i=1,nbox
!            write(6,'(f15.8)') nitrate(i)
!         end do
!         write(6,*)"phosphate"
!         do i=1,nbox
!            write(6,'(f15.8)') phosphate(i)
!         end do
!         write(6,*)"iron"
!         do i=1,nbox
!            write(6,'(f15.8)') iron(i)
!         end do
!         write(6,*)"ligand"
!         do i=1,nbox
!            write(6,'(f15.8)') ligand(i)
!         end do
!         
!         write(6,*)"P matrix"
!         do i=1,nbox
!             write(6,400) P(i,:)
!         end do
!         write(6,*)"K matrix"
!         do i=1,nbox
!             write(6,400) K(i,:)
!         end do
!         write(6,*)"R matrix"
!         do i=1,nbox
!             write(6,400) R(i,:)
!         end do
!         
!         write(6,*)"atmpco2"
!         write(6,'(f15.8)') atmpco2
        
! clean up
        call json%destroy()
#endif
       RETURN
       END SUBROUTINE MODELIO_FJSON_INPUT     
!=======================================================================

      END MODULE MOD_MODELIO