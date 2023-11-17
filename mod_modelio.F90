! -*- f90 -*-
       MODULE MOD_MODELIO
#if defined(USEJSONOUT) || defined(USEJSONIN)
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
      
#if defined(WRITESTOUTFILE)   
       CHARACTER*64 :: avg_fname             

       write (avg_fname , '(a,a)') trim(filename_avg) ,'.dat'

       iounit = 14

! open an output file and write initial values to file
       open(iounit, file=avg_fname, status="unknown")          

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
        call MODELIO_JSON_OUTPUT(                                     &
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
       open(iounit, file=avg_fname, status="old",                      &
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

       SUBROUTINE MODELIO_JSON_OUTPUT(                                &
                      filename_avg,                                    &
!                           outstep,                                    &
!                        outstepmax,                                    &
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

!       INTEGER,               intent(in) :: outstep, outstepmax

       REAL(KIND=wp), intent(in), dimension (:,:) ::                   &
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

       REAL(KIND=wp), intent(in), dimension (:) ::                     &
            tout,                                                      &
            psout,                                                     &
            atpco2out

       INTEGER(KIND=ip), intent(in), dimension (:) ::                  &
            nlout     
            

#if defined(USEJSONOUT)
       TYPE(json_core)                           :: json 
       TYPE(json_value),POINTER                  :: root
       LOGICAL                                   :: is_valid
       INTEGER                                   :: iunit 
       CHARACTER(KIND=json_CK,LEN=:),ALLOCATABLE :: error_msg
       CHARACTER*64                              :: avg_fname            

! Write filename to a variable
       write (avg_fname , '(a,a)') trim(filename_avg) ,'.json'
       
! Initialize the JSON object       
       call json%initialize(verbose=.true.,compress_vectors=.true.); &
         if (json%failed()) stop
! Create the root ouput object to store all the data
       call json%create_object(root,avg_fname); if (json%failed()) stop
!input structure:
! Would be great to write the parameters with the output, but struggling to
! find a way to write out 2-d matrices for P, K, and R.
! call json%create_object(inp,'inputs'); if (json%failed()) stop
! call json%add(p, inp); if (json%failed()) stop
!
! call json%add(inp, 't0', 0.1_wp); if (json%failed()) stop
! call json%add(inp, 'tf', 1.1_wp); if (json%failed()) stop
! call json%add(inp, 'x0', 9999.000_wp); if (json%failed()) stop
! nullify(inp)
!model output structure:
! call json%create_object(outp,'output'); if (json%failed()) stop
! call json%add(root, outp); if (json%failed()) stop

       ! Can just output directly with vectors
       call json%add(root, 'time'   , int(tout) ); if (json%failed()) stop
       call json%add(root, 'lim'    , int(nlout)); if (json%failed()) stop
       call json%add(root, 'pstar'  , psout     ); if (json%failed()) stop
       call json%add(root, 'atmpco2', atpco2out ); if (json%failed()) stop
       ! Easier to output arrays with a function
       call modelio_json_add_arr(json, root, 'theta'     , thout )
       call modelio_json_add_arr(json, root, 'salt'      , sout )
       call modelio_json_add_arr(json, root, 'carbon'    , cout )
       call modelio_json_add_arr(json, root, 'alkalinity', aout )
       call modelio_json_add_arr(json, root, 'phosphate' , pout )
       call modelio_json_add_arr(json, root, 'nitrate'   , nout )
       call modelio_json_add_arr(json, root, 'iron'      , fout )
       call modelio_json_add_arr(json, root, 'ligand'    , lout )
       call modelio_json_add_arr(json, root, 'export'    , expout )
       call modelio_json_add_arr(json, root, 'pco2'      , ocpco2out )
       ! Clean up a bit
       !nullify(outp)
       
       !validate the JSON structure
       call json%validate(root, is_valid, error_msg)
       if (.not. is_valid) then
           stop &
           'Error: root is not a valid JSON linked list: '
           write(*,*) error_msg
       end if
       
       ! Write out to file
       open(newunit=iunit, file=avg_fname, status='REPLACE')
       call json%print(root,iunit); if (json%failed()) stop
       close(iunit)
       
       ! Clean up the rest
       call json%destroy(root); if (json%failed()) stop
       
       RETURN
#endif
       END SUBROUTINE MODELIO_JSON_OUTPUT
!=======================================================================
#if defined(USEJSONOUT)
       SUBROUTINE MODELIO_JSON_ADD_ARR(json, me, variable, boxdata)
! Write array output using the json-fortran library
       IMPLICIT NONE
       TYPE(json_core),INTENT(INOUT) :: json
       TYPE(json_value),POINTER :: me, var
       CHARACTER(LEN=*),INTENT(IN) :: variable
       REAL(KIND=WP),DIMENSION(:,:),INTENT(IN) :: boxdata
       INTEGER :: i
       
       !initialize:
       nullify(var)
       
       !create the variable array before data can be added
       call json%create_array(var, trim(variable)); if (json%failed()) stop
       
       ! Add to the output object
       call json%add(me, var); if (json%failed()) stop
       
       ! Cycle through the data, adding each box/vector to the variable array
       do i=1,nbox
          call json%add(var, '', boxdata(:,i)); if (json%failed()) stop
       enddo
       
       !cleanup:
       nullify(var)
       
       RETURN
       END SUBROUTINE MODELIO_JSON_ADD_ARR
#endif
!=======================================================================

    SUBROUTINE MODELIO_JSON_INPUT(                                     &
            filename, id,                                              &
            maxyears, outputyears, outstepmax,                         &
            dx, dy, dz, depth, latitude,                               &
            Kin, Rin, Pin,                                             &
            psi_in, dif_in,                                            &
            alpha_yr, gamma_in, lt_lifein,                             &
            dldz_in, fe_input, wind_in, fopen_in,                      &
            thin, sain, cain, alin, phin, niin, fein, liin,            &
            atpco2in                                                   &
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

       INTEGER, INTENT(OUT) :: outstepmax, id  
       
       REAL(KIND=wp), INTENT(OUT) ::                                   &
            maxyears,                                                  &
            outputyears,                                               &
            gamma_in,                                                  &
            lt_lifein,                                                 &
            alpha_yr,                                                  &
            atpco2in,                                                  &
            psi_in,                                                    &
            dif_in

! Input arrays (nbox dimensionless)
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
            dldz_in,                                                   &
            wind_in,                                                   &
            fopen_in

       REAL(KIND=wp), DIMENSION(nbox,nbox), INTENT(OUT) ::             & 
            Kin,                                                       &
            Rin,                                                       &
            Pin

#if defined(USEJSONIN)
! Local variable definitions
       TYPE(json_file)           :: json
       TYPE(json_value), POINTER :: p_matrix, p_child
       TYPE(json_core)           :: core
       LOGICAL                   :: found
       INTEGER                   :: i, n_cols, n_rows, var_type
       REAL(KIND=wp)                              :: isca      
       REAL(KIND=wp), DIMENSION(:),   ALLOCATABLE :: ivec
       REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: imat


! initialize the class
       call json%initialize(); if (json%failed()) stop         &
            'error initializing json-fortran'
        
! read the file
!json_data => fson_parse("run_microCOSM_4boxvariablelt_pickup.json")
       call json%load_file(filename); if (json%failed()) stop         &
            'error: microcosm_pickup.json not found'
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
          maxyears=isca
       endif

! Output frequency
       call json%get('tout', isca, found)
       if ( .not. found ) then
          stop 1
       else
          outputyears=isca
       endif

! Number of outputs, including initial conditions
       call json%get('nout', isca, found)
       if ( .not. found ) then 
          stop 1
       else
          outstepmax=isca
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
          psi_in=isca
       endif

! Diffusive mixing strength
       call json%get('dif', isca, found)
       if ( .not. found ) then 
          stop 1
       else
          dif_in=isca
       endif

! Biological productivity max rate
       call json%get('alphabio', isca, found)
       if ( .not. found ) then 
          stop 1
       else
          alpha_yr=isca
       endif

! Fraction of biological matter that can be ligand
       call json%get('gamma', isca, found)
       if ( .not. found ) then 
          stop 1
       else
          gamma_in=isca
       endif

! Ligand degradation timescale
       call json%get('lambda', isca, found)
       if ( .not. found ) then 
          stop 1
       else
          lt_lifein=isca
       endif

!       gaovla_opt = ((gamma_in/106._wp))/lt_lifein
       
! Ligand degradation rate vertical gradient
       call json%get('dlambdadz', ivec, found)
       if ( .not. found ) then 
          stop 1
       else
          dldz_in=ivec
       endif

! Iron input rate
       call json%get('sourceFe', ivec, found)
       if ( .not. found ) then 
          stop 1
       else
          fe_input=ivec
       endif

! Wind speed
       call json%get('wind', ivec, found)
       if ( .not. found ) then 
          stop 1
       else
          wind_in=ivec
       endif

! Open ocean fraction (ice or subsurface)
       call json%get('fopen', ivec, found)
       if ( .not. found ) then 
          stop 1
       else
          fopen_in=ivec
       endif

! Initial concentrations
! Initial temperatures
       call json%get('theta', ivec, found)
       if ( .not. found ) then 
          stop 1
       else
          thin=ivec
       endif

! Initial salinities
       call json%get('salt', ivec, found)
       if ( .not. found ) then 
          stop 1
       else
          sain=ivec
       endif

! Initial carbon concentrations
       call json%get('carbon', ivec, found)
       if ( .not. found ) then 
          stop 1
       else
          cain=ivec
       endif

! Initial alkalinities
       call json%get('alkalinity', ivec, found)
       if ( .not. found ) then 
          stop 1
       else
          alin=ivec
       endif

! Initial phosphate concentrations
       call json%get('phosphate', ivec, found)
       if ( .not. found ) then 
          stop 1
       else
          phin=ivec
       endif

! Initial nitrate concentrations
       call json%get('nitrate', ivec, found)
       if ( .not. found ) then 
          stop 1
       else
          niin=ivec
       endif

! Initial iron concentrations
       call json%get('iron', ivec, found)
       if ( .not. found ) then 
          stop 1
       else
          fein=ivec
       endif

! Initial ligand concentrations
       call json%get('ligand', ivec, found)
       if ( .not. found ) then 
          stop 1
       else
          liin=ivec
       endif

! Initial atmospheric pCO2
       call json%get('atmpco2', isca, found)
       if ( .not. found ) then 
          stop 1
       else
          atpco2in=isca
       endif
      
! Read circulation matrix (Wish i could subroutine this, but fails compilation 
!   with  
! "Error: Found no matching specific binding for the call to the GENERIC 'info'")
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

        Pin=imat
        deallocate(imat)

!Read Mixing matrix (Wish i could subroutine this, but fails compilation with
! "Error: Found no matching specific binding for the call to the GENERIC 'info'")

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

        Kin=imat
        deallocate(imat)

!Read Rrmemin matrix (Wish i could subroutine this, but fails compilation with
! "Error: Found no matching specific binding for the call to the GENERIC 'info'")
        call json%info('Rremin',found,var_type,n_cols)
        if (.not. found) error stop 'error: Rremin not found'
        
        call json%info('Rremin(1)',found,var_type,n_rows)
        if (.not. found) error stop 'error: Rremin(1) not found'
        
!get a pointer to the matrix:
        call json%get('Rremin',p_matrix)
        if (.not. associated(p_matrix)) error stop 'error: Rremin not found'
        
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
 
        Rin=imat
        deallocate(imat)

! clean up
        call json%destroy()
#endif
       RETURN
       END SUBROUTINE MODELIO_JSON_INPUT     
!=======================================================================

      END MODULE MOD_MODELIO