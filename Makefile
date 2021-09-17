
.SUFFIXES:

.SUFFIXES: .o .f90 .F90

RM = /bin/rm -f

FC = gfortran
PC = f2py
FFLAGS = '-O3'
OPTIONDEFS = -DFIXATMPCO2

MODULE_OBJS_MICROCOSM = mod_precision.o mod_chemconst.o mod_chemspeciation.o mod_phsolvers.o mod_carbonchem.o mod_modelmain.o
MODULE_F90_MICROCOSM = mod_precision.f90 mod_chemconst.f90 mod_chemspeciation.f90 mod_phsolvers.F90 mod_carbonchem.f90 mod_modelmain.F90

default: $(MODULE_OBJS)

model: $(MODULE_OBJS_MICROCOSM) microcosm_model.o
	$(FC) $(OPTIONDEFS) $(FFLAGS) $(MODULE_OBJS_MICROCOSM) microcosm_model.o -o microCOSM 

pymodel: $(MODULE_OBJS_MICROCOSM) microcosm_model.o
	$(PC) $(OPTIONDEFS) --opt=$(FFLAGS) -m microCOSM -c $(MODULE_F90_MICROCOSM) microcosm_model.F90
.f90.o:;
	$(FC) $(OPTIONDEFS) -c $(FFLAGS) $*.f90 -o $*.o

.F90.o:;
	$(FC) $(OPTIONDEFS) -c $(FFLAGS) $*.F90 -o $*.o

clean:
	$(RM) *.o *.mod
