.SUFFIXES:
.SUFFIXES: .o .mod .f90 .F90 .f

RM = /bin/rm -f
CP = /bin/cp -f

FC = gfortran
PC = f2py
override FFLAGS += -g
OPTIONDEFS = -DWRITEOUTFILE

DNAD_OBJ = dnad/dnad.o

MODULE_OBJ = mod_precision.o mod_dimensions.o         \
             mod_common.o    mod_chemconst.o          \
             mod_chemspeciation.o mod_phsolvers.o     \
             mod_carbonchem.o mod_modelmain.o
             
MODEL_OBJ  = microcosm_model.o

# C preprocessing 
CPPCMD = cat $< | $(FC) $(OPTIONDEFS) -cpp -P -E $(FFLAGS)

model: $(MODULE_OBJ) $(MODEL_OBJ)
	$(FC) $(OPTIONDEFS) $(FFLAGS) $(MODULE_OBJ) $(MODEL_OBJ) -o microCOSM 

# Requires seperate pre-processing for some reason
pymodel: $(MODULE_OBJ:.o=.f) $(MODEL_OBJ:.o=.f)
	$(PC) $(OPTIONDEFS) --opt=$(FFLAGS) --f90flags=-ffree-form -m microCOSM -c $(MODULE_OBJ:.o=.f) $(MODEL_OBJ:.o=.f)

#%.o: %.F90
#	$(FC) $(OPTIONDEFS) -c $(FFLAGS) $< -o $@
	
%.f: %.F90
	$(CPPCMD) -o $@ -
%.o: %.f
	$(FC) $(OPTIONDEFS) -c $(FFLAGS) -ffree-form  $< -o $@

.PHONY : clean
clean:
	$(RM) $(MODULE_OBJ)         $(DNAD_OBJ)         $(MODEL_OBJ)
	$(RM) $(MODULE_OBJ:.o=.mod) $(DNAD_OBJ:.o=.mod) $(MODEL_OBJ:.o=.mod)
	$(RM) $(MODULE_OBJ:.o=.f)   $(DNAD_OBJ:.o=.f)   $(MODEL_OBJ:.o=.f)
