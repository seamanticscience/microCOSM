.SUFFIXES:
.SUFFIXES: .o .mod .f90 .F90 .f

RM = /bin/rm -f
CP = /bin/cp -f

FC = gfortran
PC = f2py
override FFLAGS += -g
OPTIONDEFS = -DWRITEOUTFILE

# If using JSON output, compile with json-fortran library
ifeq ($(findstring USEJSONOUT,$(OPTIONDEFS)),USEJSONOUT)
JSON_DYL = -ljsonfortran -lgfortran
JSON_DIR = /Users/jml1/GitHub/microCOSM/json-fortran/jsonfortran-gnu-8.3.0/lib/
JSON_LIB = -L$(JSON_DIR)
JSON_INC = -I$(JSON_DIR)
FFLAGS += -Wl,-rpath,/opt/local/lib/gcc-devel -Wl,-rpath,$(JSON_DIR)
else
JSON_DYL =
JSON_LIB =
JSON_INC =
endif

MODULE_OBJ = mod_precision.o mod_dimensions.o         \
             mod_common.o    mod_chemconst.o          \
             mod_chemspeciation.o mod_phsolvers.o     \
             mod_carbonchem.o mod_modelio.o           \
             mod_modelmain.o
             
MODEL_OBJ  = microcosm_model.o

# C preprocessing 
CPPCMD = cat $< | $(FC) $(OPTIONDEFS) -cpp -P -E $(FFLAGS)

model: $(MODULE_OBJ) $(MODEL_OBJ)
	$(FC) $(JSON_LIB) $(JSON_INC) $(OPTIONDEFS) $(FFLAGS) $(MODULE_OBJ) $(MODEL_OBJ) -o microCOSM $(JSON_DYL)

# Requires seperate pre-processing for some reason
pymodel: $(MODULE_OBJ:.o=.f) $(MODEL_OBJ:.o=.f)
	$(PC) $(OPTIONDEFS) --opt=$(FFLAGS) --f90flags=-ffree-form -m microCOSM -c $(MODULE_OBJ:.o=.f) $(MODEL_OBJ:.o=.f)

#%.o: %.F90
#	$(FC) $(OPTIONDEFS) -c $(FFLAGS) $< -o $@
	
%.f: %.F90
	$(CPPCMD) -o $@ -
%.o: %.f
	$(FC) $(JSON_LIB) $(JSON_INC) $(OPTIONDEFS) -c $(FFLAGS) -ffree-form -o $@ $< $(JSON_DYL)

.PHONY : clean
clean:
	$(RM) $(MODULE_OBJ)         $(MODEL_OBJ)
	$(RM) $(MODULE_OBJ:.o=.mod) $(MODEL_OBJ:.o=.mod)
	$(RM) $(MODULE_OBJ:.o=.f)   $(MODEL_OBJ:.o=.f)
