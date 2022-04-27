# Make changes here according to the needs of your compiler and MPI distribution
FLIBS           =
#FLIBS		= -L$(HOME)/mpich2/lib #This is when I compiled an MPI version myself. Try not to do that.
FFLAGS		= -fast
FXX      = mpif90    # fortran compiler to use
FLD  = $(FXX)  # Linker to use for fortran code.
#PREPROCESSOR = cpp
PREPROCESSOR = fpp
FLDFLAGS = 
#--------
# Like the previous versions, b6 can be compiled for different 
# magnetic and gravitational models, but so far l7 is the only one 
# implemented.
KETCHUPL7OBJECTS = SpecificTypes.o fvlasov.o Advection.o Compute.o \
	Dump.o Initial.o ModelL7.o
REGENERATEL7OBJECTS = SpecificTypes.o Regeneration.o Compute.o \
	Dump.o Initial.o ModelL7.o
MAXF0UPDATEOBJECTS = SpecificTypes.o Maxf0update.o Dump.o Initial.o


#--------
# We pass the fortran program files ending in .fpp through 
# the C preprocessor. Defining _SAFETY_ enables a few tests 
# that otherwise are disabled. Defining _DEBUG_ makes the 
# program print messages every now and then.
#--------

#--------
KETCHUPL7TARGET = l7/ketchup
REGENERATEL7TARGET = l7/regenerate_ketchup
MAXF0UPDATETARGET = maxf0update
#--------
#DEFINES = -D_SAFETY_ -D_DEBUG_
DEFINES = -D_SAFETY_ #-D_TIMING_
#DEFINES = -D_DEBUG_
#DEFINES = #-D_TIMING_
#--------


# make rules

all : 	$(KETCHUPL7TARGET) $(REGENERATEL7TARGET) $(MAXF0UPDATETARGET)

l7 : $(KETCHUPL7TARGET) $(REGENERATEL7TARGET) $(MAXF0UPDATETARGET)


$(KETCHUPL7TARGET): $(KETCHUPL7OBJECTS)
	$(FLD) $(KETCHUPL7OBJECTS) $(FLDFLAGS) -o $(KETCHUPL7TARGET) $(FLIBS)

$(MAXF0UPDATETARGET): $(MAXF0UPDATEOBJECTS)
	$(FLD) $(MAXF0UPDATEOBJECTS) $(FLDFLAGS) -o $(MAXF0UPDATETARGET) \
	$(FLIBS)

$(REGENERATEL7TARGET): $(REGENERATEL7OBJECTS)
	$(FLD) $(REGENERATEL7OBJECTS) $(FLDFLAGS) \
	-o $(REGENERATEL7TARGET) $(FLIBS)

clean : 
	\rm -rf *.o *.mod *~ 


# Rule used to generate object files from fortran source code:
%.o: %.f
	$(FXX) $(FFLAGS) -c $<

%.o: %.for
	$(FXX) $(FFLAGS) -c $<

%.o: %.f90
	$(FXX) $(FFLAGS) -c $<

%.o: %.F
	$(FXX) $(FFLAGS) -c $<

%.f90: %.fpp
	$(PREPROCESSOR) $(DEFINES) -C -P -E $< -o $@
#	$(PREPROCESSOR) $(DEFINES) -C -P -E $< -o $*.f95
