OS = $(shell uname)
RM = rm -f
AR = ar
ARFLAGS = rv
FC = gfortran
# Debug options for gfortran
FFLAGS = -I$(MODDIR) -I/usr/local/include -J$(MODDIR) -fbounds-check -O0 -g3 -p #-fopenmp
# # Release options for gfortran
# FFLAGS = -I$(MODDIR) -I/usr/local/include -J$(MODDIR) -O3
# # Profiling options for gfortran
# FFLAGS = -I$(MODDIR) -I/usr/local/include -J$(MODDIR) -O3 -p

# # Debug options for ifort
# FFLAGS = -I/usr/local/include -module $(MODDIR) -O0 -g3 -p -check bounds
# # Release options for ifort
# FFLAGS = -I/usr/local/include -module $(MODDIR) -fast #-fast-transcendentals
# # Profiling options for ifort
# FFLAGS = -I/usr/local/include -module $(MODDIR) -fast -p

LDFLAGS = -L. -llapack -lStreamer -lStreamerTest -lminpack #-lcgns

MODDIR = ./bindings
LIBSRCDIR = ./f90src
LIBTESTSRCDIR = ./f90src/test
MINPACKSRCDIR = ./minpack
DATDIR = ./dat
EXEDIR = ./exe
DRVDIR = .
LIBSOURCES = GeneralUtilities.f90 \
	Riemann.f90 \
	TimeAdvancementStuff.f90 \
	Godunov.f90 \
	ODE_Solvers.f90 \
	BoundaryConditionsStuff.f90 \
	FortranNormalVectors.f90 \
	cgns_interface.f90
all: lib
.POSIX:
lib: libStreamer.a
#	$(RM) $(LIBSRCDIR)/*.o
#	@echo Free object files have been removed. object files can be found in libStreamer.a
#drivers = $(
libStreamer.a: $(LIBSRCDIR)/GeneralUtilities.o \
	$(LIBSRCDIR)/Riemann.o \
	$(LIBSRCDIR)/TimeAdvancementStuff.o \
	$(LIBSRCDIR)/Godunov.o \
	$(LIBSRCDIR)/BoundaryConditionsStuff.o \
	$(LIBSRCDIR)/FortranNormalVectors.o \
#	$(LIBSRCDIR)/cgns_interface.o
#       End dependencies
	$(AR) $(ARFLAGS) $@ $?
$(LIBSRCDIR)/Godunov.o: $(LIBSRCDIR)/GeneralUtilities.o \
	$(LIBSRCDIR)/Riemann.o
$(LIBSRCDIR)/ODE_Solvers.o: $(LIBSRCDIR)/Source_functions.o
$(LIBSRCDIR)/Source_functions.o: $(LIBSRCDIR)/GeneralUtilities.o
$(LIBSRCDIR)/BoundaryConditionsStuff.o: $(LIBSRCDIR)/GeneralUtilities.o
check: lib libStreamerTest.a libminpack.a\
	$(LIBTESTSRCDIR)/test_program.f90 \
	$(DRVDIR)/Godunov_driver.f90 \
# 	End dependencies
	$(FC) $(FFLAGS) $(LDFLAGS) $(LIBTESTSRCDIR)/test_program.f90 \
	$(DRVDIR)/Godunov_driver.o  -o $(LIBTESTSRCDIR)/f90tests
	$(LIBTESTSRCDIR)/f90tests $(DATDIR)
libStreamerTest.a: $(LIBTESTSRCDIR)/GeneralUtilities_tester.o \
	$(LIBTESTSRCDIR)/Godunov_tester.o \
	$(LIBTESTSRCDIR)/Riemann_tester.o 
#	$(LIBTESTSRCDIR)/ODE_Solvers_tester.o
	$(AR) $(ARFLAGS) $@ $?
libminpack.a: $(MINPACKSRCDIR)/lmder1.o \
	$(MINPACKSRCDIR)/dpmpar.o \
	$(MINPACKSRCDIR)/lmder.o \
	$(MINPACKSRCDIR)/lmpar.o \
	$(MINPACKSRCDIR)/qrsolv.o \
	$(MINPACKSRCDIR)/enorm.o \
	$(MINPACKSRCDIR)/qrfac.o
	$(AR) $(ARFLAGS) $@ $?
$(LIBTESTSRCDIR)/Godunov_tester.o: $(DRVDIR)/Godunov_driver.o
#$(EXEDIR)/test_program: $(SRCDIR)/test_program.f90 lib
#	$(FC) $(FFLAGS) $(LDFLAGS) $< -o $@
%.o: %.f90
	$(FC) -c $(FFLAGS) $< -o $@
clean: 
	$(RM) $(LIBSRCDIR)/*.o $(LIBTESTSRCDIR)/*.o $(DRVDIR)/Godunov_driver.o \
	$(MODDIR)/*.mod libStreamer.a libStreamerTest.a *.so *.a *.mod *.o \
	minpack/*.o
