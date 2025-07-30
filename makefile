VERSION          = 3.0

EXE =  fei_${VERSION}

SOLID_DIR    = simplex
#MKLDIR   = /opt/intel/mkl/10.1.3.027

OBJ= AMODULES.o			\
	MAIN.o                          \
	CORE_GRID.o                     \
	BOUNDARY_MARKER_INIT.o		\
	GCM_MEMORY_ALLOCATE.o		\
	CORE_MEMORY_ALLOCATE.o          \
	BOUNDARY_SET_ARCLENGTH_NORM.o   \
	BOUNDARY_IBLANK_FAST.o          \
	GCM_SET_BOUNDARY.o              \
	BOUNDARY_MARKER_VEL.o           \
	BOUNDARY_MOVE.o                 \
	CORE_TIME_STEP.o                \
	CORE_FRAC_STEP_ROUTINES.o       \
	CORE_SET_SOLVERS.o              \
	CORE_SOLVE_AD.o                 \
	CORE_POISSON.o                  \
	CORE_TDMA.o                     \
	GCM_BODYINTERCEPT_VALUES.o      \
	CORE_SET_OUTER_BC.o             \
	GCM_GHOSTCELL_UPDATE.o          \
        UTIL_GENERAL.o                  \
	UTIL_DRAG_LIFT.o                \
	UTIL_PROBE.o                    \
	CORE_PARALLEL.o                 \
	CORE_MG.o                       \
	FSI_INTERFACE.o


#OBJ  =  $(SRC:.f=.o)


FC       = mpiifx
FFLAGS   = -O3 -qopenmp -fpp -xCORE-AVX2 -axCORE-AVX512
LDFLAGS  = -qmkl=parallel

#MKLLIBS      = -L$(MKLDIR)/lib/em64t  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core #-lguide
#MKLLIBS      = -L/usr/local/lib64 -llapack -lblas

#MKLLIBS      = -mkl=sequential #-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lguide
MKLLIBS      = -qmkl=parallel

GEN_LIBS     =  $(MKLLIBS) 

LDFLAGS      = $(LDFLAGS2) $(GEN_LIBS) 

SOLID_OBJS   = $(SOLID_DIR)/*.o


all: $(EXE)

%.o: %.F90
	$(FC) $(FFLAGS) -c $< -o $@

solid:  
	(cd $(SOLID_DIR); make)

${EXE}: $(OBJ) solid
	$(FC) $(OBJ) $(SOLID_OBJS) $(LDFLAGS) -o ${EXE}

clean:
	@echo "cleaning ..."
	rm -f $(OBJ) *.mod *~
	#(cd $(SOLID_DIR); make clean)
