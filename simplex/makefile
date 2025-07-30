FC       = mpiifx
FFLAGS   = -O3 -qopenmp -fpp -xCORE-AVX2 -axCORE-AVX512
LDFLAGS  = -qmkl=parallel

EXE = nonstad

all: $(EXE)

OBJS = CMODULES.o   \
       INTERCEPT_KERNELS_UNSTRUC.o   \
       CALCULATE_TRIANGLE_NORM.o  \
       main320b.o  assm310b.o dtin300b.o \
       eign310b.o  nonb324b.o nond320b.o \
       nond320b_1.o  nond320b_elas_ye_OMP.o nonl310b.o \
       nonr310b.o serv302.o solv404.o \
       stat300.o util200.o util201.o

# %.o: %.for
# 	$(F90) $(FFLAGS) $(LFLAGS) -c $< -o $@

%.o: %.for
	$(FC) $(FFLAGS) -c $< -o $@

#main402.obj: main402.f90
#	$(F90) $(FFLAGS) -c main402.f90 -I../ -o main402.obj 

$(EXE): $(OBJS)
#	${F90} -o ${EXE} $(OBJS)  ${LFLAGS}

clean:
	rm -f *.o #${EXE}
	rm -f stadyn.* *.log *.cfg
	rm -f *~


