#FC	= gfortran-mp-4.8 -O2 -g
FC	= gfortran
MC      = mpif90

LIB_DIR = /opt/local/lib
MOD_DIR = /opt/local/include

HDF5_DIR = /opt/local/include

#lapack + hdf5 flags
LDFLAGS = -framework accelerate -lcaf_mpi -lhdf5_fortran -lhdf5

CAF_FLAGS = -fcoarray=lib

# double precision
FFLAGS_OPT2 = -ggdb -g -O3 -fdefault-real-8 -fdefault-double-8\
	-ffree-line-length-none -Wuninitialized

# quadruple precision
FFLAGS_OPT4 = -ggdb -g -O3 -freal-4-real-16 -freal-8-real-16\
	-ffree-line-length-none -Wuninitialized

# double precision
FFLAGS_DEBUG2 = -ggdb  -g -fdefault-real-8 -fdefault-double-8\
#	-ffree-line-length-none -Wuninitialized

# quadruple precision
FFLAGS_DEBUG4 = -ggdb  -g -freal-4-real-16 -freal-8-real-16\
	-ffree-line-length-none -Wuninitialized

FFLAGS_DOUB = $(FFLAGS_OPT2)
FFLAGS_QUAD = $(FFLAGS_OPT4)

EXE_FILE = slugEuler
OBJS  =	driver_euler1d.o \
	read_initFile.o \
	read_pars.o \
	sim_data.o  \
	gp_data.o \
	gp_WENOinit.o \
	gp_MDinit.o \
	gp_Fluxinit.o \
	gp_MakeStencil.o \
	gp_eigens.o \
	gpM_eigens.o \
	block_data.o \
	block_init.o \
	block_finalize.o \
	sim_init.o \
	sim_initBlock.o \
	io.o \
	bc.o \
	bc_corners.o \
	DMR.o \
	eos.o\
	grid_finalize.o \
	grid_init.o \
	grid_data.o \
	primconsflux.o \
	cfl.o \
	soln_ReconEvolveAvg.o \
	soln_FOG.o \
	soln_RK2.o \
	soln_RK3.o \
	soln_RK4.o \
	reconstruction.o \
	soln_reconstruct.o \
	sim_interfaces.o \
	soln_getFlux.o \
	soln_intFlux.o \
	soln_cntrFlux.o \
	soln_gpFlux.o \
	hll.o \
	hllc.o \
	roe.o \
	averageState.o \
	soln_update.o \
	soln_PLM.o \
	slopeLimiter.o \
	soln_PPM.o \
	soln_WENO.o \
	GP.o \
	WENO.o \
	soln_GP.o \
	soln_gpWENO.o \
	linalg.o \
	eigensystem.o 
	#averageState.o \

########################################################################################
#COMPILING AND LINKING USING GENERIC SUFFIX RULE FOR F90

$(EXE_FILE) : $(OBJS)
	@$(MC) -L $(LIB_DIR) -I $(MOD_DIR) $(FFLAGS_DOUB) $(OBJS) -o $(EXE_FILE) $(LDFLAGS)
	@echo "code is now linking..."

#SOME GP ROUTINES NEED QUAD-PRECISION
GP.o: %.o : %.F90
	$(FC) $(FFLAGS_QUAD) -c $<

gp_WENOinit.o: %.o : %.F90
	$(FC) $(FFLAGS_QUAD) -c $<

gp_MDinit.o: %.o : %.F90
	$(FC) $(FFLAGS_QUAD) -c $<

gp_Fluxinit.o: %.o : %.F90
	$(FC) $(FFLAGS_QUAD) -c $<

linalg.o: %.o : %.F90
	$(FC) $(FFLAGS_QUAD) -c $<

#LET'S APPLY GENERIC SUFFIX RULE HERE FOR FORTRAN 90
.SUFFIXES : 
.SUFFIXES : .F90 .o

.F90.o:
	$(FC) $(FFLAGS_DOUB) -I $(HDF5_DIR)  -c $< $(CAF_FLAGS)

#######################################################################################
#SOME USEFUL COMMANDS
clean:
	@rm -f *.o *.mod *~ slugEuler1d

#######################################################################################
#LET'S DEFINE SOME MODULE DEPENDENCIES!
driver_euler1d.o: sim_data.o grid_data.o io.o  eos.o gp_data.o bc.o

read_pars.o     : sim_data.o grid_data.o block_data.o gp_data.o

eos.o		: grid_data.o sim_data.o
cfl.o           : grid_data.o sim_data.o eigensystem.o

grid_init.o	: grid_data.o read_initFile.o
grid_finalize.o : grid_data.o
block_finalize.o : block_data.o

hll.o		: grid_data.o primconsflux.o
hllc.o          : grid_data.o primconsflux.o
roe.o		: grid_data.o primconsflux.o eigensystem.o

io.o		: grid_data.o sim_data.o block_data.o
bc.o            : grid_data.o sim_data.o block_data.o DMR.o
bc_corners.o    : block_data.o sim_data.o grid_data.o
DMR.o           : grid_data.o sim_data.o


primconsflux.o  : grid_data.o eos.o

block_init.o    : block_data.o sim_data.o grid_data.o bc.o

sim_init.o	: sim_data.o read_initFile.o
sim_initBlock.o : sim_data.o grid_data.o primconsflux.o block_data.o bc.o

GP.o            : gp_data.o
gp_WENOinit.o   : gp_data.o linalg.o GP.o
gp_MakeStencil.o    : gp_data.o
gp_MDinit.o     : linalg.o GP.o gp_data.o grid_data.o
gp_Fluxinit.o   : linalg.o GP.o gp_data.o grid_data.o
gp_eigens.o     : gp_data.o
gpM_eigens.o    : gp_data.o
WENO.o          : gp_data.o
soln_gpWENO.o   : gp_data.o grid_data.o eigensystem.o primconsflux.o WENO.o sim_data.o

soln_gpFlux.o   : gp_data.o sim_data.o

gr_GPinit.o             : grid_data.o
sim_GPinit.o 		: grid_data.o linalg.o GP.o

soln_update.o		: grid_data.o primconsflux.o
soln_ReconEvolveAvg.o 	: grid_data.o sim_data.o primconsflux.o bc.o
soln_RK2.o              : grid_data.o primconsflux.o
soln_RK3.o              : grid_data.o primconsflux.o 
soln_RK4.o              : grid_data.o primconsflux.o

reconstruction.o        : grid_data.o sim_interfaces.o
soln_reconstruct.o 	: grid_data.o sim_data.o eigensystem.o primconsflux.o gp_data.o sim_interfaces.o reconstruction.o
soln_getFlux.o          : grid_data.o sim_data.o primconsflux.o sim_interfaces.o
soln_FOG.o		: grid_data.o eigensystem.o primconsflux.o 
soln_PLM.o		: grid_data.o sim_data.o slopeLimiter.o eigensystem.o primconsflux.o
soln_PPM.o              : grid_data.o sim_data.o slopeLimiter.o eigensystem.o primconsflux.o
soln_WENO.o             : grid_data.o sim_data.o eigensystem.o primconsflux.o WENO.o
soln_GP.o		: grid_data.o sim_data.o eigensystem.o primconsflux.o

#######################################################################################
