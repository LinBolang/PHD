#LIBS = -lblas -llapack -larpack
HOST = $(shell hostname)
ifeq ($(NERSC_HOST),edison)
LIBS = $(ARPACK) -lcuba
FC = ftn
else
ifeq ($(HOST),ubuntu)
LIBS = $(shell pkg-config --libs arpack) #-lcuba
FC      = ifort
FFLAGS	= -O3 -g -CB
else
LIBS = $(shell pkg-config --libs arpack) -lblas -llapack -larpack -lm #-lcuba
LIBS_PATH = -L/Users/xusq/Documents/BLFQ/proton/Cuba-4.2
FC      = gfortran
FFLAGS	= -O3 -g #-fcheck=all -Wall
endif
endif

#FC=/opt/intel/composer_xe_2013.3.163/bin/intel64/ifort
#FFLAGS=  -xAVX  -g -CB

OBJS=proton_with_effective_potential.o numbers.o wfproduction.o basis_info.o basis_eval.o basis_enum.o function.o olaps.o TMCv2.o hamiltonian.o integral.o hamiMatrix.o diag.o colorfactor.o distribution_functions_single.o output.o basis_enum_renorm.o renormalization.o diag_single.o hamitonian_renormalization.o
OBJS2=integration.o 
MODS=numbers.mod basis_info.mod hamiMatrix.mod colorfactor.mod

%.mod: %.f90
	$(FC) $(FFLAGS) -c  $<
%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<
pos.out:  $(MODS) $(OBJS)
	$(FC) $(FFLAGS) $(LIBS_PATH) -o $@ $(OBJS) $(LIBS)
	rm *.o
	rm *.mod
integration.out: $(MODS) $(OBJS2)
	$(FC) $(FFLAGS) $(LIBS_PATH) -o $@ $(OBJS2) $(LIBS)
	rm *.o
	rm *.mod
