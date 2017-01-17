#-*- mode: makefile;-*-  
#make =  normal compile
#make clean = Removes all objects and executable
#make opt=debug = Compiles with debug flags

CC=mpiicc
FC=mpiifort
LD=mpiifort

LIB = Source/kspace_npt_openmp.a -Wl,--start-group ${MKL_HOME}/lib/intel64/libmkl_intel_lp64.a ${MKL_HOME}/lib/intel64/libmkl_core.a ${MKL_HOME}/lib/intel64/libmkl_intel_thread.a ${MKL_HOME}/lib/intel64/libfftw2xc_single_intel.a  -Wl,--end-group -lpthread -lm -lstdc++

#INC = -I/home/atom/a/parks23/usr/bin/openblas/include -I/opt/tecplot/current/include 
#LIB = -L/home/atom/a/parks23/usr/bin/openblas/lib -lopenblas 

ifeq ($(opt),debug)
#Flags for "make opt=debug"
CFLAGS=-pg #$(GTKFLAGS)
#FFLAGS = $(INC) -JModules -fimplicit-none -pg -fbounds-check -ffpe-trap=invalid,zero,overflow -Wconversion
FFLAGS = $(INC) -g -mkl -openmp 
LDFLAGS =  -mkl -openmp $(SRCDIR)/kspace_npt_openmp.a
endif

ifeq ($(opt),pro)
#Flags for "make opt=debug"
CFLAGS=-pg #$(GTKFLAGS)
#FFLAGS = $(INC) -JModules -fimplicit-none -pg -fbounds-check -ffpe-trap=invalid,zero,overflow -Wconversion
FFLAGS = $(INC) -JModules -vec-report2
LDFLAGS =  -mkl -openmp $(SRCDIR)/kspace_npt_openmp.a
endif

ifeq ($(opt),serial)
#Flags for "make"
CFLAGS=-O3 -mkl -align array64byte -fp-model fast=2 -fimf-domain-exclusion=15 -march=native -assume buffered_io
FFLAGS=$(INC) -O3 -mkl -align array64byte -march=native -assume buffered_io
LDFLAGS=  -mkl  $(SRCDIR)/kspace_npt_openmp.a
endif

ifeq ($(opt),par)
#Flags for "make"
CFLAGS=-O3
FFLAGS=$(INC) -O3  -xHost -mkl -openmp -align array64byte -fp-model fast=2 -fimf-domain-exclusion=15 -march=native -assume buffered_io
LDFLAGS = -mkl -openmp $(SRCDIR)/kspace_npt_openmp.a
endif

SRCDIR=Source
OBJDIR=Objects

CSOURCES= #$(SRCDIR)/main.c


FSOURCES= $(SRCDIR)/global.f90  $(SRCDIR)/mod_makefolders.f90 $(SRCDIR)/mod_memory.f90 $(SRCDIR)/get_started.f90 $(SRCDIR)/mod_assign_types.f90 $(SRCDIR)/mod_readfile.f90 $(SRCDIR)/mod_average.f90 $(SRCDIR)/mod_nondimensionalize.f90 $(SRCDIR)/mod_Gaussian.f90 $(SRCDIR)/mod_build_ptta.f90 $(SRCDIR)/mod_minimg.f90 $(SRCDIR)/mod_LRC.f90 $(SRCDIR)/mod_sample.f90 $(SRCDIR)/mod_file.f90 $(SRCDIR)/mod_GroFile.f90 $(SRCDIR)/mod_datadump.f90 $(SRCDIR)/mod_print_results.f90 $(SRCDIR)/mod_buildbond_lists.f90 $(SRCDIR)/mod_pbc.f90 $(SRCDIR)/mod_matrix_multiply.f90  $(SRCDIR)/mod_init_posit.f90 $(SRCDIR)/mod_neighbor.f90 $(SRCDIR)/mod_sort.f90 $(SRCDIR)/mod_neighbor_cluster.f90 $(SRCDIR)/mod_ensemble.f90 $(SRCDIR)/mod_langevin.f90 $(SRCDIR)/mod_force_correction.f90  $(SRCDIR)/mod_force_bond.f90 $(SRCDIR)/mod_force_angle.f90 $(SRCDIR)/mod_force_dihed.f90 $(SRCDIR)/mod_force_dihed_amber.f90 $(SRCDIR)/mod_force_impro_amber.f90 $(SRCDIR)/mod_force_nonbond.f90 $(SRCDIR)/mod_shake.f90 $(SRCDIR)/mod_force.f90 $(SRCDIR)/mod_velScale.f90 $(SRCDIR)/mod_prune.f90  $(SRCDIR)/mod_qcprot.f90 $(SRCDIR)/mod_adjust.f90 $(SRCDIR)/mod_cluster.f90 $(SRCDIR)/mod_nucleation.f90 $(SRCDIR)/mod_anneal.f90  $(SRCDIR)/mod_initialize_template.f90 $(SRCDIR)/mod_initialize.f90 $(SRCDIR)/mod_restart.f90 $(SRCDIR)/mod_restart_datadump.f90 $(SRCDIR)/mod_integrate.f90 $(SRCDIR)/MD.f90

FOBJECTS= $(FSOURCES:$(SRCDIR)/%.f90=$(OBJDIR)/%.o)
COBJECTS= $(CSOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)

EXECUTABLE=new.out

all: $(CSOURCES) $(FSOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(FOBJECTS) $(COBJECTS)  
	$(FC) $(LDFLAGS) $(FOBJECTS) $(COBJECTS) -o $@ $(LIB)

$(FOBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(COBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm Modules/*.mod $(EXECUTABLE) *.modmic *.mod Objects/*.o $(FOBJECTS) $(COBJECTS)

delete:
	rm *.dat xeon.* *.gro *.xyz fort.* mymachinefile *.log beta300.*

folder:
	rm -r DATA 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60