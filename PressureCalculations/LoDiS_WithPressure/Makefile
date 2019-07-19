PROG = LODIS_all

SRCS =  module_parameters.f90 \
        module_cluster.f90 \
        module_potential.f90 \
        module_enforce.f90 \
        module_distance.f90 \
        module_meta.f90 \
        module_mgo.f90\
        module_env.f90\
        module_sticky.f90\
        module_function_var.f90\
        main_MD.f90 \
        read_input.f90 \
        mgo_read.f90 \
        mgo_coord_no.f90 \
        mgo_force.f90\
        geometry_calc.f90\
	substrate_force_choice.f90 \
	substrate_hex_force.f90 \
	substrate_dsq_force.f90 \
        mgo_read_mot.f90 \
        mgo_force_mot.f90 \
        mgo_coord_no_mot.f90 \
        force_env.f90\
        env_coord_no.f90\
        env_coord_no2.f90\
        bim_cn.f90 bim_cn_light.f90\
        d_com.f90 d_com_light.f90\
        sticky_force.f90 \
        choice_pot.f90 read_rgl.f90\
	acf_routine.f90\
	dft_routine.f90\
	suppress.f90\
        preparo.f90 \
        dlaran.f90 init.f90 bigvoi.f90 voisin.f90\
        force_CHOICE.f90 force_rgl.f90 \
	force_lj1.f90 force_gir.f90 force_par.f90 \
        collective.f90 coord_no.f90 coord_no_light.f90 force_meta.f90 pair_distances.f90\
        secneigh_no.f90 secneigh_no_light.f90 stackfault_no.f90 stackfault_no_light.f90\
        comneighfun_no.f90 comneighfun_no_light.f90\
        therma.f90 gauss.f90 vel_prep.f90 vel.f90 \
        time.f90 thermodynamics.f90 growth.f90 \
        output_calor.f90 output_quench.f90 output_growth.f90 output_coalescence.f90 \
        output_meta.f90 subnomi.f90 \

OBJS = 	module_parameters.o \
        module_cluster.o \
        module_potential.o \
        module_enforce.o \
        module_distance.o \
        module_meta.o \
        module_mgo.o \
        module_env.o\
        module_sticky.o \
        module_function_var.o\
        main_MD.o \
        read_input.o \
        mgo_read.o \
        mgo_coord_no.o \
        mgo_force.o\
        geometry_calc.o\
	substrate_force_choice.o \
	substrate_hex_force.o \
	substrate_dsq_force.o \
        mgo_read_mot.o \
        mgo_force_mot.o \
        mgo_coord_no_mot.o \
        force_env.o\
        env_coord_no.o\
        env_coord_no2.o\
        bim_cn.o bim_cn_light.o\
        d_com.o d_com_light.o\
        sticky_force.o \
        choice_pot.o \
        read_rgl.o \
	acf_routine.o \
	dft_routine.o \
	suppress.o \
        preparo.o \
        dlaran.o init.o bigvoi.o voisin.o \
        force_CHOICE.o force_rgl.o force_lj1.o \
        force_gir.o force_par.o \
        collective.o coord_no.o coord_no_light.o force_meta.o pair_distances.o\
        secneigh_no.o secneigh_no_light.o stackfault_no.o stackfault_no_light.o\
        comneighfun_no.o comneighfun_no_light.o\
        therma.o gauss.o vel_prep.o vel.o \
        time.o thermodynamics.o growth.o \
        output_calor.o output_quench.o output_growth.o output_coalescence.o \
        output_meta.o subnomi.o \

##________________________________________________________________________________
## ifort 
##________________________________________________________________________________
##---> LIBS on Henry:
#MKLROOT=/opt/apps/intel/composer_xe_2011_sp1.7.256/mkl/
# LIBS =  -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a \
#   $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm 
##====
##---> LIBS on Ruby:
# MKLROOT=/opt/intel/mkl/10.2.1.017/
# LIBS =  $(MKLROOT)/lib/em64t/libmkl_solver_lp64_sequential.a -Wl,--start-group $(MKLROOT)/lib/em64t/libmkl_intel_lp64.a \
#   $(MKLROOT)/lib/em64t/libmkl_sequential.a $(MKLROOT)/lib/em64t/libmkl_core.a -Wl,--end-group -lpthread
##====
# FFLAGS =
# F90 = ifort
 #F90FLAGS = -C -O0 -g -traceback                          # DEBUG #-ffpe-trap=zero
 #LDFLAGS  = -C -O0 -g -traceback                          # DEBUG #-ffpe-trap=zero

# F90FLAGS = -C -O2 -real-size 64 -double-size 64 -openmp  # PERFORMANCE - OpenMP
# LDFLAGS  = -C -O2 -openmp                                # PERFORMANCE - OpenMP

 #F90FLAGS = -C -O0 -g -real-size 64 -double-size 64 -openmp  # valgrind
 #LDFLAGS  = -C -O0 -g -openmp                                # valgrind

 #F90FLAGS = -C -O2 -real-size 64  -double-size 64 -g -pg  # PROFILING
 #LDFLAGS  = -C -O2 -g -pg                                 # PROFILING
#================================================================================

##________________________________________________________________________________
## gfortran 
##________________________________________________________________________________
##---> LIBS on ADA:
# MKLROOT=/opt/apps/intel/mkl11.3/mkl/
# LIBS =   -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_gf_lp64.a $(MKLROOT)/lib/intel64/libmkl_gnu_thread.a \
#   $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -ldl -lpthread -lm
##====
##---> LIBS on Ubuntu et al. 
 #LIBS = /lib64/ld-linux-x86-64.so.2      # for gfortran ubuntu 12.04
#LIBS = /usr/lib/lapack/liblapack.so.3  # for gfortran scientific linux 64
 #LIBS = /usr/lib/liblapack.so.3.0.3    # for gfortran scientific linux 32
 LIBS = -llapack # for my Mac OSX laptop 
##====
# FFLAGS =
F90 = gfortran   
# F90 = /opt/local/bin/gfortran-mp-6 # on my Mac OSX laptop or often under /usr/local/bin
 # F90 = /usr/local/bin/gfortran-mp-6 
 #F90FLAGS = -C -O2 -g -fdefault-real-8 -fdefault-double-8 -fopenmp # OpenMP PERF
 #LDFLAGS  = -C -O2 -g -fopenmp                                     # OpenMP PERF
  F90FLAGS = -C -O2 -g -fbounds-check -Wall -Wextra -Warray-temporaries -Wconversion -fdefault-real-8 -fdefault-double-8 -fbacktrace 
  LDFLAGS  = $(F90FLAGS) -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan  -fopenmp
#
# #F90FLAGS = -C -O0 -g -fdefault-real-8 -fdefault-double-8 -fopenmp # valgrind
# #LDFLAGS  = -C -O0 -g -fopenmp                                     # valgrind
##================================================================================

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS)  $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c  $<

