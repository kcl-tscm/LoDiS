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
        force_mgo.f90\
        mgo_read_mot.f90 \
        force_mgo_mot.f90 \
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
        therma.f90 gauss.f90 vel_prep.f90 vel.f90 dft_test.f90\
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
        force_mgo.o\
        mgo_read_mot.o \
        force_mgo_mot.o \
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
        therma.o gauss.o vel_prep.o vel.o dft_test.o\
        time.o thermodynamics.o growth.o \
        output_calor.o output_quench.o output_growth.o output_coalescence.o \
        output_meta.o subnomi.o \

##================================================================================
# Write your Fortran compiler and its libraries here including -fopenmp. 
#
##================================================================================

#For example, uncomment (delete the #) the line below:
#F90=gfortran -fopenmp
#LIBS=


all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS)  $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c  $<

