OBJ = obj
SRC = src
MOD = mod
EXEC = xsMC+rf.${MACHINE}

ifeq ($(MACHINE),dlghp)

	HOME = /home/dg6/code

	F77	= ${HOME}/openmpi/gnu_64/bin/mpif90 -J${MOD}/ -g -fbacktrace #-march=core2 -O3  -fdefault-real-8
	F90	= ${HOME}/openmpi/gnu_64/bin/mpif90 -J${MOD}/ -g -fbounds-check -fbacktrace #-march=core2 -O3  -fdefault-real-8
	WARN = -Wall 
	DISLIN = -I ${HOME}/dislin/dislin64/gf -L ${HOME}/dislin/dislin64 -ldislin -DUSE_DISLIN
	NETCDF = -I ${HOME}/netcdf/netcdf_gnu64/include -L ${HOME}/netcdf/netcdf_gnu64/lib -lnetcdf
	PNETCDF = -I $(HOME)/pNetCdf/pnetcdf_gnu64/include $(HOME)/pNetCdf/pnetcdf_gnu64/lib/libpnetcdf.a

	#F77 = ${HOME}/openmpi/pgi/bin/mpif77 -g -Mbounds -traceback -Kieee
	#F90 = ${HOME}/openmpi/pgi/bin/mpif90 -g -Mbounds -traceback -module ${MOD}/ -Kieee
	#WARN = 
	#DISLIN = 
	#NETCDF = -I ${HOME}/netcdf/netcdf_pgi/include -L ${HOME}/netcdf/netcdf_pgi/lib -lnetcdf
	#PNETCDF = -I $(HOME)/pNetCdf/pnetcdf_pgi/include $(HOME)/pNetCdf/pnetcdf_pgi/lib/libpnetcdf.a

endif
ifeq ($(MACHINE),franklin)

	#F77 = f77 -g	
	#F90 = ftn -g -J${MOD}/
	#WARN = -fbounds-check -Wall 
	F77 = f77 -g -Kieee -traceback -Mbounds -fast
	F90 = ftn -g -module ${MOD}/ -Kieee -traceback -Mbounds -fast
	WARN = 
	#F77 = f77  
	#F90 = ftn -module ${MOD}/ 
	#WARN = -g -C 

endif
ifeq ($(MACHINE),jaguar)
	#F90 = ftn -g -J${MOD}/
	#WARN = -fbounds-check -Wall -fbacktrace
	F77 = f77 -g -Mbounds -traceback -Kieee 
	F90 = ftn -g -Mbounds -traceback -module ${MOD}/ -Kieee
	WARN = 
	NETCDF = ${NETCDF_FLIB}
	PNETCDF = ${PNETCDF_LIB} -L ${PNETCDF_DIR}
endif

OBJECTS = ${OBJ}/eqdsk.o \
${OBJ}/dlg.o \
${OBJ}/fitpack.o \
${OBJ}/gc_terms.o \
${OBJ}/interp.o \
${OBJ}/gc_integrate.o \
${OBJ}/read_particle_list.o \
${OBJ}/init_mpi.o \
${OBJ}/write_f_rzvv.o \
${OBJ}/read_namelist.o \
${OBJ}/beselI.o \
${OBJ}/constants.o \
${OBJ}/ranlib.o \
${OBJ}/luxury.o \
${OBJ}/read_ql.o \
${OBJ}/communications.o \
${OBJ}/read_mchoi.o \
${OBJ}/netcdf_check.o \
${OBJ}/powerAbsGrid.o \
${OBJ}/erf_external.o \
${OBJ}/collision_frequencies.o \
${OBJ}/beselJ.o \
${OBJ}/airy_functions.o

ifeq (${MACHINE},dlghp)
sMC-rf: ${SRC}/sMC-rf.f90 ${OBJECTS}
	${F90} ${SRC}/sMC-rf.f90 -o ${EXEC} ${OBJECTS} ${WARN} ${NETCDF} ${DISLIN} ${PNETCDF}
else
sMC-rf: ${SRC}/sMC-rf.f90 ${OBJECTS}
	${F90} ${SRC}/sMC-rf.f90 -o ${EXEC} ${OBJECTS} ${WARN} ${NETCDF} ${PNETCDF}
endif

${MOD}/write_f_rzvv.mod: ${SRC}/write_f_rzvv.f90 ${OBJ}/write_f_rzvv.o
${OBJ}/write_f_rzvv.o: ${SRC}/write_f_rzvv.f90 ${MOD}/powerAbsGrid.mod
	${F90} -c ${SRC}/write_f_rzvv.f90 -o ${OBJ}/write_f_rzvv.o ${WARN} ${NETCDF}

${MOD}/gc_integrate.mod: ${SRC}/gc_integrate.F90 ${OBJ}/gc_integrate.o 

ifeq ($(MACHINE),dlghp)
${OBJ}/gc_integrate.o: ${SRC}/gc_integrate.F90 ${MOD}/interp.mod ${MOD}/gc_terms.mod ${MOD}/read_namelist.mod ${MOD}/init_mpi.mod ${MOD}/luxury.mod ${MOD}/communications.mod ${MOD}/read_mchoi.mod ${MOD}/powerAbsGrid.mod ${MOD}/erf_external.mod ${MOD}/collision_frequencies.mod ${MOD}/bessJ_mod.mod ${OBJ}/airy_functions.o
	${F90} -c ${SRC}/gc_integrate.F90 -o ${OBJ}/gc_integrate.o ${WARN} ${DISLIN} 
else
${OBJ}/gc_integrate.o: ${SRC}/gc_integrate.F90 ${MOD}/interp.mod ${MOD}/gc_terms.mod ${MOD}/read_namelist.mod ${MOD}/init_mpi.mod ${MOD}/luxury.mod ${MOD}/communications.mod ${MOD}/read_mchoi.mod ${MOD}/powerAbsGrid.mod ${MOD}/erf_external.mod ${MOD}/collision_frequencies.mod ${MOD}/bessJ_mod.mod ${OBJ}/airy_functions.o
	${F90} -c ${SRC}/gc_integrate.F90 -o ${OBJ}/gc_integrate.o ${WARN} 
endif

${MOD}/communications.mod: ${SRC}/communications.f90 ${OBJ}/communications.o 
${OBJ}/communications.o : ${SRC}/communications.f90 ${MOD}/init_mpi.mod ${MOD}/read_ql.mod
	${F90} -c ${SRC}/communications.f90 -o ${OBJ}/communications.o ${WARN} 

${MOD}/read_ql.mod: ${SRC}/read_ql.F90 ${OBJ}/read_ql.o ${MOD}/init_mpi.mod ${MOD}/dlg.mod
${OBJ}/read_ql.o : ${SRC}/read_ql.F90
	${F90} -c ${SRC}/read_ql.F90 -o ${OBJ}/read_ql.o ${WARN} ${PNETCDF} ${NETCDF} 

${MOD}/read_particle_list.mod: ${SRC}/read_particle_list.f90 ${OBJ}/read_particle_list.o ${MOD}/init_mpi.mod
${OBJ}/read_particle_list.o: ${SRC}/read_particle_list.f90 ${MOD}/read_namelist.mod ${MOD}/dlg.mod ${MOD}/powerAbsGrid.mod
	${F90} -c ${SRC}/read_particle_list.f90 -o ${OBJ}/read_particle_list.o ${WARN} ${NETCDF}

${MOD}/init_mpi.mod: ${SRC}/init_mpi.f90 ${OBJ}/init_mpi.o
${OBJ}/init_mpi.o: ${SRC}/init_mpi.f90 ${MOD}/read_namelist.mod
	${F90} -c ${SRC}/init_mpi.f90 -o ${OBJ}/init_mpi.o ${WARN}

${MOD}/interp.mod: ${SRC}/interp.f90 ${OBJ}/interp.o 
${OBJ}/interp.o: ${SRC}/interp.f90 ${MOD}/gc_terms.mod ${MOD}/eqdsk.mod
	${F90} -c ${SRC}/interp.f90 -o ${OBJ}/interp.o ${WARN}

${MOD}/gc_terms.mod: ${SRC}/gc_terms.f90 ${OBJ}/gc_terms.o  
${OBJ}/gc_terms.o: ${SRC}/gc_terms.f90 ${MOD}/eqdsk.mod ${MOD}/constants.mod
	${F90} -c ${SRC}/gc_terms.f90 -o ${OBJ}/gc_terms.o ${WARN}

${MOD}/powerAbsGrid.mod: ${SRC}/powerAbsGrid.f90 ${OBJ}/powerAbsGrid.o  
${OBJ}/powerAbsGrid.o: ${SRC}/powerAbsGrid.f90 ${MOD}/eqdsk.mod 
	${F90} -c ${SRC}/powerAbsGrid.f90 -o ${OBJ}/powerAbsGrid.o ${WARN} ${NETCDF}

${MOD}/eqdsk.mod: ${SRC}/eqdsk.f90 ${OBJ}/eqdsk.o  
${OBJ}/eqdsk.o: ${SRC}/eqdsk.f90 ${MOD}/dlg.mod 
	${F90} -c ${SRC}/eqdsk.f90 -o ${OBJ}/eqdsk.o ${WARN}

${MOD}/read_mchoi.mod: ${SRC}/read_mchoi.f90 ${OBJ}/read_mchoi.o 
${OBJ}/read_mchoi.o: ${SRC}/read_mchoi.f90 ${MOD}/read_namelist.mod ${MOD}/netcdf_check.mod
	${F90} -c ${SRC}/read_mchoi.f90 -o ${OBJ}/read_mchoi.o ${WARN} ${NETCDF}

${MOD}/constants.mod: ${SRC}/constants.f90 ${OBJ}/constants.o
${OBJ}/constants.o: ${SRC}/constants.f90 ${MOD}/read_namelist.mod
	${F90} -c ${SRC}/constants.f90 -o ${OBJ}/constants.o ${WARN}

${MOD}/read_namelist.mod: ${SRC}/read_namelist.f90 ${OBJ}/read_namelist.o
${OBJ}/read_namelist.o: ${SRC}/read_namelist.f90
	${F90} -c ${SRC}/read_namelist.f90 -o ${OBJ}/read_namelist.o ${WARN}

${MOD}/netcdf_check.mod: ${SRC}/netcdf_check.f90 ${OBJ}/netcdf_check.o
${OBJ}/netcdf_check.o: ${SRC}/netcdf_check.f90
	${F90} -c ${SRC}/netcdf_check.f90 -o ${OBJ}/netcdf_check.o ${WARN} ${NETCDF}

${MOD}/dlg.mod: ${SRC}/dlg.f90 ${OBJ}/dlg.o
${OBJ}/dlg.o: ${SRC}/dlg.f90
	${F90} -c ${SRC}/dlg.f90 -o ${OBJ}/dlg.o ${WARN} ${NETCDF}

${OBJ}/beselI.o: ${SRC}/bessel/beselI.f90
	${F90} -c ${SRC}/bessel/beselI.f90 -o ${OBJ}/beselI.o ${WARN}

${MOD}/bessJ_mod.mod: ${SRC}/bessel/bessJ_mod.f90 ${OBJ}/beselJ.o
${OBJ}/beselJ.o: ${SRC}/bessel/bessJ_mod.f90
	${F90} -c ${SRC}/bessel/bessJ_mod.f90 -o ${OBJ}/beselJ.o ${WARN}

${MOD}/collision_frequencies.mod: ${SRC}/collision_frequencies.f90 ${OBJ}/collision_frequencies.o 
${OBJ}/collision_frequencies.o: ${SRC}/collision_frequencies.f90 ${MOD}/erf_external.mod ${MOD}/constants.mod
	${F90} -c ${SRC}/collision_frequencies.f90 -o ${OBJ}/collision_frequencies.o ${WARN}

${MOD}/erf_external.mod: ${SRC}/erf.f90 ${OBJ}/erf_external.o
${OBJ}/erf_external.o: ${SRC}/erf.f90
	${F90} -c ${SRC}/erf.f90 -o ${OBJ}/erf_external.o ${WARN}

${OBJ}/fitpack.o: ${SRC}/fitpack.f
	${F77} -c ${SRC}/fitpack.f -o ${OBJ}/fitpack.o 

${OBJ}/ranlib.o: ${SRC}/ranlib.f
	${F77} -c ${SRC}/ranlib.f -o ${OBJ}/ranlib.o 

${MOD}/luxury.mod: ${SRC}/luxury.f90 ${OBJ}/luxury.o
${OBJ}/luxury.o: ${SRC}/luxury.f90
	${F90} -c ${SRC}/luxury.f90 -o ${OBJ}/luxury.o

${OBJ}/airy_functions.o: ${SRC}/airy/airy_functions.f90 ${SRC}/airy/airy_real ${SRC}/airy/airy_complex ${SRC}/airy/airy_head ${SRC}/airy/airy_parameters
	$(F90) -c -o ${OBJ}/airy_functions.o ${SRC}/airy/airy_functions.f90

clean:
	rm ${OBJ}/*.o ${MOD}/*.mod ${EXEC} sMC-rf.o
