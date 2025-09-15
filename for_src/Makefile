
include ../make.inc

all_files = config.o main_module.o timing_module.o main.o integrate.o \
            diagnose.o congrad.o parallel_mpi.o restart.o advection.o \
            tendency.o particles.o


all: $(all_files)
	$(F90)  $(all_files) $(F90FLAGS) $(CDFFLAGS) $(MPIFLAGS)  -o model.x

clean :
	 rm -f *.o *.mod 

timing_module.o    : timing_module.f90
	 	$(F90) $(F90FLAGS) -c timing_module.f90
        
main_module.o    : main_module.f90
	 	$(F90) $(F90FLAGS) $(MPIFLAGS) -c main_module.f90
main.o         : main.f90 main_module.o timing_module.o
	 	$(F90) $(F90FLAGS) $(MPIFLAGS)  -c main.f90
integrate.o: integrate.f90 main_module.o timing_module.o
		$(F90) $(F90FLAGS) $(MPIFLAGS) -c integrate.f90
tendency.o: tendency.f90 main_module.o timing_module.o
		$(F90) $(F90FLAGS) $(MPIFLAGS) -c tendency.f90        
advection.o: advection.f90 main_module.o timing_module.o
		$(F90) $(F90FLAGS) $(MPIFLAGS) -c advection.f90        
parallel_mpi.o: parallel_mpi.f90 main_module.o 
		$(F90) $(F90FLAGS) $(MPIFLAGS) -c parallel_mpi.f90		
restart.o: restart.f90 main_module.o 
		$(F90) $(F90FLAGS) $(MPIFLAGS) -c restart.f90		
congrad.o: congrad.f90 main_module.o 
		$(F90) $(F90FLAGS) $(MPIFLAGS) -c congrad.f90
diagnose.o: diagnose.F90 main_module.o
		$(F90) $(F90FLAGS) $(CDFFLAGS) $(MPIFLAGS) -c diagnose.F90 
particles.o: particles.f90 main_module.o
		$(F90) $(F90FLAGS) $(CDFFLAGS) $(MPIFLAGS) -c particles.f90         
config.o: config.f90 main_module.o particles.o
		$(F90) $(F90FLAGS) $(MPIFLAGS)  $(CDFFLAGS) -c config.f90

