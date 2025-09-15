
OSX = "" 

#########################################################################
# platform dependent options
#########################################################################

ifdef OSX
 F90 =  mpif90
 CDFFLAGS = -I/usr/local/include -lnetcdff -lnetcdf #-Dwith_netcdf_parallel
 MPIFLAGS =
 F90FLAGS =  -Wall -O3 -fbounds-check -finit-real=nan -fallow-argument-mismatch
endif


ifdef LEVANTE
 cdf = nf-config
 F90 = mpifort
 CDFFLAGS = -I`$(cdf) --includedir`  `$(cdf) --flibs` 
 CDFFLAGS = ${CDFFLAGS} -Wl,-rpath,/sw/spack-levante/netcdf-fortran-4.5.3-k6xq5g/lib  
 CDFFLAGS = ${CDFFLAGS} #-Dwith_netcdf_parallel
 MPIFLAGS =
 F90FLAGS =
endif




