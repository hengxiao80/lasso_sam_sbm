# Makefile for various platforms
# Execute using Build csh-script only!
# Used together with Perl scripts in SRC/SCRIPT 
# (C) 2005 Marat Khairoutdinov
# $Id: Makefile 1595 2014-10-30 21:36:30Z dschanen@uwm.edu $
#------------------------------------------------------------------
# uncomment to disable timers:
#
#NOTIMERS=-DDISABLE_TIMERS
#-----------------------------------------------------------------

SAM = SAM_$(ADV_DIR)_$(SGS_DIR)_$(RAD_DIR)_$(MICRO_DIR)

### Platform independent preprocessor flags (UW Milwaukee ) ###
# SAM without CLUBB (revert to the original code)
CPPFLAGS = -DNETCDF

# REAL_TYPE is here simply so the user doesn't have to remove the CLUBB external directory
#CPPFLAGS = -DCLUBB_REAL_TYPE=8 

# SAM-CLUBB w/ modified radiation and microphysics restart.
#CPPFLAGS = -DCLUBB -DALTOCU_RAD -DMICRO_RESTART -DCLUBB_REAL_TYPE=8 -DUWM_STATS -DUWM_MISC

# SAM-CLUBB with UWM modifications
#CPPFLAGS = -DCLUBB -DCLUBB_REAL_TYPE=8 -DNETCDF -DUWM_STATS -DUWM_MISC

# SAM-CLUBB with UWM modifications and PNNL stats
#CPPFLAGS = -DCLUBB -DCLUBB_REAL_TYPE=8 -DNETCDF -DUWM_STATS -DUWM_MISC -DPNNL_STATS

### End platform independent preprocessor flags ###

### Platform independent files we want warnings on (UW Milwaukee code) ###
WARNABLE = clubb_sgs.o clubbvars.o clubb_silhs_vars.o stat_clubb.o sgs.o
### End platform independent warnings ###

# Determine platform 
PLATFORM := $(shell uname -s)

#------------------------------------------------------------------------
# Cray (on hopper), added by Minghuai Wang (Minghuai.Wang@pnnl.gov)
#
PLATFORM := Cray
ifeq ($(PLATFORM),Cray)
FF77 = mpif77 -c -extend-source 132
FF90 = mpif90 -c
CC = mpicc -c -DLINUX

LD = mpif90

INC_NETCDF   := $(NETCDF)/include
LIB_NETCDF   := $(NETCDF)/lib
INC_MPI      := /sw/cirrus/spack-envs/base/opt/linux-rhel8-x86_64/intel-20.0.4/openmpi-4.1.4-64cfg5aklxpi6qkfrwxfreamixupiexa/include
LIB_MPI      := /sw/cirrus/spack-envs/base/opt/linux-rhel8-x86_64/intel-20.0.4/openmpi-4.1.4-64cfg5aklxpi6qkfrwxfreamixupiexa/lib

##pgi compiler
#FFLAGS =  -DNO_R16 -i4 -Mdalign -Mextend -byteswapio
#FFLAGS =  -DNO_R16 -i4 -Mdalign -Mextend

## flags for debug
#FFLAGS += -g -Ktrap=fp -Mbounds -Kieee

##intel compiler
# FFLAGS = -O0 -g -traceback -check bounds -check uninit -ftrapuv
FFLAGS = -O3

##FFLAGS += -r8

# for SAM_CLUBB (noCLUBB)
FFLAGS += -I${INC_MPI} -I${INC_NETCDF} ${CPPFLAGS}
# for single moment with CLUBB
#FFLAGS += -I${INC_MPI} -I${INC_NETCDF} -DCLUBB -DCLUBB_REAL_TYPE=4 -DNETCDF -DCLUBB_SAM -Dsam1mom -DUWM_STATS -DUWM_MISC
# for Morrison microphysics with CLUBB
#FFLAGS += -I${INC_MPI} -I${INC_NETCDF} -DCLUBB -DCLUBB_REAL_TYPE=8 -DNETCDF -DCLUBB_SAM -DUWM_STATS -DUWM_MISC -DPNNL_STATS

#FFLAGS += -I${INC_MPI} -I${INC_NETCDF} -DNETCDF 
LDFLAGS = -L${LIB_MPI} -L${LIB_NETCDF} -lnetcdf -lnetcdff


###+++mhwang add for clubb (acml library on hopper)
#  ACML = $(ACML_DIR)/pgi64_mp/lib
#  ACML = $(ACML_DIR)/ifort64_mp/lib

#  LDFLAGS += -L${ACML} -Wl,-rpath,${ACML} -lacml_mp -lacml_mv
#  LDFLAGS += -L${ACML} -Wl,-rpath,${ACML} -lacml_mp
LDFLAGS += -mkl ## mwswong: intel compiler
#CFLAGS += -mkl
###---mhwang
endif


#----------------------------------
# AIX (tested only on IBM SP)
#

ifeq ($(PLATFORM),AIX)

#INC_MPI      := /usr/local/include
#LIB_MPI      := /usr/local/lib
#INC_NETCDF   := /usr/local/include
#LIB_NETCDF   := /usr/local/lib


FF77 = mpxlf90_r -c -qsuffix=f=f -qfixed=132
FF90 = mpxlf90_r -c -qsuffix=f=f90
CC = cc -c -DAIX
FFLAGS = -c -O3 -qstrict -qmaxmem=-1 -qarch=auto -qspillsize=5000 -Q
#FFLAGS = -c -qinitauto=FF -g -qflttrap=zerodivide:enable -qflttrap=ov:zero:inv:en
LD = mpxlf90_r
LDFLAGS = -bmaxdata:512000000 -bmaxstack:256000000

endif

#------------------------------------------------------------------------
# SGI
#------------------------------------------------------------------------

ifeq ($(PLATFORM),IRIX64)

INC_MPI      := /usr/local/include
LIB_MPI      := /usr/local/lib
INC_NETCDF   := /usr/local/include
LIB_NETCDF   := /usr/local/lib

FF77 = f90 -c -fixedform  -extend_source
FF90 = f90 -c -freeform
CC = cc -c -DIRIX64
FFLAGS = -O3 
#FFLAGS = -g -DEBUG:subscript_check=ON:trap_uninitialized=ON 
FFLAGS += -I${INC_MPI} -I${INC_NETCDF} ${CPPFLAGS}
LD = f90 
LDFLAGS = -L${LIB_MPI} -L${LIB_NETCDF} -lmpi -lnetcdf

endif
#----------------------------------------------------------------------
# Linux, Intel Compiler
#

#ifeq ($(PLATFORM),Linux)
#
#LIB_MPI = /usr/local/pkg/iopenmpi/lib
#INC_MPI = /usr/local/pkg/iopenmpi/include
#INC_NETCDF = /nfs/user08/marat/local/include
#LIB_NETCDF = /nfs/user08/marat/local/lib
#
#
#FF77 = /usr/local/pkg/iopenmpi/bin/mpif90 -c -fixed -extend_source
#FF90 = /usr/local/pkg/iopenmpi/bin/mpif90 -c
#CC = mpicc -c -DLINUX
#
#
#FFLAGS = -O3 
##FFLAGS = -g -ftrapuv -check all
#
#FFLAGS += -I${INC_MPI} -I${INC_NETCDF}
#LD = /usr/local/pkg/iopenmpi/bin/mpif90
#LDFLAGS = -L${LIB_NETCDF} -lnetcdf
#
#endif
#----------------------------------------------------------------------
# Linux, XLF compiler, Bluegene at San Diego SC
#
#ifeq ($(PLATFORM),Linux)

#INC_NETCDF   := /usr/local/apps/V1R3/netcdf-3.6.0-p1/include
#LIB_NETCDF   := /usr/local/apps/V1R3/netcdf-3.6.0-p1/lib

#FF77 = mpxlf90  -qarch=440 -qsuffix=f=f -qfixed=132
#FF90 = mpxlf90  -qarch=440 -qsuffix=f=f90
#CC = mpcc -c -DLinux
#FFLAGS = -c -O3 -qtune=440 -qstrict -qmaxmem=-1 -qspillsize=5000 -Q
##FFLAGS = -c -qinitauto=FF -g -qflttrap=zerodivide:enable -qflttrap=ov:zero:inv:en
#FFLAGS +=  -I${INC_NETCDF}
#LD = mpxlf90
#LDFLAGS = -L${LIB_NETCDF} -lnetcdf

#endif

#----------------------------------------------------------------------
# Linux, XLF compiler, Bluegene at BNL
#

#ifeq ($(PLATFORM),Linux)
ifeq ($(PLATFORM),NOTNOT) # Change to line above to use

ifeq ($(HOST),fen)
# Compute nodes:
INC_NETCDF   := /gpfs/home2/marat/local/netcdf/include
LIB_NETCDF   := /gpfs/home2/marat/local/netcdf/lib
FF77 = mpixlf90 -qarch=440 -qsuffix=f=f -qfixed=132
FF90 = mpixlf90 -qarch=440 -qsuffix=f=f90
CC = mpicc -c -DLinux
FFLAGS = -c -O3 -g -qtune=440 -qstrict -qmaxmem=-1 -qspillsize=5000 -Q
#FFLAGS = -c -qinitauto=FF -g -qflttrap=zerodivide:enable -qflttrap=ov:zero:inv:en
FFLAGS +=  -I${INC_NETCDF} ${CPPFLAGS}
LD = mpixlf90  
LDFLAGS = -L${LIB_NETCDF} -L/bgl/local/lib/ -L/opt/ibmmath/lib -lnetcdf -llapack.rts -lesslbg

else
# front end nodes and vis nodes
INC_NETCDF   := /gpfs/home2/marat/netcdf-vis-3.6.2/include
LIB_NETCDF   := /gpfs/home2/marat/netcdf-vis-3.6.2/lib
LIB_MPI = /usr/local/mpich2-1.0.5p4/lib
INC_MPI = /usr/local/mpich2-1.0.5p4/include
FF77 = xlf90 -qsuffix=f=f -qfixed=132
FF90 = xlf90 -qsuffix=f=f90
CC = cc -c -DLinux
FFLAGS = -c -O3  -qstrict -qmaxmem=-1 -qspillsize=5000 -Q
#FFLAGS = -c -qinitauto=FF -g -qflttrap=zerodivide:enable -qflttrap=ov:zero:inv:en
FFLAGS +=  -I${INC_NETCDF} -I${INC_MPI} ${CPPFLAGS}
LD = xlf90 
LDFLAGS = -L${LIB_NETCDF} -L${LIB_MPI} -L/gpfs/home1/slatest/blas/BLAS  -L/gpfs/home1/slatest/vislapack-2/lapack-3.2 -lmpich -lnetcdf -llapack_RS6K -lblas_LINUX 
endif


endif
#----------------------------------------------------------------------
# Linux, Portland Group Compiler
#
# This is the standard compiler for our nightly plots.  The GNU Fortran compiler
# is no longer used for these because some cases showed strange results when
# running in a high resolution LES mode. This may be a RHEL/CentOS specific issue.

# The next 3 lines must be left unaltered so that nightly tests will work
# {nightly pgi} DO NOT REMOVE THIS LINE
ifeq ($(PLATFORM),Linux)
#ifeq ($(PLATFORM),NOTNOT) # Change to line above to use

PGI_MPI = /opt/pgi/linux86-64/2013/mpi/mpich

LIB_MPI = $(PGI_MPI)/lib
INC_MPI = $(PGI_MPI)/include
LIB_NETCDF = /usr/local/netcdf-pgi/lib
INC_NETCDF = /usr/local/netcdf-pgi/include
#LIB_LAPACK = -llapack -lblas
LIB_LAPACK = -lacml

FF77 = $(PGI_MPI)/bin/mpif90 -c -Mfixed
FF90 = $(PGI_MPI)/bin/mpif90 -c -Mfreeform
CC = gcc -c -DLINUX

#FFLAGS = -g -C -Ktrap=fp -Kieee -Mcache_align
FFLAGS = -g -C -Munroll -O2 -Kieee -Mcache_align -Mvect -Ktrap=fp -Mstandard

FFLAGS += -I${INC_MPI} -I${INC_NETCDF} ${CPPFLAGS}
LD = $(PGI_MPI)/bin/mpif90 $(FFLAGS)
LDFLAGS = -L${LIB_MPI} -lmpich -L${LIB_NETCDF} -lnetcdf
LDFLAGS += ${LIB_LAPACK}

endif
#----------------------------------------------------------------------
# Linux, Sun Studio Compiler / MPICH (UW Milwaukee)
#
# This configuration is less well tested than PGI or Intel on Linux

# The next 3 lines must be left unaltered so that nightly tests will work
# {nightly sunstudio} DO NOT REMOVE THIS LINE
#ifeq ($(PLATFORM),Linux)
ifeq ($(PLATFORM),NOTNOT) # Change to line above to use

SUN_MPI = /usr/local/mpi/mpich2-1.4.1p1-oracle

LIB_MPI = $(SUN_MPI)/lib
INC_MPI = $(SUN_MPI)/include
LIB_NETCDF = /usr/local/netcdf-sun64/lib
INC_NETCDF = /usr/local/netcdf-sun64/include

FF77 = $(SUN_MPI)/bin/mpif90 -c 
FF90 = $(SUN_MPI)/bin/mpif90 -c
CC = suncc -m64 -c -DSUNOS

# Optimizing
#FFLAGS = -m64 -dalign -g -xO3 -xarch=native -ftrap=%none -DNO_LAPACK_ISNAN
# Debugging
FFLAGS = -traceback=common -m64 -dalign -g -fns=no -stackvar -xcheck=init_local -ftrap=common -DNO_LAPACK_ISNAN
# Regular SAM
#FFLAGS += -I${INC_MPI} -I${INC_NETCDF}
# For SAM-CLUBB
FFLAGS += -I${INC_MPI} -I${INC_NETCDF} ${CPPFLAGS} 
LD = $(SUN_MPI)/bin/mpif90 -g -m64 -dalign -ftrap=common
LDFLAGS = -L${LIB_MPI} -R$(LIB_MPI) -L${LIB_NETCDF} -lnetcdf 
# Sun performance library
LDFLAGS += -xlic_lib=sunperf

#WARNINGS = -w1 -ansi # The w1 flag enables warnings and errors, but not cautions

#Clubb standalone
WARNINGS = -w3 -ansi # This is the same as used in CLUBB standalone

#WARNINGS = -w4 -ansi # Comments, notes and cautions.

endif

#----------------------------------------------------------------------
# Linux, G95 Compiler / MPICH (UW Milwaukee)
#
# Note: This compiler is not regularly updated any longer, so it should be used
# with caution.

# The next 3 lines must be left unaltered so that nightly tests will work
# {nightly g95} DO NOT REMOVE THIS LINE
#ifeq ($(PLATFORM),Linux)
ifeq ($(PLATFORM),NOTNOT) # Change to line above to use

G95_MPI = /usr/local/mpi/mpich2-1.1-g95-ssm

LIB_MPI = $(G95_MPI)/lib
INC_MPI = $(G95_MPI)/include
LIB_NETCDF = /usr/local/netcdf-g95/lib
INC_NETCDF = /usr/local/netcdf-g95/include

FF77 = $(G95_MPI)/bin/mpif90 -c 
FF90 = $(G95_MPI)/bin/mpif90 -c 
CC = gcc -c -DLINUX

# Optimizing with array bounds checking.
# Note that as of gcc4, sse is the default floating point unit for x64 and 
# all math will use the IEEE-754 floating point standard
FFLAGS = -g -fbounds-check -O3 -msse3 

# Regular SAM
#LD = g95
#FFLAGS += -I${INC_MPI} -I${INC_NETCDF}
#LDFLAGS = -L${LIB_MPI} -lmpi -lmpi_f77 -L${LIB_NETCDF} -lnetcdf
# For SAM_CLUBB
FFLAGS += -I${INC_MPI} -I${INC_NETCDF} ${CPPFLAGS}
LD = $(G95_MPI)/bin/mpif90
LDFLAGS = -L${LIB_MPI} -L${LIB_NETCDF} -lnetcdf ${LIB_LAPACK}
# ATLAS BLAS / LAPACK
LDFLAGS += -L/usr/local/lib -llapack -L/usr/local/atlas/lib -lf77blas -lcblas -latlas 

WARNINGS = -Wall -Wextra -Wno=142,165,167 -pedantic

endif
#----------------------------------------------------------------------
# Linux, GNU Fortran Compiler / MPICH (UW Milwaukee)
#
# The standard gfortran that comes with RHEL/CentOS has shown some strange
# results in high resolution LES mode, so this is no longer our standard 
# compiler at UW Milwaukee.

# The next 3 lines must be left unaltered so that nightly tests will work
# {nightly gfortran} DO NOT REMOVE THIS LINE
#ifeq ($(PLATFORM),Linux)
ifeq ($(PLATFORM),NOTNOT) # Change to line above to use

MPI = /usr/local/mpi/mpich2-1.4.1p1-gfortran

LIB_MPI = $(MPI)/lib
INC_MPI = $(MPI)/include
LIB_NETCDF = /usr/local/netcdf-gfortran/lib
INC_NETCDF = /usr/local/netcdf-gfortran/include

FF77 = $(MPI)/bin/mpif90 -c 
FF90 = $(MPI)/bin/mpif90 -c 
CC = gcc -c -DLINUX

# Optimizing with array bounds checking.
# Note that as of gcc4, sse is the default floating point unit for x64 and 
# all math will use the IEEE-754 floating point standard
#FFLAGS = -march=native -g -fbounds-check -O2 -mieee-fp

# Debugging
FFLAGS = -g -fbounds-check -mieee-fp -finit-real=nan -finit-integer=-99999 -finit-logical=false -fall-intrinsics

# For CLEX9, et al. (restart Morrison microphysics, simplified radiation)
#FFLAGS += -I${INC_MPI} -I${INC_NETCDF} -DCLUBB -DALTOCU_RAD -DMICRO_RESTART -DCLUBB_REAL_TYPE=8 -DUWM_STATS -DUWM_MISC
# For standard SAM_CLUBB
FFLAGS += -I${INC_MPI} -I${INC_NETCDF} ${CPPFLAGS}

LD = $(MPI)/bin/mpif90 ${FFLAGS}
LDFLAGS = -L${LIB_MPI} -L${LIB_NETCDF} -lnetcdf ${LIB_LAPACK}
# AMD Core math library v4
ACML = /opt/acml4.4.0/gfortran64/lib
LDFLAGS += -L${ACML} -Wl,-rpath ${ACML} -lacml -lacml_mv

WARNINGS = -Wall -Wconversion -Wunderflow -Wcharacter-truncation -std=f95 -pedantic

endif

#----------------------------------------------------------------------
# Linux, Intel Fortran Compiler / MPICH (UW Milwaukee)
#
# This configuration is tested to some degree at UW-Milwaukee, and the Intel
# Fortran compiler has been used to put the CLUBB PBL into the WRF model.
# UWM uses Intel Fortran 13.0.1.117

# The next 3 lines must be left unaltered so that nightly tests will work
# {nightly intel} DO NOT REMOVE THIS LINE
#ifeq ($(PLATFORM),Linux)
ifeq ($(PLATFORM),NOTNOT) # Change to line above to use

INTEL_MPI = /usr/local/mpi/mpich2-1.4.1p1-intel

LIB_MPI = $(INTEL_MPI)/lib
INC_MPI = $(INTEL_MPI)/include 
LIB_NETCDF = /usr/local/netcdf-intel64/lib
INC_NETCDF = /usr/local/netcdf-intel64/include
LIB_LAPACK = -mkl=sequential

FF77 = $(INTEL_MPI)/bin/mpif90 -xHost -c 
FF90 = $(INTEL_MPI)/bin/mpif90 -xHost -c 
CC = $(INTEL_MPI)/bin/mpicc -c -DLINUX

# Optimizing
# -fp-model precise is probably slightly slower, but provides more accurate results
#FFLAGS = -O3 -ipo0 
# Debugging
FFLAGS = -debug full -CB -g -traceback -check bounds -check uninit -ftrapuv -fpe0 -ftz 

# Regular SAM
#FFLAGS += -I${INC_MPI} -I${INC_NETCDF}
# For SAM_CLUBB
FFLAGS += -I${INC_MPI} -I${INC_NETCDF} ${CPPFLAGS}

LD = $(INTEL_MPI)/bin/mpif90 $(FFLAGS)
#LD = ifort $(FFLAGS)
LDFLAGS = -L${LIB_MPI} -L${LIB_NETCDF} -lnetcdf ${LIB_LAPACK}

WARNINGS = -warn -warn notruncated_source
endif

#------------------------------------------------------------------------------
# Linux, Intel Fortran Compiler / hd1 (UWM CLUSTER)
#
# This configuration is used to run SAM_CLUBB on hd1

# Note, you must load the correct modules. This is done in the
# submission script (SCRIPTS/HD1). However,for clairty, they are: 

#ifeq ($(PLATFORM),Linux)
ifeq ($(PLATFORM),NOTNOT) # Change to line above to use

INTEL_MPI = /sharedapps/ICC/11.1/openmpi/1.6.5

LIB_MPI = $(INTEL_MPI)/lib
INC_MPI = $(INTEL_MPI)/include

LIB_NETCDF = /sharedapps/ICC/11.1/netcdf/4.3.0/lib
INC_NETCDF = /sharedapps/ICC/11.1/netcdf/4.3.0/include

LIB_ZLIB = /sharedapps/ICC/11.1/zlib/1.2.8/lib
INC_ZLIB = /sharedapps/ICC/11.1/zlib/1.2.8/include

LIB_HDF5 = /sharedapps/ICC/11.1/hdf5/1.8.12/lib
INC_HDF5 = /sharedapps/ICC/11.1/hdf5/1.8.12/include

LIB_LAPACK = -mkl=sequential


FF77 = ifort -xHost -c
FF90 = ifort -xHost -c
CC   = icc -c -DLINUX

# Optimizing
# -fp-model precise is probably slightly slower, but provides more accurate results
#FFLAGS = -O3 -ipo0 
# Debugging
FFLAGS = -debug full -CB -g -traceback -check bounds -check uninit -ftrapuv

# Regular SAM
#FFLAGS += -I${INC_MPI} -I${INC_NETCDF}
# For SAM_CLUBB
FFLAGS += -I${INC_MPI} -I${INC_NETCDF} -I${INC_ZLIB} -I${INC_HDF5} ${CPPFLAGS} -xHost -ip -no-prec-div -static-intel

# For SAM
#FFLAGS += -I${INC_MPI} -I${INC_NETCDF} -I${INC_ZLIB} -I${INC_HDF5} -DNETCDF -DUWM_STATS -DUWM_MISC -xHost -ip -no-prec-div -static-intel

LD = mpif90 $(FFLAGS)
LDFLAGS = -L${LIB_MPI} -L${LIB_NETCDF} -lnetcdff -lnetcdf ${LIB_LAPACK} -L${LIB_ZLIB} -L${LIB_HDF5}

WARNINGS = -warn -warn notruncated_source
endif

#----------------------------------------------------------------------
# Linux, Pathscale Fortran Compiler / Open MPI
# This configuration is largely untested at UW-Milwaukee, 
# since we don't have a license to use the Pathscale compiler.

#ifeq ($(PLATFORM),Linux)
ifeq ($(PLATFORM),NOTNOT) # Change to line above to use

PATHSCALE_MPI = /usr/local/openmpi-1.3.1-pathscale

LIB_MPI = $(PATHSCALE_MPI)/lib
INC_MPI = $(PATHSCALE_MPI)/include 
LIB_NETCDF = /usr/local/netcdf-pathscale/lib
INC_NETCDF = /usr/local/netcdf-pathscale/include
ACMLPATH    := /opt/acml4.2.0/pathscale64/lib
LIB_LAPACK     := -L$(ACMLPATH) -Wl,-rpath,$(ACMLPATH) -lacml -lacml_mv

FF77 = $(PATHSCALE_MPI)/bin/mpif77 -march=auto -c 
FF90 = $(PATHSCALE_MPI)/bin/mpif90 -march=auto -c 
CC = $(PATHSCALE_MPI)/bin/mpicc -march=auto -c -DLINUX

# Optimizing
#FFLAGS = -g -ffortran-bounds-check -O3 
FFLAGS = -g -ffortran-bounds-check

# Regular SAM
#FFLAGS += -I${INC_MPI} -I${INC_NETCDF}

# For SAM_CLUBB
FFLAGS += -I${INC_MPI} -I${INC_NETCDF} ${CPPFLAGS} 

LD = $(PATHSCALE_MPI)/bin/mpif90 -march=auto -O3
LDFLAGS = -L${LIB_MPI} -L${LIB_NETCDF} -lnetcdf ${LIB_LAPACK}

WARNINGS = -Wall -ansi 
endif

#----------------------------------------------------------------------
# Linux, Open64 Fortran Compiler / MPICH (UW Milwaukee steele.math.uwm.edu)
#
# Like Intel Fortran, this is moderately well tested with SAM-CLUBB
# While this compiler performs pretty well, the SE Linux policies may 
# require some adjustment in order to compile and run with Open64.

#ifeq ($(PLATFORM),Linux)
ifeq ($(PLATFORM),NOTNOT) # Change to line above to use

OPEN64_MPI = /usr/local/openmpi-1.4.1-open64

LIB_MPI = $(OPEN64_MPI)/lib
INC_MPI = $(OPEN64_MPI)/include 
LIB_NETCDF = /usr/local/netcdf-open64/lib
INC_NETCDF = /usr/local/netcdf-open64/include
ACMLPATH    := /opt/acml4.3.0/open64_64/lib
LIB_LAPACK     := -L$(ACMLPATH) -Wl,-rpath,$(ACMLPATH) -lacml -lacml_mv

FF77 = $(OPEN64_MPI)/bin/mpif90 -c 
FF90 = $(OPEN64_MPI)/bin/mpif90 -c 
CC = $(OPEN64_MPI)/bin/mpicc -c -DLINUX

# Optimizing
FFLAGS = -g -C -O2 -m3dnow -msse3 # Newer AMD processors
#FFLAGS = -g -C -O2 -msse3 # Newer Intel processors

FFLAGS += -I${INC_MPI} -I${INC_NETCDF} ${CPPFLAGS}

LD = $(OPEN64_MPI)/bin/mpif90 ${FFLAGS}
LDFLAGS = -L${LIB_MPI} -L${LIB_NETCDF} -lnetcdf ${LIB_LAPACK}

WARNINGS = -Wall -ansi 
endif


#--------------------------------------------
# Apple Mac OS X (Darwin) (Absoft Fortran)
#

#ifeq ($(PLATFORM),Darwin)

#INC_NETCDF   := /usr/local/absoft/include
#LIB_NETCDF   := /usr/local/absoft/lib
#INC_MPI      := /usr/local/absoft/include
#LIB_MPI       := /usr/local/absoft/lib

#FF77 = f90 -c -f fixed
#FF90 = f90 -c -f free
#CC = cc -c -DMACOSX

#FFLAGS = -O3 -noconsole -nowdir -YEXT_NAMES=LCS -s -YEXT_SFX=_  -z4
#LD = f90
#LDFLAGS = -L${LIB_MPI} -L${LIB_NETCDF} -lmpich -lnetcdf

#endif

#--------------------------------------------
# Apple Mac OS X (Darwin) (NAG Fortran)
#

#ifeq ($(PLATFORM),Darwin)

#INC_NETCDF   := /usr/local/nag/include
#LIB_NETCDF   := /usr/local/nag/lib
#INC_MPI      := /usr/local/nag/include
#LIB_MPI       := /usr/local/nag/lib

#FF77 = f95 -c -fixed -kind=byte
#FF90 = f95 -c -free -kind=byte
#CC = cc -c -DMACOSX

#FFLAGS =      # don't use any optimization -O* option! Will crash!
##FFLAGS = -gline -C=all -C=undefined  # use for debugging
#FFLAGS += -I$(SAM_SRC)/$(RAD_DIR) -I${INC_MPI} -I${INC_NETCDF}
#LD = f95 
##LDFLAGS = -L${LIB_MPI} -L${LIB_NETCDF} -lmpich -lnetcdf 
#LDFLAGS =  

#endif


#----------------------------------
# Apple Mac OS X (Darwin) (XLF compiler)
#

#ifeq ($(PLATFORM),Darwin)

#INC_NETCDF   := /usr/local/xlf/include
#LIB_NETCDF   := /usr/local/xlf/lib
#INC_MPI      := /usr/local/xlf/include
#LIB_MPI       := /usr/local/xlf/lib

#FF77 = xlf90 -c -qsuffix=f=f -qfixed=132
#FF90 = xlf90 -c -qsuffix=f=f90
#CC = cc -c -DMACOSX 
#FFLAGS = -c -O3 -qstrict -qmaxmem=-1 -qarch=auto -qspillsize=5000 -Q
##FFLAGS = -c -qinitauto=FF -g -qflttrap=zerodivide:enable -qflttrap=ov:zero:inv:en
#FFLAGS += -I$(SAM_SRC)/$(RAD_DIR) -I${INC_MPI} -I${INC_NETCDF}
#LD = xlf90
#LDFLAGS = -L${LIB_MPI} -L${LIB_NETCDF} -lmpi -lnetcdf

#endif

#----------------------------------
# Apple Mac OS X (Darwin) (Intel compiler)
#

ifeq ($(PLATFORM),Darwin)

INC_MPI      := /usr/local/include
LIB_MPI       := /usr/local/lib
INC_NETCDF      := /usr/local/include
LIB_NETCDF       := /usr/local/lib

FF77 = mpif90 -c -cm -nowarn -fixed -extend_source
FF90 = mpif90 -c -cm -nowarn -free
CC =  mpicc -c -DLINUX  

FFLAGS = -O3 -pad
#FFLAGS = -g -ftrapuv -check all

FFLAGS += -I${INC_MPI} -I${INC_NETCDF} ${CPPFLAGS}

LD = mpif90
LDFLAGS =  -L${LIB_NETCDF} -lnetcdf




endif

#----------------------------------
# Sun Solaris 10 (Sun Fortran 8.3 compiler / Cluster Tools 6)
#
# Options chosen for UW Milwaukee SunFire X4600 (8 processor SMP) machine,
# using Sun Cluster Tools 6 for MPI.
#  Note: this is no longer tested or maintained at UWM

ifeq ($(PLATFORM),SunOS)

INC_MPI      := /opt/SUNWhpc/include/amd64
LIB_MPI      := /opt/SUNWhpc/lib/amd64
INC_NETCDF   := /usr/local/netcdf-sun64/include
LIB_NETCDF   := /usr/local/netcdf-sun64/lib

FF77 = mpf90 -xtarget=native64 -g -c
FF90 = mpf90 -xtarget=native64 -g -c
CC = mpcc -xtarget=native64 -g -c -DSUNOS

# Optimizing
#FFLAGS = -xO3 -xvector=simd -dalign -ftrap=%none

# Debugging
FFLAGS = -ftrap=common -C -dalign -stackvar -xcheck=init_local

FFLAGS += -I${INC_MPI} -I${INC_NETCDF} ${CPPFLAGS}

# Optimizing
#LD = mpf90 -g -xO3 -xvector=simd -xarch=native64 -dalign -ftrap=%none

# Debugging
LD = mpf90 -g -xarch=native64 -dalign -ftrap=common -C

# Regular SAM
#LDFLAGS = -xlibmopt -L${LIB_NETCDF} -lnetcdf -L${LIB_MPI} -R${LIB_MPI} -lmpi_mt
# For SAM_CLUBB
LDFLAGS = -xlibmopt -L${LIB_MPI} -R${LIB_MPI} -xlic_lib=sunperf -L${LIB_NETCDF} -lnetcdf -lmpi_mt 

WARNINGS = -w3
endif

#----------------------------------
#----------------------------------------------
# you don't need to edit below this line


#compute the search path
dirs := . $(shell cat Filepath)
VPATH    := $(foreach dir,$(dirs),$(wildcard $(dir))) 

.SUFFIXES:
.SUFFIXES: .f .f90 .F90 .F .c .o



all: $(SAM_DIR)/$(SAM)
#ifneq (,$(findstring DCLUBB, $(FFLAGS) ) ) 	
#	@perl $(UTIL_CLUBB)/CLUBBStandardsCheck.pl $(SAM_SRC)/SGS_CLUBB/clubb_sgs.F90;
#	@perl $(UTIL_CLUBB)/CLUBBStandardsCheck.pl $(SAM_SRC)/SGS_CLUBB/clubb_silhs_vars.F90;
#	@perl $(UTIL_CLUBB)/CLUBBStandardsCheck.pl $(SAM_SRC)/SGS_CLUBB/clubbvars.F90;
#	@perl $(UTIL_CLUBB)/CLUBBStandardsCheck.pl $(SAM_SRC)/SGS_CLUBB/stat_clubb.F90;
#endif
SOURCES   := $(shell cat Srcfiles)

Depends: Srcfiles Filepath
	$(SAM_SRC)/SCRIPT/mkDepends Filepath Srcfiles > $@

Srcfiles: Filepath
	$(SAM_SRC)/SCRIPT/mkSrcfiles > $@

OBJS      := $(addsuffix .o, $(basename $(SOURCES))) 

$(SAM_DIR)/$(SAM): $(OBJS)
	$(LD) -o $@ $(OBJS) $(LDFLAGS)

.F90.o:
	${FF90}  ${FFLAGS} $(if $(filter $@, $(WARNABLE) ), $(WARNINGS) ) $<;
.f90.o:
	${FF90}  ${FFLAGS} $(if $(filter $@, $(WARNABLE) ), $(WARNINGS) ) $<;
.f.o:
	${FF77}  ${FFLAGS} $<
.F.o:
	${FF77}  ${FFLAGS} $<
.c.o:
	${CC}  ${CFLAGS} -I$(SAM_SRC)/TIMING $(NOTIMERS) $<

include Depends
