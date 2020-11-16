
MAKE = make
F90 = gfortran
F90_OPTS = -O3 -ffast-math -funroll-loops -fopenmp  -m64 -I/usr/include/mkl
F77 = gfortran
F77_OPTS = -O3 -ffast-math -funroll-loops -fopenmp  -m64 -I/usr/include/mkl
AR = ar
LIB_SYS =  -L/usr/lib/x86_64-linux-gnu -Wl,--no-as-needed -lmkl_rt -lpthread -lm -ldl
# LAPACK and BLAS libraries
LIB_LPK = -llapack -lblas
LIB_FFT =
SRC_OMP =

#-------------------------------------------------------------------------------
#SRC_MPI = mpi_stub.f90
# To enable MPI parallelism the MPI version of the Fortran compiler must be
# used. This is usually mpif90. Uncomment the following lines and run
# 'make clean' followed by 'make'.
F90 = mpif90
F77 = mpif90
SRC_MPI =
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#SRC_MKL = mkl_stub.f90
# To enable MKL parallelism, link with the MKL library then uncomment the
# following line and run 'make clean' followed by 'make'.
SRC_MKL =
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
SRC_OBLAS = oblas_stub.f90
# To enable OpenBLAS parallelism, link with the OpenBLAS library then uncomment
# the following line and run 'make clean' followed by 'make'.
#SRC_OBLAS =
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
SRC_BLIS = blis_stub.f90
# To enable BLIS parallelism, link with the BLIS library then uncomment the
# following line and run 'make clean' followed by 'make'.
#SRC_BLIS =
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#SRC_libxc = libxcifc_stub.f90
# To enable Libxc first download and compile version 5.x of the library. Next
# copy the files libxcf90.a and libxc.a to the elk/src directory. Then uncomment
# the following lines and run 'make clean' followed by 'make'.
LIB_libxc = libxcf90.a libxc.a
SRC_libxc = libxcf90.f90 libxcifc.f90
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#SRC_FFT = zfftifc.f90
# To use a different FFT library, copy the relevant library or include files to
# the elk/src directory (eg. mkl_dfti.f90), uncomment the appropriate lines
# below and run 'make clean' followed by 'make'.
#SRC_FFT = zfftifc_fftw.f90
#LIB_FFT = libfftw3.a
SRC_FFT = mkl_dfti.f90 zfftifc_mkl.f90
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#SRC_W90S = w90_stub.f90
# To enable the Wannier90 library copy libwannier.a to the elk/src directory.
# Then uncomment the following lines and run 'make clean' followed by 'make'.
SRC_W90S =
LIB_W90 = libwannier.a
#-------------------------------------------------------------------------------
